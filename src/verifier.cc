/***************************************************************************
 *            verifier.cc
 *
 *  Copyright 2011  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <string>

#include "verifier.h"

namespace Ariadne {

Verifier::Verifier(
		const HybridReachabilityAnalyser& outer_analyser,
		const HybridReachabilityAnalyser& lower_analyser) :
			_outer_analyser(outer_analyser.clone()),
			_lower_analyser(lower_analyser.clone()),
			_settings(new VerificationSettings()),
			_safety_coarse_outer_approximation(new OuterApproximationCache()),
			_dominating_coarse_outer_approximation(new OuterApproximationCache()),
			_dominated_coarse_outer_approximation(new OuterApproximationCache())
{
	_outer_analyser->verb_tab_prefix = verb_tab_prefix + verifier_max_verbosity_level_used;
	_lower_analyser->verb_tab_prefix = verb_tab_prefix + verifier_max_verbosity_level_used;
}

Verifier::Verifier(const HybridReachabilityAnalyser& analyser) :
			_outer_analyser(analyser.clone()),
			_lower_analyser(analyser.clone()),
			_settings(new VerificationSettings()),
			_safety_coarse_outer_approximation(new OuterApproximationCache()),
			_dominating_coarse_outer_approximation(new OuterApproximationCache()),
			_dominated_coarse_outer_approximation(new OuterApproximationCache())
{
	_outer_analyser->verb_tab_prefix = verb_tab_prefix + verifier_max_verbosity_level_used;
	_lower_analyser->verb_tab_prefix = verb_tab_prefix + verifier_max_verbosity_level_used;
}

Verifier::~Verifier()
{
}


tribool
Verifier::
safety(SafetyVerificationInput& verInput) const
{
	RealConstantSet constants;

	return _safety_nosplitting(verInput,constants);
}


tribool
Verifier::
_safety(
		SafetyVerificationInput& verInput,
		const RealConstant& constant) const
{
	HybridAutomaton& system = verInput.getSystem();

	Real originalParameterValue = system.accessible_constant_value(constant.name());

	system.substitute(constant);
	tribool result = safety(verInput);
	system.substitute(constant,originalParameterValue);

	return result;
}

tribool
Verifier::
_safety(
		SafetyVerificationInput& verInput,
		const RealConstant& constant,
		const Float& value) const
{
	const RealConstant modifiedParameter(constant.name(),Interval(value));

	return _safety(verInput, modifiedParameter);
}


tribool
Verifier::
_safety_nosplitting(
		SafetyVerificationInput& verInput,
		const RealConstantSet& constants) const
{
	ARIADNE_LOG(2,"\n");
	ARIADNE_LOG(2,"Iterative verification...\n");

	HybridAutomaton& system = verInput.getSystem();

	if (_settings->plot_results)
		_plot_dirpath_init(system.name());

	_resetAndChooseInitialSafetySettings(system,verInput.getDomain(),constants);

	_outer_analyser->settings().maximum_grid_depth = _outer_analyser->settings().lowest_maximum_grid_depth;
	_lower_analyser->settings().maximum_grid_depth = _lower_analyser->settings().lowest_maximum_grid_depth;

	while (_outer_analyser->settings().maximum_grid_depth
			<= _outer_analyser->settings().highest_maximum_grid_depth ||
			_lower_analyser->settings().maximum_grid_depth
			<= _lower_analyser->settings().highest_maximum_grid_depth) {

		ARIADNE_LOG(2, "Outer depth " << _outer_analyser->settings().maximum_grid_depth << ", " <<
				"lower depth " << _lower_analyser->settings().maximum_grid_depth << "\n");

		tribool result = _safety_once(system,verInput.getInitialSet(),verInput.getSafetyConstraint(),constants);

    	_outer_analyser->settings().maximum_grid_depth++;
    	_lower_analyser->settings().maximum_grid_depth++;

		if (!indeterminate(result))
			return result;
    }

	// Return indeterminate
	return indeterminate;
}


tribool
Verifier::
_safety_once(
		SystemType& system,
		const HybridImageSet& initial_set,
		const HybridConstraintSet& safety_constraint,
		const RealConstantSet& constants) const
{
		ARIADNE_LOG(3, "Verification...\n");

		if (_safety_proving_once(system,initial_set,safety_constraint,constants)) {
			ARIADNE_LOG(3, "Safe.\n");
			return true;
		}

		if (_safety_disproving_once(system,initial_set,safety_constraint,constants)) {
			ARIADNE_LOG(3, "Unsafe.\n");
			return false;
		}

		ARIADNE_LOG(3, "Indeterminate.\n");
		return indeterminate;
}


bool
Verifier::
_safety_proving_once(
		SystemType& system,
		const HybridImageSet& initial_set,
		const HybridConstraintSet& safety_constraint,
		const RealConstantSet& constants) const
{
	bool result;

	ARIADNE_LOG(4,"Proving...\n");
	if (!_grid_depth_is_within_bounds_under(UPPER_SEMANTICS)) {
		ARIADNE_LOG(4,"Not proved.\n");
		return false;
	}

	RealConstantSet original_constants = system.accessible_constants();

	system.substitute(constants,_settings->use_param_midpoints_for_proving);

	ARIADNE_LOG(4,"Setting parameters for this proving iteration...\n");

	_tuneIterativeStepSettings(system,_safety_coarse_outer_approximation->get(),_safety_reachability_restriction,UPPER_SEMANTICS);

	ARIADNE_LOG(4,"Performing outer reachability analysis...\n");

	try {

		HybridGridTreeSet reach = _outer_analyser->outer_chain_reach(
				system,initial_set,DIRECTION_FORWARD,_safety_reachability_restriction);

		_update_safety_cached_reachability_with(reach);

		result = definitely(covers(safety_constraint,reach));

		ARIADNE_LOG(5, "The reachable set " << (!result ? "does not satisfy":"satisfies") << " the safety constraint.\n");

		if (_settings->plot_results)
			_plot_reach(reach,"outer",_outer_analyser->settings().maximum_grid_depth);

		// We refine only if we have no result from the initial reach set and we have a restriction on the reachability
		if (!result && !_safety_reachability_restriction.empty() && _settings->enable_fb_refinement_for_proving)
			result = _forward_backward_refinement_check(system,initial_set,safety_constraint,_safety_reachability_restriction);

	} catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(5, "The outer reached region is partially out of the domain (skipped).\n");
		result = false;
	} catch (ReachUnsatisfiesConstraintException ex) {
		ARIADNE_LOG(5, "The outer reached region is not inside the safe region (skipped).\n");
		result = false;
	}

	system.substitute(original_constants);

	ARIADNE_LOG(4, (result ? "Proved.\n" : "Not proved.\n") );

	return result;
}


bool
Verifier::
_forward_backward_refinement_check(
		SystemType& system,
		const HybridImageSet& initial_set,
		const HybridConstraintSet& constraint_set,
		const HybridGridTreeSet& reachability) const
{
	bool result = false;

	ARIADNE_LOG(5, "Performing forward-backward refinement...\n");

	HybridGridTreeSet residual_reachability = reachability;

	int& maximum_grid_depth = _outer_analyser->settings().maximum_grid_depth;

	/* Given the initial set I, forward reachability F, and unsafe region U:

	(if the settings can be refined)
	1. Refine the accuracy settings;
	2. Evaluate the final unsafe set H = F ∩ U ; if H = ∅, return TRUE;
	3. If the accuracy is over the maximum, return FALSE;
	4. Obtain the backward reachability B under the constraint set F;
	5. If B = F, return FALSE;
	(if the settings can be refined)
	6. Refine the accuracy settings;
	7. Evaluate the new initial set I = B ∩ I; if I = ∅, return TRUE;
	8. If the accuracy is over the maximum, return FALSE;
	9. Starting from I, obtain the forward reachability F under the constraint set B;
	10. If F = B, return FALSE;

	11. Restart from 1.

	*/

	while (true) {

		if (maximum_grid_depth >= _outer_analyser->settings().highest_maximum_grid_depth) {
			result = false;
			break;
		}

		++maximum_grid_depth;
		_outer_analyser->forward_backward_refine_evolution_settings();

		HybridGridTreeSet backward_initial = _reachability_refinement_starting_set(system,initial_set,constraint_set,
				residual_reachability,DIRECTION_BACKWARD);

		if (backward_initial.empty()) {
			result = true;
			break;
		}

		ARIADNE_LOG(5,"Retrieving backward reachability...\n");

		residual_reachability = _outer_analyser->outer_chain_reach(system,backward_initial,
				DIRECTION_BACKWARD,residual_reachability);

		ARIADNE_LOG(6,"Residual reachability size: " << residual_reachability.size() << "\n");

		if (_settings->plot_results)
			_plot_reach(residual_reachability,"backward",maximum_grid_depth);

		HybridGridTreeSet forward_initial = _reachability_refinement_starting_set(system,initial_set,constraint_set,
				residual_reachability,DIRECTION_FORWARD);

		if (forward_initial.empty()) {
			result = true;
			break;
		}

		ARIADNE_LOG(5,"Retrieving forward reachability...\n");

		residual_reachability = _outer_analyser->outer_chain_reach(system,forward_initial,
				DIRECTION_FORWARD,residual_reachability);

		ARIADNE_LOG(6,"Residual reachability size: " << residual_reachability.size() << "\n");

		if (_settings->plot_results)
			_plot_reach(residual_reachability,"forward",maximum_grid_depth);
	}

	_safety_reachability_restriction = residual_reachability;

	ARIADNE_LOG(5, (result ? "Successful.\n" : "Failed.\n"));

	return result;
}


HybridGridTreeSet
Verifier::
_reachability_refinement_starting_set(
		SystemType& system,
		const HybridImageSet& initial_set,
		const HybridConstraintSet& constraint_set,
		const HybridGridTreeSet& reachability_restriction,
		EvolutionDirection direction) const
{
	int& maximum_grid_depth = _outer_analyser->settings().maximum_grid_depth;

	HybridGridTreeSet result(*_outer_analyser->settings().grid);

	if (direction == DIRECTION_FORWARD) {
		result.adjoin_outer_approximation(initial_set,maximum_grid_depth);
		result.mince(maximum_grid_depth);
		result.restrict(reachability_restriction);
	} else {
		result = reachability_restriction;
		result.mince(maximum_grid_depth);
		result.remove(covered_cells(result,constraint_set));
	}

	ARIADNE_LOG(6,"Starting set size: " << result.size() << "\n");

	if (_settings->plot_results)
		_plot_reach(result,(direction == DIRECTION_FORWARD ? "initial" : "final"),maximum_grid_depth);

	return result;
}


void
Verifier::
_update_safety_cached_reachability_with(const HybridGridTreeSet& reach) const
{
	if (_safety_coarse_outer_approximation->is_set())
		_safety_reachability_restriction = reach;
	else
		_safety_coarse_outer_approximation->set(reach);
}


bool
Verifier::
_safety_disproving_once(
		HybridAutomaton& system,
		const HybridImageSet& initial_set,
		const HybridConstraintSet& safety_constraint,
		const RealConstantSet& constants) const
{
	ARIADNE_LOG(4,"Disproving...\n");
	if (!_grid_depth_is_within_bounds_under(LOWER_SEMANTICS)) {
		ARIADNE_LOG(4,"Not disproved.\n");
		return false;
	}

	bool result = false;

	RealConstantSet original_constants = system.accessible_constants();

	system.substitute(constants,_settings->use_param_midpoints_for_disproving);

	ARIADNE_LOG(4,"Setting parameters for this disproving iteration...\n");

	_tuneIterativeStepSettings(system,_safety_coarse_outer_approximation->get(),_safety_reachability_restriction,LOWER_SEMANTICS);

	ARIADNE_LOG(4,"Performing lower reachability analysis and getting disprove data...\n");

	try {
		std::pair<HybridGridTreeSet,DisproveData> reachAndDisproveData =
				_lower_analyser->lower_chain_reach(system,initial_set,safety_constraint,_safety_reachability_restriction);
		const HybridGridTreeSet& reach = reachAndDisproveData.first;
		const DisproveData& disproveData = reachAndDisproveData.second;

		ARIADNE_LOG(5,"Epsilon: " << disproveData.getEpsilon() << "\n");

		HybridGridTreeSet reachability_restriction;
		reachability_restriction.adjoin_outer_approximation(_lower_analyser->settings().domain_bounds,_lower_analyser->settings().maximum_grid_depth);

		HybridGridTreeSet possibly_safe_cells = _lower_analyser->possibly_feasible_cells(
				reach,safety_constraint,disproveData.getEpsilon(),reachability_restriction);

		result = (possibly_safe_cells.size() < reach.size());

		if (_settings->plot_results)
			_plot_reach(reach,"lower",_lower_analyser->settings().maximum_grid_depth);

	} catch (ReachUnsatisfiesConstraintException ex) {
		ARIADNE_LOG(5, "The lower reached region is partially outside the safe region (skipped).\n");
		result = true;
	}

	system.substitute(original_constants);

	ARIADNE_LOG(4, (result ? "Disproved.\n" : "Not disproved.\n") );

	return result;
}


bool
Verifier::
_grid_depth_is_within_bounds_under(Semantics semantics) const
{
	const DiscreteEvolutionSettings& settings = (semantics == UPPER_SEMANTICS ? _outer_analyser->settings() : _lower_analyser->settings());

	if (settings.maximum_grid_depth < settings.lowest_maximum_grid_depth) {
		ARIADNE_LOG(4,"Will skip verification since the depth is lower than the lowest allowed.\n");
		return false;
	}
	if (settings.maximum_grid_depth > settings.highest_maximum_grid_depth) {
		ARIADNE_LOG(4,"Will skip verification since the depth is higher than the highest allowed.\n");
		return false;
	}

	return true;
}


std::list<ParametricOutcome>
Verifier::
parametric_safety(
		SafetyVerificationInput& verInput,
		const RealConstantSet& params) const
{
	ARIADNE_ASSERT_MSG(params.size() > 0, "Provide at least one parameter.");

	std::list<ParametricOutcome> result;

	std::list<RealConstantSet> splittings = maximally_split_parameters(params,_settings->maximum_parameter_depth);
	uint i=0;
	for (std::list<RealConstantSet>::const_iterator splitting_it = splittings.begin();
													 splitting_it != splittings.end();
													 ++splitting_it)
	{
		ARIADNE_LOG(1,"<Split parameters set #" << ++i << "/" << splittings.size() << ">\n");
		RealConstantSet current_params = *splitting_it;
		ARIADNE_LOG(1,"Parameter values: " << current_params << " ");
		tribool outcome = _safety_nosplitting(verInput,current_params);
		ARIADNE_LOG(1,"Outcome: " << pretty_print(outcome) << "\n");
		result.push_back(ParametricOutcome(current_params,outcome));
	}

	return result;
}

tribool
Verifier::
dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated) const
{
	const RealConstantSet dominatingConstants;
	return _dominance(dominating,dominated,dominatingConstants);
}

tribool
Verifier::
_dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealConstant& constant) const
{
	HybridAutomaton& system = dominating.getSystem();

	Real original_value = system.accessible_constant_value(constant.name());

	system.substitute(constant);
	tribool result = dominance(dominating,dominated);
	system.substitute(constant,original_value);

	return result;
}

tribool
Verifier::
_dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealConstant& constant,
		const Float& value) const
{
	const RealConstant modifiedConstant(constant.name(),Interval(value));

	return _dominance(dominating,dominated,modifiedConstant);
}


std::list<ParametricOutcome>
Verifier::
parametric_dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealConstantSet& dominating_params) const
{
	ARIADNE_ASSERT_MSG(dominating_params.size() > 0, "Provide at least one parameter.");

	std::list<ParametricOutcome> result;

	std::list<RealConstantSet> splittings = maximally_split_parameters(dominating_params,_settings->maximum_parameter_depth);
	uint i=0;
	for (std::list<RealConstantSet>::const_iterator splitting_it = splittings.begin(); splitting_it != splittings.end(); ++splitting_it) {
		ARIADNE_LOG(1,"<Split parameters set #" << ++i << "/" << splittings.size() << ">\n");
		RealConstantSet current_params = *splitting_it;
		ARIADNE_LOG(1,"Parameter values: " << current_params << " ");
		tribool outcome = _dominance(dominating,dominated,current_params);
		ARIADNE_LOG(1,"Outcome: " << pretty_print(outcome) << "\n");
		result.push_back(ParametricOutcome(current_params,outcome));
	}

	return result;
}


tribool
Verifier::_dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealConstantSet& constants) const
{
	ARIADNE_ASSERT(dominating.getProjection().size() == dominated.getProjection().size());

	ARIADNE_LOG(1, "Dominance checking...\n");

	if (_settings->plot_results)
		_plot_dirpath_init(dominating.getSystem().name() + "&" + dominated.getSystem().name());

	_resetAndChooseInitialDominanceSettings();

	int initial_depth = max(_outer_analyser->settings().lowest_maximum_grid_depth,
							_lower_analyser->settings().lowest_maximum_grid_depth);
	int final_depth = min(_outer_analyser->settings().highest_maximum_grid_depth,
						  _lower_analyser->settings().highest_maximum_grid_depth);
    for (int depth = initial_depth; depth <= final_depth; ++depth)
	{
    	_outer_analyser->settings().maximum_grid_depth = depth;
    	_lower_analyser->settings().maximum_grid_depth = depth;

		ARIADNE_LOG(2, "Depth " << depth << "\n");

		if (_dominance_proving_once(dominating, dominated, constants)) {
			ARIADNE_LOG(3, "Dominates.\n");
			return true;
		}

		if (_dominance_disproving_once(dominating, dominated, constants)) {
			ARIADNE_LOG(3, "Does not dominate.\n");
			return false;
		}
    }

	ARIADNE_LOG(3, "Indeterminate.\n");
	return indeterminate;
}

bool
Verifier::_dominance_proving_once(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealConstantSet& constants) const
{
	ARIADNE_LOG(3,"Proving...\n");

	bool result;

	const RealConstantSet& original_constants = dominating.getSystem().accessible_constants();

	dominating.getSystem().substitute(constants,_settings->use_param_midpoints_for_proving);

	try {

		Box shrinked_dominated_bounds = _dominance_shrinked_lower_bounds(dominated,constants,DOMINATED_SYSTEM);

		HybridBoxes shrinked_dominated_bounds_on_dominating_space = Ariadne::project(shrinked_dominated_bounds,
				dominating.getProjection(),dominating.getSystem().state_space());

		Box dominating_bounds = _dominance_outer_bounds(
				dominating,shrinked_dominated_bounds_on_dominating_space,constants,DOMINATING_SYSTEM);

		result = inside(dominating_bounds,shrinked_dominated_bounds);

		ARIADNE_LOG(4, "The outer reached region of the dominating system is " << (!result ? "not ": "") <<
					   "inside the projected shrinked lower reached region of the dominated system.\n");

	} catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominating system is partially out of the domain (skipped).\n");
		result = false;
	} catch (ReachUnsatisfiesConstraintException ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominating system is not inside " +
				"the projected shrinked lower reached region of the dominated system (skipped).\n");
		result = false;
	}

	ARIADNE_LOG(3, (result ? "Proved.\n" : "Not proved.\n") );

	dominating.getSystem().substitute(original_constants);

	return result;
}

bool
Verifier::_dominance_disproving_once(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
	  	const RealConstantSet& constants) const
{
	ARIADNE_LOG(3,"Disproving...\n");

	bool result;

	const RealConstantSet& original_constants = dominating.getSystem().accessible_constants();

	dominating.getSystem().substitute(constants,_settings->use_param_midpoints_for_disproving);

	try {

		Box shrinked_dominating_bounds = _dominance_shrinked_lower_bounds(dominating,constants,DOMINATING_SYSTEM);

		HybridBoxes shrinked_dominating_bounds_on_dominated_space = Ariadne::project(shrinked_dominating_bounds,
				dominated.getProjection(),dominated.getSystem().state_space());

		Box dominated_bounds = _dominance_outer_bounds(
				dominated,shrinked_dominating_bounds_on_dominated_space,constants,DOMINATED_SYSTEM);

		result = !inside(shrinked_dominating_bounds,dominated_bounds);

	} catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominated system is partially out of the domain (skipped).\n");
		result = false;
	} catch (ReachEnclosesTargetException ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominated system encloses " +
				"the projected shrinked lower reached region of the dominated system (skipped).\n");
		result = true;
	}

	ARIADNE_LOG(3, (result ? "Disproved.\n" : "Not disproved.\n") );

	dominating.getSystem().substitute(original_constants);

	return result;
}


Box
Verifier::
_dominance_shrinked_lower_bounds(
		DominanceVerificationInput& verInfo,
		const RealConstantSet& constants,
		DominanceSystem dominanceSystem) const
{
	string descriptor = (dominanceSystem == DOMINATING_SYSTEM ? "dominating" : "dominated");
	HybridGridTreeSet outer_approximation = (dominanceSystem == DOMINATING_SYSTEM ?
			_dominating_coarse_outer_approximation->get() : _dominated_coarse_outer_approximation->get());
	HybridGridTreeSet reachability_restriction = (dominanceSystem == DOMINATING_SYSTEM ?
			_dominating_reachability_restriction : _dominated_reachability_restriction);

	ARIADNE_LOG(4,"Choosing the settings for the lower reached region of the " << descriptor << " system...\n");

	RealConstantSet emptyLockedConstants;
	_chooseDominanceSettings(verInfo,emptyLockedConstants,outer_approximation,reachability_restriction,LOWER_SEMANTICS);

	ARIADNE_LOG(4,"Getting the lower reached region of the " << descriptor << " system...\n");

	HybridGridTreeSet reach;
	DisproveData disproveData(verInfo.getSystem().state_space());
	make_lpair<HybridGridTreeSet,DisproveData>(reach,disproveData) =
			_lower_analyser->lower_chain_reach(verInfo.getSystem(),verInfo.getInitialSet(),reachability_restriction);

	// We must shrink the lower approximation of the system, but underapproximating in terms of rounding
	HybridBoxes shrinked_bounds = Ariadne::shrink_in(disproveData.getReachBounds(),disproveData.getEpsilon());

	Box projected_shrinked_bounds = Ariadne::project(shrinked_bounds,verInfo.getProjection());

	ARIADNE_LOG(5,"Epsilon: " << disproveData.getEpsilon() << "\n");
	ARIADNE_LOG(5,"Projected shrinked " << descriptor << " bounds: " << projected_shrinked_bounds << "\n");

	if (_settings->plot_results)
		_plot_dominance(reach,dominanceSystem,LOWER_SEMANTICS);

	return projected_shrinked_bounds;
}


Box
Verifier::
_dominance_outer_bounds(
		DominanceVerificationInput& verInput,
		HybridBoxes& lower_bounds_on_this_space,
		const RealConstantSet& constants,
		DominanceSystem dominanceSystem) const
{
	string descriptor = (dominanceSystem == DOMINATING_SYSTEM ? "dominating" : "dominated");
	OuterApproximationCache& outer_approximation_cache = (dominanceSystem == DOMINATING_SYSTEM ?
			*_dominating_coarse_outer_approximation : *_dominated_coarse_outer_approximation);
	HybridGridTreeSet& reachability_restriction = (dominanceSystem == DOMINATING_SYSTEM ?
			_dominating_reachability_restriction : _dominated_reachability_restriction);

	ARIADNE_LOG(4,"Choosing the settings for the outer reached region of the " << descriptor << " system...\n");

	_chooseDominanceSettings(verInput,constants,outer_approximation_cache.get(),reachability_restriction,UPPER_SEMANTICS);

	ARIADNE_LOG(4,"Getting the outer reached region of the " << descriptor << " system...\n");

	HybridGridTreeSet reach = _outer_analyser->outer_chain_reach(
			verInput.getSystem(),verInput.getInitialSet(),DIRECTION_FORWARD,reachability_restriction);

	if (!outer_approximation_cache.is_set()) {
		outer_approximation_cache.set(reach);
	} else {
		reachability_restriction = reach;
	}

	Box projected_bounds = Ariadne::project(reach.bounding_box(),verInput.getProjection());

	ARIADNE_LOG(5,"Projected " << descriptor << " bounds: " << projected_bounds << "\n");

	if (_settings->plot_results)
		_plot_dominance(reach,dominanceSystem,UPPER_SEMANTICS);

	return projected_bounds;
}


void
Verifier::
_resetAndChooseInitialSafetySettings(
		const HybridAutomaton& system,
		const HybridBoxes& domain,
		const RealConstantSet& locked_constants) const
{
	_safety_coarse_outer_approximation->reset();

	_safety_reachability_restriction = HybridGridTreeSet();

	_chooseInitialSafetySettings(system,domain,locked_constants,UPPER_SEMANTICS);
	_chooseInitialSafetySettings(system,domain,locked_constants,LOWER_SEMANTICS);
}


void
Verifier::
_chooseInitialSafetySettings(
		const HybridAutomaton& system,
		const HybridBoxes& domain,
		const RealConstantSet& locked_constants,
		Semantics semantics) const
{
	ARIADNE_LOG(3,"Choosing the initial settings of the " << (semantics == UPPER_SEMANTICS ? "outer " : "lower ") << "analyser...\n");
	DiscreteEvolutionSettings& settings = (semantics == UPPER_SEMANTICS ?
											   _outer_analyser->settings() : _lower_analyser->settings());

	settings.domain_bounds = domain;
	ARIADNE_LOG(4, "Domain: " << domain << "\n");
	settings.lock_to_grid_time = getLockToGridTime(system,domain);
	ARIADNE_LOG(4, "Lock to grid time: " << settings.lock_to_grid_time << "\n");
	settings.locked_constants = locked_constants;
	ARIADNE_LOG(4, "Locked constants: " << locked_constants << "\n");
}


void
Verifier::
_tuneIterativeStepSettings(
		const HybridAutomaton& system,
		const HybridGridTreeSet& hgts_domain,
		const HybridGridTreeSet& reachability_restriction,
		Semantics semantics) const
{
	HybridReachabilityAnalyser& analyser = (semantics == UPPER_SEMANTICS ? *_outer_analyser : *_lower_analyser);

	ARIADNE_LOG(5, "Derivatives evaluation source: " << (hgts_domain.empty() ? "Domain box" : "Outer approximation") << "\n");

	HybridFloatVector hmad = getHybridMaximumAbsoluteDerivatives(system,hgts_domain,analyser.settings().domain_bounds);
	ARIADNE_LOG(5, "Derivatives bounds: " << hmad << "\n");
	analyser.settings().grid = boost::shared_ptr<HybridGrid>(
			new HybridGrid(getHybridGrid(hmad,analyser.settings().domain_bounds)));
	ARIADNE_LOG(5, "Grid lengths: " << analyser.settings().grid->lengths() << "\n");

	ARIADNE_LOG(5, "Use reachability restriction: " << pretty_print(!reachability_restriction.empty()) << "\n");

	analyser.tuneEvolverSettings(system,hmad,analyser.settings().maximum_grid_depth,semantics);
}

void
Verifier::
_resetAndChooseInitialDominanceSettings() const
{
	_dominating_coarse_outer_approximation->reset();
	_dominated_coarse_outer_approximation->reset();

	_dominating_reachability_restriction = HybridGridTreeSet();
	_dominated_reachability_restriction = HybridGridTreeSet();
}

void
Verifier::
_chooseDominanceSettings(
		const DominanceVerificationInput& verInput,
		const RealConstantSet& locked_constants,
		const HybridGridTreeSet& outer_reach,
		const HybridGridTreeSet& outer_approx_constraint,
		Semantics semantics) const
{
	HybridReachabilityAnalyser& analyser = (semantics == UPPER_SEMANTICS ? *_outer_analyser : *_lower_analyser);

	analyser.settings().domain_bounds = verInput.getDomain();
	ARIADNE_LOG(5, "Domain: " << analyser.settings().domain_bounds << "\n");

	_tuneIterativeStepSettings(verInput.getSystem(),outer_reach,outer_approx_constraint,semantics);

	analyser.settings().lock_to_grid_time = getLockToGridTime(verInput.getSystem(),analyser.settings().domain_bounds);
	ARIADNE_LOG(5, "Lock to grid time: " << analyser.settings().lock_to_grid_time << "\n");
	analyser.settings().locked_constants = locked_constants;
	ARIADNE_LOG(5, "Locked constants: " << analyser.settings().locked_constants << "\n");
}

void
Verifier::
_plot_dirpath_init(std::string basename) const
{
	time_t mytime;
	time(&mytime);
	string foldername = basename+"-png";

	mkdir(foldername.c_str(),0777);
	string timestring = asctime(localtime(&mytime));
	timestring.erase(std::remove(timestring.begin(), timestring.end(), '\n'), timestring.end());
	foldername = foldername+"/"+timestring;
	mkdir(foldername.c_str(),0777);

	_plot_dirpath = foldername;
}


void
Verifier::
_plot_reach(
		const HybridGridTreeSet& reach,
		string base_filename,
		int accuracy) const
{
	char mgd_char[10];
	sprintf(mgd_char,"%i",accuracy);
	base_filename.append(mgd_char);
	plot(_plot_dirpath,base_filename,reach);
}


void
Verifier::
_plot_dominance(
		const HybridGridTreeSet& reach,
		DominanceSystem dominanceSystem,
		Semantics semantics) const
{
	int maximum_grid_depth = _outer_analyser->settings().maximum_grid_depth;

	string system_descr = (dominanceSystem == DOMINATING_SYSTEM ? "dominating" : "dominated");
	string verification_descr = (semantics == UPPER_SEMANTICS ? "pos" : "neg");

	char mgd_char[10];
	sprintf(mgd_char,"%i",maximum_grid_depth);
	string filename = system_descr + "-" + verification_descr + "-";
	filename.append(mgd_char);
	plot(_plot_dirpath,filename,reach);
}


std::string
pretty_print(tribool value)
{
	if (definitely(value))
		return "True";
	if (!possibly(value))
		return "False";
	return "Indeterminate";
}


std::list<RealConstantSet>
maximally_split_parameters(
		const RealConstantSet& params,
		const uint& maximum_parameter_depth)
{
	std::list<RealConstantSet> source;
	std::list<RealConstantSet> destination;
	destination.push_back(params);

	for (RealConstantSet::const_iterator param_it = params.begin(); param_it != params.end(); ++param_it)
	{
		for (uint i=0; i<maximum_parameter_depth; i++) {
			source.clear();
			source.insert<std::list<RealConstantSet>::const_iterator>(source.begin(),destination.begin(),destination.end());
			destination.clear();

			while (!source.empty()) {
				RealConstantSet currentParams = source.back();
				source.pop_back();

				RealConstantSet newConfigurationLeft = currentParams;
				RealConstantSet newConfigurationRight = currentParams;
				newConfigurationLeft.erase(*param_it);
				newConfigurationRight.erase(*param_it);

				const Real& currentInterval = currentParams.find(*param_it)->value();
				newConfigurationLeft.insert(RealConstant(param_it->name(),
						Interval(currentInterval.lower(),currentInterval.midpoint())));
				newConfigurationRight.insert(RealConstant(param_it->name(),
						Interval(currentInterval.midpoint(),currentInterval.upper())));

				destination.push_back(newConfigurationLeft);
				destination.push_back(newConfigurationRight);
			}
		}
	}

	return destination;
}

}
