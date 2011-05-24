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


Verifier::Verifier(const HybridReachabilityAnalyser& analyser) :
			_analyser(analyser.clone()),
			_settings(new VerificationSettings()),
			_safety_coarse_outer_approximation(new OuterApproximationCache()),
			_dominating_coarse_outer_approximation(new OuterApproximationCache()),
			_dominated_coarse_outer_approximation(new OuterApproximationCache())
{
	_analyser->verb_tab_prefix = verb_tab_prefix + verifier_max_verbosity_level_used;
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

	_analyser->settings().maximum_grid_depth = _analyser->settings().lowest_maximum_grid_depth;

	while (_analyser->settings().maximum_grid_depth <= _analyser->settings().highest_maximum_grid_depth) {

		ARIADNE_LOG(2, "Depth " << _analyser->settings().maximum_grid_depth << "\n");

		tribool result = _safety_once(system,verInput.getInitialSet(),verInput.getSafetyConstraint(),constants);

    	_analyser->settings().maximum_grid_depth++;

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

	int& maximum_grid_depth = _analyser->settings().maximum_grid_depth;

	RealConstantSet original_constants = system.accessible_constants();

	system.substitute(constants,_settings->use_param_midpoints_for_proving);

	ARIADNE_LOG(4,"Setting parameters for this proving iteration...\n");

	bool equal_grid_for_all_locations = false;
	_tuneIterativeStepSettings(system,_safety_coarse_outer_approximation->get(),
			equal_grid_for_all_locations,UPPER_SEMANTICS);

	ARIADNE_LOG(5, "Using reachability restriction: " << pretty_print(!_safety_reachability_restriction.empty()) << "\n");

	ARIADNE_LOG(4,"Performing outer reachability analysis...\n");

	try {

		/* Given the initial set I, forward reachability F, and unsafe region U:

		1. Evaluate the final unsafe set H = F ∩ U ; if H = ∅, return TRUE;
		2. If the accuracy is over the maximum, return FALSE;
		3. Obtain the backward reachability B under the constraint set F;
		4. If B = F, return FALSE;
		5. Evaluate the new initial set I = B ∩ I; if I = ∅, return TRUE;
		6. If the accuracy is over the maximum, return FALSE;
		7. Starting from I, obtain the forward reachability F under the constraint set B;
		8. If F = B, return FALSE;

		*/

		if (!_safety_reachability_restriction.empty() && _settings->enable_backward_refinement_for_testing_inclusion) {

			HybridGridTreeSet backward_initial = _reachability_refinement_starting_set(system,initial_set,
					safety_constraint,_safety_reachability_restriction,DIRECTION_BACKWARD);

			if (backward_initial.empty()) {
				ARIADNE_LOG(4, "The initial set for backward reachability is empty.\n");
				system.substitute(original_constants);
				return true;
			}

			ARIADNE_LOG(5,"Retrieving backward reachability...\n");

			HybridGridTreeSet backward_reach = _analyser->outer_chain_reach(system,backward_initial,
					DIRECTION_BACKWARD,_safety_reachability_restriction);

			_safety_reachability_restriction = backward_reach;

			ARIADNE_LOG(6,"Reachability size: " << backward_reach.size() << "\n");

			if (_settings->plot_results)
				_plot_reach(backward_reach,"backward",maximum_grid_depth);

		}

		HybridGridTreeSet forward_initial;
		if (!_safety_reachability_restriction.empty()) {
			forward_initial = _reachability_refinement_starting_set(system,initial_set,
				safety_constraint,_safety_reachability_restriction,DIRECTION_FORWARD);
		} else {
			forward_initial.adjoin_outer_approximation(initial_set,maximum_grid_depth);
			forward_initial.mince(maximum_grid_depth);
		}

		if (forward_initial.empty()) {
			ARIADNE_LOG(4, "The initial set for forward reachability is empty.\n");
			system.substitute(original_constants);
			return true;
		}

		ARIADNE_LOG(5,"Retrieving forward reachability...\n");

		HybridGridTreeSet forward_reach = _analyser->outer_chain_reach(system,forward_initial,
				DIRECTION_FORWARD,_safety_reachability_restriction);

		ARIADNE_LOG(6,"Reachability size: " << forward_reach.size() << "\n");

		if (_settings->plot_results)
			_plot_reach(forward_reach,"forward",maximum_grid_depth);

		_update_safety_cached_reachability_with(forward_reach);

		// Additional simplified check, also useful when only forward reachability is available
		result = definitely(covers(safety_constraint,forward_reach));

	} catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(4, "The outer reached region is partially out of the domain (skipped).\n");
		result = false;
	} catch (ReachUnsatisfiesConstraintException ex) {
		ARIADNE_LOG(4, "The outer reached region is not inside the safe region (skipped).\n");
		result = false;
	}

	system.substitute(original_constants);

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
	int& maximum_grid_depth = _analyser->settings().maximum_grid_depth;

	HybridGridTreeSet result(*_analyser->settings().grid);

	if (direction == DIRECTION_FORWARD) {
		result.adjoin_outer_approximation(initial_set,maximum_grid_depth);
		result.mince(maximum_grid_depth);
		result.restrict(reachability_restriction);
	} else {
		result = reachability_restriction;
		result.mince(maximum_grid_depth);
		result.remove(definitely_covered_cells(result,constraint_set));
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

	bool result = false;

	RealConstantSet original_constants = system.accessible_constants();

	system.substitute(constants,_settings->use_param_midpoints_for_disproving);

	ARIADNE_LOG(4,"Setting parameters for this disproving iteration...\n");

	bool equal_grid_for_all_locations = false;
	_tuneIterativeStepSettings(system,_safety_coarse_outer_approximation->get(),
			equal_grid_for_all_locations,LOWER_SEMANTICS);

	ARIADNE_LOG(5, "Using reachability restriction: " << pretty_print(!_safety_reachability_restriction.empty()) << "\n");

	ARIADNE_LOG(4,"Performing lower reachability analysis...\n");

	try {

		HybridGridTreeSet reach;
		HybridFloatVector epsilon;
		make_lpair<HybridGridTreeSet,HybridFloatVector>(reach,epsilon) = _analyser->lower_reach_and_epsilon(
				system,initial_set,safety_constraint,_safety_reachability_restriction);

		ARIADNE_LOG(5, "Epsilon: " << epsilon << "\n");

		if (_settings->plot_results)
			_plot_reach(reach,"lower",_analyser->settings().maximum_grid_depth);

	} catch (ReachUnsatisfiesConstraintException ex) {
		ARIADNE_LOG(5, "The lower reached region is partially outside the safe region (skipped).\n");
		result = true;
	}

	system.substitute(original_constants);

	ARIADNE_LOG(4, (result ? "Disproved.\n" : "Not disproved.\n") );

	return result;
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

	_resetAndChooseInitialDominanceSettings(dominating,dominated);

    for (int depth = _analyser->settings().lowest_maximum_grid_depth;
    		depth <= _analyser->settings().highest_maximum_grid_depth; ++depth)
	{
    	_analyser->settings().maximum_grid_depth = depth;

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
		GridTreeSet flattened_dominated_lower_reach;
		Vector<Float> flattened_epsilon;
		make_lpair<GridTreeSet,Vector<Float> >(flattened_dominated_lower_reach,flattened_epsilon) =
				_dominance_flattened_lower_reach_and_epsilon(dominated,constants,DOMINATED_SYSTEM);

		GridTreeSet flattened_dominating_outer_reach = _dominance_flattened_outer_reach(dominating,constants,DOMINATING_SYSTEM);

		result = definitely(covers(flattened_dominated_lower_reach,flattened_dominating_outer_reach,flattened_epsilon));

		ARIADNE_LOG(4, "The outer reached region of the dominating system is " << (!result ? "not ": "") <<
					   "inside the projected shrinked lower reached region of the dominated system.\n");

	} catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominating system is partially out of the domain (skipped).\n");
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
		GridTreeSet flattened_dominating_lower_reach;
		Vector<Float> flattened_epsilon;
		make_lpair<GridTreeSet,Vector<Float> >(flattened_dominating_lower_reach,flattened_epsilon) =
				_dominance_flattened_lower_reach_and_epsilon(dominating,constants,DOMINATING_SYSTEM);

		GridTreeSet flattened_dominated_outer_reach = _dominance_flattened_outer_reach(dominated,constants,DOMINATED_SYSTEM);

		result = definitely(!inside(flattened_dominating_lower_reach,flattened_dominated_outer_reach,
				flattened_epsilon,_analyser->settings().maximum_grid_depth));

	} catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominated system is partially out of the domain (skipped).\n");
		result = false;
	}

	ARIADNE_LOG(3, (result ? "Disproved.\n" : "Not disproved.\n") );

	dominating.getSystem().substitute(original_constants);

	return result;
}


std::pair<GridTreeSet,Vector<Float> >
Verifier::
_dominance_flattened_lower_reach_and_epsilon(
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

	const Vector<uint>& projection = verInfo.getProjection();

	HybridGridTreeSet reach;
	HybridFloatVector epsilon;
	make_lpair<HybridGridTreeSet,HybridFloatVector>(reach,epsilon) = _analyser->lower_reach_and_epsilon(
			verInfo.getSystem(),verInfo.getInitialSet(),reachability_restriction);

	GridTreeSet flattened_reach = flatten_and_project_down(reach,projection);

	HybridFloatVector::const_iterator epsilon_it = epsilon.begin();
	Vector<Float> flattened_epsilon(projection.size());
	for (HybridFloatVector::const_iterator epsilon_it = epsilon.begin(); epsilon_it != epsilon.end(); ++epsilon_it) {
		for (uint i=0; i<projection.size(); ++i)
			flattened_epsilon[i] = max(flattened_epsilon[i],epsilon_it->second[projection[i]]);
	}

	ARIADNE_LOG(5,"Flattened epsilon: " << flattened_epsilon << "\n");

	return std::pair<GridTreeSet,Vector<Float> >(flattened_reach,flattened_epsilon);
}


GridTreeSet
Verifier::
_dominance_flattened_outer_reach(
		DominanceVerificationInput& verInput,
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

	HybridGridTreeSet reach = _analyser->outer_chain_reach(
			verInput.getSystem(),verInput.getInitialSet(),DIRECTION_FORWARD,reachability_restriction);

	if (!outer_approximation_cache.is_set()) {
		outer_approximation_cache.set(reach);
	} else {
		reachability_restriction = reach;
	}

	if (_settings->plot_results)
		_plot_dominance(reach,dominanceSystem,UPPER_SEMANTICS);

	return flatten_and_project_down(reach,verInput.getProjection());
}


void
Verifier::
_resetAndChooseInitialSafetySettings(
		const HybridAutomaton& system,
		const HybridBoxes& domain,
		const RealConstantSet& locked_constants) const
{
	static const bool equal_grid_for_all_locations = false;

	_safety_coarse_outer_approximation->reset();

	_safety_reachability_restriction = HybridGridTreeSet();

	_chooseInitialSafetySettings(system,domain,locked_constants);

	if (_settings->enable_domain_enforcing) {
		std::pair<HybridGridTreeSet,HybridGridTreeSet> reach_pair =
				_getCoarseOuterApproximationAndReachabilityRestriction(system,domain,equal_grid_for_all_locations);
		_safety_coarse_outer_approximation->set(reach_pair.first);
		_safety_reachability_restriction = reach_pair.second;
	}
}


std::pair<HybridGridTreeSet,HybridGridTreeSet>
Verifier::
_getCoarseOuterApproximationAndReachabilityRestriction(
		const HybridAutomaton& system,
		const HybridBoxes& domain,
		bool equal_grid_for_all_locations) const
{
	static const int accuracy = 2;

	HybridFloatVector hmad = getHybridMaximumAbsoluteDerivatives(system,HybridGridTreeSet(),domain);
	HybridGrid coarse_grid = getHybridGrid(hmad,domain,equal_grid_for_all_locations);

	HybridGridTreeSet coarse_outer_approximation(coarse_grid);
	coarse_outer_approximation.adjoin_outer_approximation(domain,accuracy);

	hmad = getHybridMaximumAbsoluteDerivatives(system,coarse_outer_approximation,domain);
	HybridGrid fine_grid = getHybridGrid(hmad,domain,equal_grid_for_all_locations);
	HybridGridTreeSet reachability_restriction(fine_grid);
	reachability_restriction.adjoin_outer_approximation(domain,accuracy);

	return std::pair<HybridGridTreeSet,HybridGridTreeSet>(coarse_outer_approximation,reachability_restriction);
}


void
Verifier::
_chooseInitialSafetySettings(
		const HybridAutomaton& system,
		const HybridBoxes& domain,
		const RealConstantSet& locked_constants) const
{
	ARIADNE_LOG(3,"Choosing the initial settings of the analyser...\n");
	DiscreteEvolutionSettings& settings = _analyser->settings();

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
		bool equal_grid_for_all_locations,
		Semantics semantics) const
{
	ARIADNE_LOG(5, "Derivatives evaluation source: " << (hgts_domain.empty() ? "Domain box" : "Outer approximation") << "\n");

	HybridFloatVector hmad = getHybridMaximumAbsoluteDerivatives(system,hgts_domain,_analyser->settings().domain_bounds);
	ARIADNE_LOG(5, "Derivatives bounds: " << hmad << "\n");
	_analyser->settings().grid = boost::shared_ptr<HybridGrid>(
			new HybridGrid(getHybridGrid(hmad,_analyser->settings().domain_bounds,equal_grid_for_all_locations)));
	ARIADNE_LOG(5, "Grid lengths: " << _analyser->settings().grid->lengths() << "\n");

	_analyser->tuneEvolverSettings(system,hmad,_analyser->settings().maximum_grid_depth,semantics);
}

void
Verifier::
_resetAndChooseInitialDominanceSettings(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated) const
{
	static const bool equal_grid_for_all_locations = true;

	_dominating_coarse_outer_approximation->reset();
	_dominated_coarse_outer_approximation->reset();

	_dominating_reachability_restriction = HybridGridTreeSet();
	_dominated_reachability_restriction = HybridGridTreeSet();

	if (_settings->enable_domain_enforcing) {
		std::pair<HybridGridTreeSet,HybridGridTreeSet> dominating_reach_pair =
				_getCoarseOuterApproximationAndReachabilityRestriction(dominating.getSystem(),
						dominating.getDomain(),equal_grid_for_all_locations);
		_dominating_coarse_outer_approximation->set(dominating_reach_pair.first);
		_dominating_reachability_restriction = dominating_reach_pair.second;

		std::pair<HybridGridTreeSet,HybridGridTreeSet> dominated_reach_pair =
				_getCoarseOuterApproximationAndReachabilityRestriction(dominated.getSystem(),
						dominated.getDomain(),equal_grid_for_all_locations);
		_dominated_coarse_outer_approximation->set(dominated_reach_pair.first);
		_dominated_reachability_restriction = dominated_reach_pair.second;
	}
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
	_analyser->settings().domain_bounds = verInput.getDomain();
	ARIADNE_LOG(5, "Domain: " << _analyser->settings().domain_bounds << "\n");

	bool equal_grid_for_all_locations = true;
	_tuneIterativeStepSettings(verInput.getSystem(),outer_reach,equal_grid_for_all_locations,semantics);

	ARIADNE_LOG(5, "Use reachability restriction: " << pretty_print(!outer_approx_constraint.empty()) << "\n");

	_analyser->settings().lock_to_grid_time = getLockToGridTime(verInput.getSystem(),_analyser->settings().domain_bounds);
	ARIADNE_LOG(5, "Lock to grid time: " << _analyser->settings().lock_to_grid_time << "\n");
	_analyser->settings().locked_constants = locked_constants;
	ARIADNE_LOG(5, "Locked constants: " << _analyser->settings().locked_constants << "\n");
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
	int maximum_grid_depth = _analyser->settings().maximum_grid_depth;

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
