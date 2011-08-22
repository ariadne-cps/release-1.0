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


const unsigned int VERIFIER_CHILD_OFFSET = 5;

Verifier::Verifier() :
			_settings(new VerificationSettings()),
			_safety_coarse_outer_approximation(new OuterApproximationCache()),
			_dominating_coarse_outer_approximation(new OuterApproximationCache()),
			_dominated_coarse_outer_approximation(new OuterApproximationCache())
{
    this->charcode = "v";
    this->child_tab_offset = VERIFIER_CHILD_OFFSET;
}



Verifier::~Verifier()
{
}

tribool
Verifier::
safety(SafetyVerificationInput& verInput) const
{
	RealParameterSet params;

	return _safety_nosplitting(verInput,params);
}


tribool
Verifier::
_safety(
		SafetyVerificationInput& verInput,
		const RealParameter& param) const
{
	SystemType& system = verInput.getSystem();

	Real originalParameterValue = system.parameter_value(param.name());

	system.substitute(param);
	tribool result = safety(verInput);
	system.substitute(RealParameter(param.name(),originalParameterValue));

	return result;
}

tribool
Verifier::
_safety(
		SafetyVerificationInput& verInput,
		const RealParameter& param,
		const Float& value) const
{
	const RealParameter modifiedParameter(param.name(),Interval(value));

	return _safety(verInput, modifiedParameter);
}


tribool
Verifier::
_safety_nosplitting(
		SafetyVerificationInput& verInput,
		const RealParameterSet& params) const
{
	ARIADNE_LOG(2,"Iterative verification...");

	SystemType& system = verInput.getSystem();

	if (_settings->plot_results)
		_plot_dirpath_init(system.name());

	_reset_safety_state();

	time_t initial_time = time(NULL);
	for (int accuracy = 0; time(NULL) - initial_time < _settings->time_limit_for_outcome; ++accuracy) {

		ARIADNE_LOG(2, "Accuracy " << accuracy);

		tribool result = _safety_once(verInput,accuracy,params);

		if (!indeterminate(result))
			return result;
    }

	return indeterminate;
}


tribool
Verifier::
_safety_once(SafetyVerificationInput& verInput,
		const unsigned int& accuracy,
		const RealParameterSet& params) const
{
    ARIADNE_LOG(3, "Verification...");

    if (_safety_proving_once(verInput,accuracy,params)) {
        ARIADNE_LOG(3, "Safe.");
        return true;
    }

    if (_safety_disproving_once(verInput,accuracy,params)) {
        ARIADNE_LOG(3, "Unsafe.");
        return false;
    }

    ARIADNE_LOG(3, "Indeterminate.");
    return indeterminate;
}


bool
Verifier::
_safety_proving_once(
        SafetyVerificationInput& verInput,
		const unsigned int& accuracy,
		const RealParameterSet& params) const
{
	bool result;

	SystemType& sys = verInput.getSystem();
	const HybridImageSet& initial_set = verInput.getInitialSet();
	const HybridConstraintSet& safety_constraint = verInput.getSafetyConstraint();

	RealParameterSet original_params = sys.parameters();

	sys.substitute_all(params,_settings->use_param_midpoints_for_proving);

    ARIADNE_LOG(4,"Creating the analyser...");

    const bool EQUAL_GRID_FOR_ALL_LOCATIONS = false;
    AnalyserPtrType analyser = _get_tuned_analyser(verInput,parameters_identifiers(params),
            _safety_coarse_outer_approximation->get(),safety_constraint,accuracy,EQUAL_GRID_FOR_ALL_LOCATIONS,UPPER_SEMANTICS);

	ARIADNE_LOG(5, "Using reachability restriction: " << tribool_pretty_print(!_safety_reachability_restriction.empty()));

	ARIADNE_LOG(4,"Performing outer reachability analysis...");

	try {

		ARIADNE_LOG(5, "Parameters: " << sys.parameters());

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

			HybridGridTreeSet backward_initial = analyser->initial_cell_set(initial_set,safety_constraint,
			        _safety_reachability_restriction,DIRECTION_BACKWARD);

			if (backward_initial.empty()) {
				ARIADNE_LOG(4, "The initial cell set for backward reachability is empty.");
				sys.substitute_all(original_params);
				return true;
			} else {
	            ARIADNE_LOG(5,"Backward initial set size: " << backward_initial.size());
	            if (_settings->plot_results)
	                _plot_reach(backward_initial,"final",accuracy);
			}

			ARIADNE_LOG(5,"Retrieving backward reachability...");

			HybridGridTreeSet backward_reach = analyser->outer_chain_reach(backward_initial,
					_safety_reachability_restriction,DIRECTION_BACKWARD);

			_safety_reachability_restriction = backward_reach;

			ARIADNE_LOG(6,"Reachability size: " << backward_reach.size());

			if (_settings->plot_results)
				_plot_reach(backward_reach,"backward",accuracy);

		}

		HybridGridTreeSet forward_initial;
		if (!_safety_reachability_restriction.empty()) {
            forward_initial = analyser->initial_cell_set(initial_set,safety_constraint,
                    _safety_reachability_restriction,DIRECTION_FORWARD);
		} else {
			forward_initial.adjoin_outer_approximation(initial_set,accuracy);
			forward_initial.mince(accuracy);
		}

		if (forward_initial.empty()) {
			ARIADNE_LOG(4, "The initial cell set for forward reachability is empty.");
			sys.substitute_all(original_params);
			return true;
		} else {
            ARIADNE_LOG(5,"Forward initial set size: " << forward_initial.size());
            if (_settings->plot_results)
                _plot_reach(forward_initial,"initial",accuracy);
        }

		ARIADNE_LOG(5,"Retrieving forward reachability...");

		HybridGridTreeSet forward_reach = analyser->outer_chain_reach(forward_initial,
				_safety_reachability_restriction,DIRECTION_FORWARD);

		ARIADNE_LOG(6,"Reachability size: " << forward_reach.size());

		if (_settings->plot_results)
			_plot_reach(forward_reach,"forward",accuracy);

		_update_safety_cached_reachability_with(forward_reach);

		// Additional simplified check, also useful when only forward reachability is available
		result = definitely(covers(safety_constraint,forward_reach));

	} catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(4, "The outer reached region is partially out of the domain (skipped).");
		result = false;
	} catch (ReachUnsatisfiesConstraintException ex) {
		ARIADNE_LOG(4, "The outer reached region is not inside the safe region (skipped).");
		result = false;
	}

	sys.substitute_all(original_params);

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
        SafetyVerificationInput& verInput,
        const unsigned int& accuracy,
		const RealParameterSet& params) const
{
	bool result = false;

    SystemType& sys = verInput.getSystem();
    const HybridImageSet& initial_set = verInput.getInitialSet();
    const HybridConstraintSet& safety_constraint = verInput.getSafetyConstraint();

	RealParameterSet original_params = sys.parameters();

	sys.substitute_all(params,_settings->use_param_midpoints_for_disproving);

    ARIADNE_LOG(4,"Creating the analyser...");

    const bool EQUAL_GRID_FOR_ALL_LOCATIONS = false;
    AnalyserPtrType analyser = _get_tuned_analyser(verInput,parameters_identifiers(params),
            _safety_coarse_outer_approximation->get(),safety_constraint,accuracy,EQUAL_GRID_FOR_ALL_LOCATIONS,LOWER_SEMANTICS);

    ARIADNE_LOG(5, "Using reachability restriction: " << tribool_pretty_print(!_safety_reachability_restriction.empty()));

	ARIADNE_LOG(4,"Performing lower reachability analysis...");

	try {
		HybridGridTreeSet reach;
		HybridFloatVector epsilon;
		make_lpair<HybridGridTreeSet,HybridFloatVector>(reach,epsilon) = analyser->lower_chain_reach_and_epsilon(
				initial_set,_safety_reachability_restriction);

		ARIADNE_LOG(5, "Epsilon: " << epsilon);

		if (_settings->plot_results)
			_plot_reach(reach,"lower",accuracy);

	} catch (ReachUnsatisfiesConstraintException ex) {
		ARIADNE_LOG(5, "The lower reached region is partially outside the safe region (skipped).");
		result = true;
	}

	sys.substitute_all(original_params);

	ARIADNE_LOG(4, (result ? "Disproved." : "Not disproved.") );

	return result;
}


std::list<ParametricOutcome>
Verifier::
parametric_safety(
		SafetyVerificationInput& verInput,
		const RealParameterSet& params) const
{
	ARIADNE_ASSERT_MSG(params.size() > 0, "Provide at least one parameter.");

	std::list<ParametricOutcome> result;

	std::list<RealParameterSet> splittings = maximally_split_parameters(params,_settings->maximum_parameter_depth);
	uint i=0;
	for (std::list<RealParameterSet>::const_iterator splitting_it = splittings.begin();
													 splitting_it != splittings.end();
													 ++splitting_it)
	{
	    RealParameterSet current_params = *splitting_it;
		ARIADNE_LOG(1,"Split parameters set #" << ++i << "/" << splittings.size() << ": values " << current_params);
		tribool outcome = _safety_nosplitting(verInput,current_params);
		ARIADNE_LOG(1,"Outcome: " << tribool_pretty_print(outcome));
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
	const RealParameterSet dominatingConstants;
	return _dominance(dominating,dominated,dominatingConstants);
}

tribool
Verifier::
_dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealParameter& param) const
{
	SystemType& system = dominating.getSystem();

	Real original_value = system.parameter_value(param.name());

	system.substitute(param);
	tribool result = dominance(dominating,dominated);
	system.substitute(RealParameter(param.name(),original_value));

	return result;
}

tribool
Verifier::
_dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealParameter& param,
		const Float& value) const
{
	const RealParameter modifiedConstant(param.name(),Interval(value));

	return _dominance(dominating,dominated,modifiedConstant);
}


std::list<ParametricOutcome>
Verifier::
parametric_dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealParameterSet& dominating_params) const
{
	ARIADNE_ASSERT_MSG(dominating_params.size() > 0, "Provide at least one parameter.");

	std::list<ParametricOutcome> result;

	std::list<RealParameterSet> splittings = maximally_split_parameters(dominating_params,_settings->maximum_parameter_depth);
	uint i=0;
	for (std::list<RealParameterSet>::const_iterator splitting_it = splittings.begin(); splitting_it != splittings.end(); ++splitting_it) {
	    RealParameterSet current_params = *splitting_it;
	    ARIADNE_LOG(1,"Split parameters set #" << ++i << "/" << splittings.size() << ": values " << current_params);
		tribool outcome = _dominance(dominating,dominated,current_params);
		ARIADNE_LOG(1,"Outcome: " << tribool_pretty_print(outcome));
		result.push_back(ParametricOutcome(current_params,outcome));
	}

	return result;
}


tribool
Verifier::_dominance(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealParameterSet& params) const
{
	ARIADNE_ASSERT(dominating.getProjection().size() == dominated.getProjection().size());

	ARIADNE_LOG(1, "Dominance checking...");

	if (_settings->plot_results)
		_plot_dirpath_init(dominating.getSystem().name() + "&" + dominated.getSystem().name());

	_reset_dominance_state();

	time_t initial_time = time(NULL);
    for (int accuracy = 0; time(NULL) - initial_time < _settings->time_limit_for_outcome; ++accuracy)
	{
		ARIADNE_LOG(2, "Accuracy " << accuracy);

		if (_dominance_proving_once(dominating, dominated, params, accuracy)) {
			ARIADNE_LOG(3, "Dominates.");
			return true;
		}

		if (_dominance_disproving_once(dominating, dominated, params, accuracy)) {
			ARIADNE_LOG(3, "Does not dominate.");
			return false;
		}
    }

	ARIADNE_LOG(3, "Indeterminate.");
	return indeterminate;
}

bool
Verifier::_dominance_proving_once(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealParameterSet& params,
		const unsigned int& accuracy) const
{
	ARIADNE_LOG(3,"Proving...");

	bool result;

	const RealParameterSet& original_constants = dominating.getSystem().parameters();

	dominating.getSystem().substitute_all(params,_settings->use_param_midpoints_for_proving);

	try {
		GridTreeSet flattened_dominated_lower_reach;
		Vector<Float> flattened_epsilon;
		make_lpair<GridTreeSet,Vector<Float> >(flattened_dominated_lower_reach,flattened_epsilon) =
				_dominance_flattened_lower_reach_and_epsilon(dominated,params,DOMINATED_SYSTEM,accuracy);

		GridTreeSet flattened_dominating_outer_reach = _dominance_flattened_outer_reach(dominating,params,DOMINATING_SYSTEM,accuracy);

		result = definitely(covers(flattened_dominated_lower_reach,flattened_dominating_outer_reach,flattened_epsilon));

		ARIADNE_LOG(4, "The outer reached region of the dominating system is " << (!result ? "not ": "") <<
					   "inside the projected shrinked lower reached region of the dominated system.");

	} catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominating system is partially out of the domain (skipped).");
		result = false;
	}

	ARIADNE_LOG(3, (result ? "Proved." : "Not proved.") );

	dominating.getSystem().substitute_all(original_constants);

	return result;
}

bool
Verifier::_dominance_disproving_once(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealParameterSet& params,
		const unsigned int& accuracy) const
{
	ARIADNE_LOG(3,"Disproving...");

	bool result;

	const RealParameterSet& original_params = dominating.getSystem().parameters();

	dominating.getSystem().substitute_all(params,_settings->use_param_midpoints_for_disproving);

	try {
		GridTreeSet flattened_dominating_lower_reach;
		Vector<Float> flattened_epsilon;
		make_lpair<GridTreeSet,Vector<Float> >(flattened_dominating_lower_reach,flattened_epsilon) =
				_dominance_flattened_lower_reach_and_epsilon(dominating,params,DOMINATING_SYSTEM,accuracy);

		GridTreeSet flattened_dominated_outer_reach = _dominance_flattened_outer_reach(dominated,params,DOMINATED_SYSTEM,accuracy);

		result = definitely(!inside(flattened_dominating_lower_reach,flattened_dominated_outer_reach,
				flattened_epsilon,accuracy));

	} catch (ReachOutOfDomainException ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominated system is partially out of the domain (skipped).");
		result = false;
	}

	ARIADNE_LOG(3, (result ? "Disproved." : "Not disproved.") );

	dominating.getSystem().substitute_all(original_params);

	return result;
}


std::pair<GridTreeSet,Vector<Float> >
Verifier::
_dominance_flattened_lower_reach_and_epsilon(
		DominanceVerificationInput& verInput,
		const RealParameterSet& params,
		DominanceSystem dominanceSystem,
		const unsigned int& accuracy) const
{
	string descriptor = (dominanceSystem == DOMINATING_SYSTEM ? "dominating" : "dominated");
	HybridGridTreeSet outer_approximation = (dominanceSystem == DOMINATING_SYSTEM ?
			_dominating_coarse_outer_approximation->get() : _dominated_coarse_outer_approximation->get());
	HybridGridTreeSet reachability_restriction = (dominanceSystem == DOMINATING_SYSTEM ?
			_dominating_reachability_restriction : _dominated_reachability_restriction);
	Set<Identifier> locked_params_ids = (dominanceSystem == DOMINATING_SYSTEM ? parameters_identifiers(params) : Set<Identifier>());
	HybridConstraintSet dominance_constraint; // No constraint is enforceable

	ARIADNE_LOG(4,"Creating the analyser for the " << descriptor << " system...");

    const bool EQUAL_GRID_FOR_ALL_LOCATIONS = true;
    AnalyserPtrType analyser = _get_tuned_analyser(verInput,locked_params_ids,
            outer_approximation,dominance_constraint,accuracy,EQUAL_GRID_FOR_ALL_LOCATIONS,LOWER_SEMANTICS);

	ARIADNE_LOG(5, "Using reachability restriction: " << tribool_pretty_print(!reachability_restriction.empty()));

	ARIADNE_LOG(4,"Getting its lower reached region...");

	const Vector<uint>& projection = verInput.getProjection();

	HybridGridTreeSet reach;
	HybridFloatVector epsilon;
	make_lpair<HybridGridTreeSet,HybridFloatVector>(reach,epsilon) = analyser->lower_chain_reach_and_epsilon(
			verInput.getInitialSet(),reachability_restriction);

    if (_settings->plot_results)
        _plot_dominance(reach,dominanceSystem,accuracy,LOWER_SEMANTICS);

	GridTreeSet flattened_reach = flatten_and_project_down(reach,projection);

	HybridFloatVector::const_iterator epsilon_it = epsilon.begin();
	Vector<Float> flattened_epsilon(projection.size());
	for (HybridFloatVector::const_iterator epsilon_it = epsilon.begin(); epsilon_it != epsilon.end(); ++epsilon_it) {
		for (uint i=0; i<projection.size(); ++i)
			flattened_epsilon[i] = max(flattened_epsilon[i],epsilon_it->second[projection[i]]);
	}

	ARIADNE_LOG(5,"Flattened epsilon: " << flattened_epsilon);

	return std::pair<GridTreeSet,Vector<Float> >(flattened_reach,flattened_epsilon);
}


GridTreeSet
Verifier::
_dominance_flattened_outer_reach(
		DominanceVerificationInput& verInput,
		const RealParameterSet& params,
		DominanceSystem dominanceSystem,
		const unsigned int& accuracy) const
{
	string descriptor = (dominanceSystem == DOMINATING_SYSTEM ? "dominating" : "dominated");
	OuterApproximationCache& outer_approximation_cache = (dominanceSystem == DOMINATING_SYSTEM ?
			*_dominating_coarse_outer_approximation : *_dominated_coarse_outer_approximation);
	HybridGridTreeSet& reachability_restriction = (dominanceSystem == DOMINATING_SYSTEM ?
			_dominating_reachability_restriction : _dominated_reachability_restriction);
	Set<Identifier> locked_params_ids = (dominanceSystem == DOMINATING_SYSTEM ? parameters_identifiers(params) : Set<Identifier>());
	HybridConstraintSet dominance_constraint; // No constraint is enforceable

    ARIADNE_LOG(4,"Creating the analyser for the " << descriptor << " system...");

    const bool EQUAL_GRID_FOR_ALL_LOCATIONS = true;
    AnalyserPtrType analyser = _get_tuned_analyser(verInput,locked_params_ids,
            outer_approximation_cache.get(),dominance_constraint,accuracy,EQUAL_GRID_FOR_ALL_LOCATIONS,UPPER_SEMANTICS);

    ARIADNE_LOG(5, "Using reachability restriction: " << tribool_pretty_print(!reachability_restriction.empty()));

	ARIADNE_LOG(4,"Getting its outer reached region...");

	HybridGridTreeSet reach = analyser->outer_chain_reach(verInput.getInitialSet(),reachability_restriction,DIRECTION_FORWARD);

	if (!outer_approximation_cache.is_set()) {
		outer_approximation_cache.set(reach);
	} else {
		reachability_restriction = reach;
	}

	if (_settings->plot_results)
		_plot_dominance(reach,dominanceSystem,accuracy,UPPER_SEMANTICS);

	return flatten_and_project_down(reach,verInput.getProjection());
}

Verifier::AnalyserPtrType
Verifier::
_get_tuned_analyser(
        const VerificationInput& verInput,
        const Set<Identifier>& locked_params_ids,
        const HybridGridTreeSet& outer_approx,
        const HybridConstraintSet& constraint_set,
        int accuracy,
        bool EQUAL_GRID_FOR_ALL_LOCATIONS,
        Semantics semantics) const
{
    const SystemType& sys = verInput.getSystem();
    const HybridBoxes& domain = verInput.getDomain();

    AnalyserPtrType analyser(new HybridReachabilityAnalyser(sys));

    analyser->verbosity = this->verbosity - this->child_tab_offset;
    analyser->tab_offset = this->tab_offset + this->child_tab_offset;

    analyser->tune_settings(sys,domain,locked_params_ids,outer_approx,constraint_set,
            accuracy,EQUAL_GRID_FOR_ALL_LOCATIONS,semantics);

    return analyser;
}


void
Verifier::
_reset_safety_state() const
{
	_safety_coarse_outer_approximation->reset();

	_safety_reachability_restriction = HybridGridTreeSet();
}


void
Verifier::
_reset_dominance_state() const
{
    _dominating_coarse_outer_approximation->reset();
    _dominated_coarse_outer_approximation->reset();

    _dominating_reachability_restriction = HybridGridTreeSet();
    _dominated_reachability_restriction = HybridGridTreeSet();
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
		int accuracy,
		Semantics semantics) const
{
	string system_descr = (dominanceSystem == DOMINATING_SYSTEM ? "dominating" : "dominated");
	string verification_descr = (semantics == UPPER_SEMANTICS ? "pos" : "neg");

	char mgd_char[10];
	sprintf(mgd_char,"%i",accuracy);
	string filename = system_descr + "-" + verification_descr + "-";
	filename.append(mgd_char);
	plot(_plot_dirpath,filename,reach);
}


std::list<RealParameterSet>
maximally_split_parameters(
		const RealParameterSet& params,
		const uint& maximum_parameter_depth)
{
	std::list<RealParameterSet> source;
	std::list<RealParameterSet> destination;
	destination.push_back(params);

	for (RealParameterSet::const_iterator param_it = params.begin(); param_it != params.end(); ++param_it)
	{
		for (uint i=0; i<maximum_parameter_depth; i++) {
			source.clear();
			source.insert<std::list<RealParameterSet>::const_iterator>(source.begin(),destination.begin(),destination.end());
			destination.clear();

			while (!source.empty()) {
				RealParameterSet currentParams = source.back();
				source.pop_back();

				RealParameterSet newConfigurationLeft;
				RealParameterSet newConfigurationRight;
				;
				for (RealParameterSet::const_iterator current_it = currentParams.begin(); current_it != currentParams.end(); ++current_it) {
					if (current_it->name() == param_it->name()) {
						Real currentInterval = current_it->value();
						newConfigurationLeft.insert(RealParameter(current_it->name(),
								Interval(currentInterval.lower(),currentInterval.midpoint())));
						newConfigurationRight.insert(RealParameter(current_it->name(),
								Interval(currentInterval.midpoint(),currentInterval.upper())));
					} else {
						newConfigurationLeft.insert(*current_it);
						newConfigurationRight.insert(*current_it);
					}
				}

				destination.push_back(newConfigurationLeft);
				destination.push_back(newConfigurationRight);
			}
		}
	}

	return destination;
}

}
