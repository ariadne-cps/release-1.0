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
#include <string>

#include "verifier.h"
#include "logging.h"
#include "interruptible.h"

namespace Ariadne {


Verifier::Verifier()
    : _settings(new VerifierSettings())
{
    this->charcode = "v";
}


Verifier::~Verifier()
{
}


tribool
Verifier::
safety(SafetyVerificationInput& verInput) const
{
	RealParameterSet params;

	tribool result = _safety_nosplitting(verInput,params);

    ARIADNE_LOG(1,"Outcome: " << tribool_pretty_print(result));

    return result;
}


tribool
Verifier::
_safety_nosplitting(
		SafetyVerificationInput& verInput,
		const RealParameterSet& params) const
{
	ARIADNE_LOG(2,"Safety checking...");

	_reset_start_time();

	SystemType& system = verInput.getSystem();

	if (_settings->plot_results)
		_plotter_ptr.reset(new PlotHelper(system.name()));

	_init_safety_restriction(verInput);

	try {
		int accuracy = 0;
		while (true) {

			ARIADNE_LOG(2, "Accuracy " << accuracy);

			_safety_restriction->refine_at(accuracy);

			tribool result = _safety_once(verInput,params);

			if (!indeterminate(result))
				return result;

			++accuracy;
		}
	} catch (TimeoutException& ex) {
		ARIADNE_LOG(2, "Time available for verification expired.");
    }

	return indeterminate;
}


tribool
Verifier::
_safety_once(SafetyVerificationInput& verInput,
		const RealParameterSet& params) const
{
	try {
		if (_safety_proving_once(verInput,params)) {
			ARIADNE_LOG(3, "Safe.");
			return true;
		}
	} catch (SqrtNumericException& ex) {
		ARIADNE_LOG(3, "WARNING: caught " << ex.what() << " exception, skipping.");
	}
	try {
		if (_safety_disproving_once(verInput,params)) {
			ARIADNE_LOG(3, "Unsafe.");
			return false;
		}
	} catch (SqrtNumericException& ex) {
		ARIADNE_LOG(3, "WARNING: caught " << ex.what() << " exception, skipping.");
	}

    ARIADNE_LOG(3, "Indeterminate.");
    return indeterminate;
}


bool
Verifier::
_safety_proving_once(
        SafetyVerificationInput& verInput,
		const RealParameterSet& params) const
{
	bool result;

    ARIADNE_LOG(3, "Proving...");

	SystemType& sys = verInput.getSystem();

	RealParameterSet original_params = sys.parameters();

	sys.substitute_all(params,_settings->use_param_midpoints_for_proving);

	ARIADNE_LOG(4,"Performing outer reachability analysis...");

	try {

	    ARIADNE_LOG(5,"Computing the initial set for forward reachability...");

	    SetApproximationType forward_initial = _safety_restriction->outer_intersection_with(verInput.getInitialSet());

	    if (forward_initial.empty()) {
	        throw EmptyInitialCellSetException("The initial cell set for forward reachability is empty (skipped).");
	    } else {
	        ARIADNE_LOG(6,"Forward initial set size: " << forward_initial.size());
	        if (_settings->plot_results)
	            _plotter_ptr->plot(forward_initial,"initial",_safety_restriction->accuracy());
	    }

		result = _safety_proving_once_forward_analysis(verInput.getSystem(),forward_initial,verInput.getSafetyConstraint(),params);

		if (!indeterminate(result) && _settings->enable_backward_refinement_for_safety_proving)
		    _safety_proving_once_backward_refinement(verInput.getSystem(),forward_initial,verInput.getSafetyConstraint(),params);

	} catch (ReachOutOfDomainException& ex) {
		ARIADNE_LOG(5, "The outer reached region is partially out of the domain (skipped).");
		result = false;
	} catch (ReachUnsatisfiesConstraintException& ex) {
		ARIADNE_LOG(5, "The outer reached region is not inside the safe region (skipped).");
		result = false;
	} catch (EmptyInitialCellSetException& ex) {
        ARIADNE_LOG(5, ex.what());
        result = true;
    }

	sys.substitute_all(original_params);

    ARIADNE_LOG(4, (result ? "Proved." : "Not proved.") );

	return result;
}


bool
Verifier::
_safety_proving_once_forward_analysis(
		const SystemType& sys,
		const SetApproximationType& initial_set,
		const HybridConstraintSet& safety_constraint,
        const RealParameterSet& params) const
{
    const unsigned ANALYSER_TAB_OFFSET = 5;

    ARIADNE_LOG(5,"Creating the analyser for forward reachability...");

    AnalyserPtrType analyser = _get_tuned_analyser(sys,parameters_identifiers(params),_safety_restriction,
    		safety_constraint,ANALYSER_TAB_OFFSET);

    ARIADNE_LOG(5,"Retrieving forward reachability...");

    SetApproximationType forward_reach = analyser->outer_chain_reach(initial_set,DIRECTION_FORWARD);

    ARIADNE_LOG(6,"Reachability size: " << forward_reach.size());

    if (_settings->plot_results)
        _plotter_ptr->plot(forward_reach,"forward",_safety_restriction->accuracy());

    _safety_restriction->update_with(forward_reach);

    return definitely(covers(safety_constraint,forward_reach));
}



void
Verifier::
_safety_proving_once_backward_refinement(
		const SystemType& sys,
		const SetApproximationType& initial_set,
		const HybridConstraintSet& safety_constraint,
        const RealParameterSet& params) const
{
    const unsigned ANALYSER_TAB_OFFSET = 5;

    ARIADNE_LOG(5,"Computing the initial set for backward reachability...");

    SetApproximationType backward_initial = _safety_restriction->outer_difference_from(safety_constraint);

    if (backward_initial.empty()) {
        throw EmptyInitialCellSetException("The initial cell set for backward reachability is empty (skipped).");
    } else {
        ARIADNE_LOG(6,"Backward initial set size: " << backward_initial.size());
        if (_settings->plot_results)
            _plotter_ptr->plot(backward_initial,"final",_safety_restriction->accuracy());
    }

    ARIADNE_LOG(5,"Creating the analyser...");

    AnalyserPtrType analyser = _get_tuned_analyser(sys,parameters_identifiers(params),_safety_restriction,
    		safety_constraint,ANALYSER_TAB_OFFSET);

    ARIADNE_LOG(5,"Retrieving backward reachability...");

    SetApproximationType backward_reach = analyser->outer_chain_reach(backward_initial,DIRECTION_BACKWARD);

    _safety_restriction->update_with(backward_reach);

    ARIADNE_LOG(6,"Reachability size: " << backward_reach.size());

    SetApproximationType reached_initial = initial_set;
    _safety_restriction->apply_to(reached_initial);
    if (reached_initial.empty()) {
        throw EmptyInitialCellSetException("The initial set is now outside the reachability restriction (skipped).");
    }

    if (_settings->plot_results)
        _plotter_ptr->plot(backward_reach,"backward",_safety_restriction->accuracy());
}


bool
Verifier::
_safety_disproving_once(
        SafetyVerificationInput& verInput,
		const RealParameterSet& params) const
{
    const unsigned ANALYSER_TAB_OFFSET = 5;

    ARIADNE_LOG(3, "Disproving...");

    SystemType& sys = verInput.getSystem();
    const HybridBoundedConstraintSet& initial_set = verInput.getInitialSet();
    const HybridConstraintSet& safety_constraint = verInput.getSafetyConstraint();

	RealParameterSet original_params = sys.parameters();

	sys.substitute_all(params,_settings->use_param_midpoints_for_disproving);

	ARIADNE_LOG(4,"Performing lower reachability analysis...");

    ARIADNE_LOG(5,"Creating the analyser for forward reachability...");

    AnalyserPtrType analyser = _get_tuned_analyser(verInput.getSystem(),parameters_identifiers(params),_safety_restriction,
            safety_constraint,ANALYSER_TAB_OFFSET);

	bool result = false;
	try {

	    ARIADNE_LOG(5,"Retrieving forward reachability...");

		SetApproximationType reach;
		HybridFloatVector epsilon;
		make_lpair<SetApproximationType,HybridFloatVector>(reach,epsilon) = analyser->lower_chain_reach_and_epsilon(initial_set);

		if (_settings->plot_results)
			_plotter_ptr->plot(reach,"lower",_safety_restriction->accuracy());

	} catch (ReachUnsatisfiesConstraintException& ex) {
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

	ARIADNE_LOG(2, "Dominance checking...");

	_reset_start_time();

	if (_settings->plot_results)
		_plotter_ptr.reset(new PlotHelper(dominating.getSystem().name() + "+" + dominated.getSystem().name()));

	_init_dominance_restriction(dominating,DOMINATING_SYSTEM);
	_init_dominance_restriction(dominated,DOMINATED_SYSTEM);

	try {
		int accuracy = 0;

		while (true) {
			ARIADNE_LOG(2, "Accuracy " << accuracy);

			_dominating_restriction->refine_at(accuracy);
			_dominated_restriction->refine_at(accuracy);

			try {
				if (_dominance_proving_once(dominating, dominated, params)) {
					ARIADNE_LOG(2, "Dominates.");
					return true;
				}
			} catch (SqrtNumericException& ex) {
				ARIADNE_LOG(3, "WARNING: caught " << ex.what() << ", skipping.");
			}

			try {
				if (_dominance_disproving_once(dominating, dominated, params)) {
					ARIADNE_LOG(2, "Does not dominate.");
					return false;
				}
			} catch (SqrtNumericException& ex) {
				ARIADNE_LOG(3, "WARNING: caught " << ex.what() << ", skipping.");
			}

			++accuracy;
		}
    } catch (TimeoutException& ex) {
		ARIADNE_LOG(2, "Time available for verification expired.");
    }

	ARIADNE_LOG(2, "Indeterminate.");
	return indeterminate;
}

bool
Verifier::_dominance_proving_once(
		DominanceVerificationInput& dominating,
		DominanceVerificationInput& dominated,
		const RealParameterSet& params) const
{
	ARIADNE_LOG(3,"Proving...");

	bool result;

	const RealParameterSet& original_constants = dominating.getSystem().parameters();

	dominating.getSystem().substitute_all(params,_settings->use_param_midpoints_for_proving);

	try {
		DenotableSetType flattened_dominated_lower_reach;
		Vector<Float> flattened_epsilon;
		make_lpair<DenotableSetType,Vector<Float> >(flattened_dominated_lower_reach,flattened_epsilon) =
				_dominance_flattened_lower_reach_and_epsilon(dominated,params,DOMINATED_SYSTEM);

		DenotableSetType flattened_dominating_outer_reach = _dominance_flattened_outer_reach(dominating,params,DOMINATING_SYSTEM);

		result = definitely(covers(flattened_dominated_lower_reach,flattened_dominating_outer_reach,flattened_epsilon));

		ARIADNE_LOG(4, "The outer reached region of the dominating system is " << (!result ? "not ": "") <<
					   "inside the projected shrinked lower reached region of the dominated system.");

	} catch (ReachOutOfDomainException& ex) {
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
		const RealParameterSet& params) const
{
	ARIADNE_LOG(3,"Disproving...");

	bool result;

	const RealParameterSet& original_params = dominating.getSystem().parameters();

	dominating.getSystem().substitute_all(params,_settings->use_param_midpoints_for_disproving);

	try {
		DenotableSetType flattened_dominating_lower_reach;
		Vector<Float> flattened_epsilon;
		make_lpair<DenotableSetType,Vector<Float> >(flattened_dominating_lower_reach,flattened_epsilon) =
				_dominance_flattened_lower_reach_and_epsilon(dominating,params,DOMINATING_SYSTEM);

		DenotableSetType flattened_dominated_outer_reach = _dominance_flattened_outer_reach(dominated,params,DOMINATED_SYSTEM);

		result = definitely(!inside(flattened_dominating_lower_reach,flattened_dominated_outer_reach,
				flattened_epsilon,_dominated_restriction->accuracy()));

	} catch (ReachOutOfDomainException& ex) {
		ARIADNE_LOG(4,"The outer reached region of the dominated system is partially out of the domain (skipped).");
		result = false;
	}

	ARIADNE_LOG(3, (result ? "Disproved." : "Not disproved.") );

	dominating.getSystem().substitute_all(original_params);

	return result;
}


std::pair<DenotableSetType,Vector<Float> >
Verifier::
_dominance_flattened_lower_reach_and_epsilon(
		DominanceVerificationInput& verInput,
		const RealParameterSet& params,
		DominanceSystem dominanceSystem) const
{
    const unsigned ANALYSER_TAB_OFFSET = 4;

	string descriptor = (dominanceSystem == DOMINATING_SYSTEM ? "dominating" : "dominated");
	ReachabilityRestrictionPtr& restriction = (dominanceSystem == DOMINATING_SYSTEM ?
			_dominating_restriction : _dominated_restriction);
	Set<Identifier> locked_params_ids = (dominanceSystem == DOMINATING_SYSTEM ? parameters_identifiers(params) : Set<Identifier>());
	HybridConstraintSet dominance_constraint(restriction->space()); // No constraint is enforceable

	ARIADNE_LOG(4,"Creating the analyser for the " << descriptor << " system...");

    AnalyserPtrType analyser = _get_tuned_analyser(verInput.getSystem(),locked_params_ids,restriction,
    		dominance_constraint,ANALYSER_TAB_OFFSET);

	ARIADNE_LOG(4,"Getting its lower reached region...");

	const Vector<uint>& projection = verInput.getProjection();

	SetApproximationType reach;
	HybridFloatVector epsilon;
	make_lpair<SetApproximationType,HybridFloatVector>(reach,epsilon) = analyser->lower_chain_reach_and_epsilon(verInput.getInitialSet());

    if (_settings->plot_results)
        _plot_dominance(reach,dominanceSystem,restriction->accuracy(),LOWER_SEMANTICS);

    ARIADNE_LOG(4,"Flattening the reached region to non-hybrid and projecting to the common variables space...");

	DenotableSetType flattened_reach = flatten_and_project_down(reach,projection);

	HybridFloatVector::const_iterator epsilon_it = epsilon.begin();
	Vector<Float> flattened_epsilon(projection.size());
	for (HybridFloatVector::const_iterator epsilon_it = epsilon.begin(); epsilon_it != epsilon.end(); ++epsilon_it) {
		for (uint i=0; i<projection.size(); ++i)
			flattened_epsilon[i] = max(flattened_epsilon[i],epsilon_it->second[projection[i]]);
	}

	ARIADNE_LOG(5,"Flattened epsilon: " << flattened_epsilon);

	return std::pair<DenotableSetType,Vector<Float> >(flattened_reach,flattened_epsilon);
}


DenotableSetType
Verifier::
_dominance_flattened_outer_reach(
		DominanceVerificationInput& verInput,
		const RealParameterSet& params,
		DominanceSystem dominanceSystem) const
{
    const unsigned ANALYSER_TAB_OFFSET = 4;

	string descriptor = (dominanceSystem == DOMINATING_SYSTEM ? "dominating" : "dominated");
	ReachabilityRestrictionPtr& restriction = (dominanceSystem == DOMINATING_SYSTEM ?
			_dominating_restriction : _dominated_restriction);
	Set<Identifier> locked_params_ids = (dominanceSystem == DOMINATING_SYSTEM ? parameters_identifiers(params) : Set<Identifier>());
	HybridConstraintSet dominance_constraint(restriction->space()); // No constraint is enforceable

    ARIADNE_LOG(4,"Creating the analyser for the " << descriptor << " system...");

    AnalyserPtrType analyser = _get_tuned_analyser(verInput.getSystem(),locked_params_ids,restriction,
    		dominance_constraint,ANALYSER_TAB_OFFSET);

	ARIADNE_LOG(4,"Getting its outer reached region...");

	SetApproximationType reach = analyser->outer_chain_reach(verInput.getInitialSet(),DIRECTION_FORWARD);

	restriction->update_with(reach);

	if (_settings->plot_results)
		_plot_dominance(reach,dominanceSystem,restriction->accuracy(),UPPER_SEMANTICS);

	ARIADNE_LOG(4,"Flattening the reached region to non-hybrid and projecting to the common variables space...");

	DenotableSetType flattened_reach = flatten_and_project_down(reach,verInput.getProjection());

	return flattened_reach;
}

Verifier::AnalyserPtrType
Verifier::
_get_tuned_analyser(
        const SystemType& sys,
        const Set<Identifier>& locked_params_ids,
        const ReachabilityRestrictionPtr& restriction,
        const HybridConstraintSet& constraint_set,
        unsigned ADD_TAB_OFFSET) const
{
    AnalyserPtrType analyser(new HybridReachabilityAnalyser(sys,*restriction));

    analyser->verbosity = this->verbosity - ADD_TAB_OFFSET;
    analyser->tab_offset = this->tab_offset + ADD_TAB_OFFSET;

    analyser->ttl = this->_remaining_time();

    analyser->tune_settings(locked_params_ids,constraint_set);

    return analyser;
}


void
Verifier::
_init_safety_restriction(const SafetyVerificationInput& verInput) const
{
	const bool EQUAL_GRID_FOR_ALL_LOCATIONS = false;

	_init_restriction(verInput,_safety_restriction,EQUAL_GRID_FOR_ALL_LOCATIONS);
}



void
Verifier::
_init_dominance_restriction(
		const DominanceVerificationInput& verInput,
		DominanceSystem dominanceSystem) const
{
	const bool EQUAL_GRID_FOR_ALL_LOCATIONS = true;

	ReachabilityRestrictionPtr& restriction = (dominanceSystem == DOMINATING_SYSTEM ?
			_dominating_restriction : _dominated_restriction);

	_init_restriction(verInput,restriction,EQUAL_GRID_FOR_ALL_LOCATIONS);
}


void
Verifier::
_init_restriction(
		const VerificationInput& verInput,
		ReachabilityRestrictionPtr& restriction,
		bool EQUAL_GRID_FOR_ALL_LOCATIONS) const
{
	const int ACCURACY = 0;

	const SystemType& sys = verInput.getSystem();
    const HybridBoxes& domain = verInput.getDomain();

    HybridFloatVector hmad = getHybridMidpointAbsoluteDerivatives(sys,domain);
    HybridGrid grid = getHybridGrid(hmad,domain,EQUAL_GRID_FOR_ALL_LOCATIONS);

    restriction.reset(new ReachabilityRestriction(domain,grid,ACCURACY));
}

void
Verifier::
_plot_dominance(
		const SetApproximationType& reach,
		DominanceSystem dominanceSystem,
		int accuracy,
		Semantics semantics) const
{
	string system_descr = (dominanceSystem == DOMINATING_SYSTEM ? "dominating" : "dominated");
	string verification_descr = (semantics == UPPER_SEMANTICS ? "pos" : "neg");

	string filename = system_descr + "-" + verification_descr + "-";

    _plotter_ptr->plot(reach,filename,accuracy);
}


VerifierSettings::VerifierSettings() :
        plot_results(false),
        maximum_parameter_depth(3),
        use_param_midpoints_for_proving(false),
        use_param_midpoints_for_disproving(true),
        enable_backward_refinement_for_safety_proving(true)
{ }

std::ostream&
operator<<(std::ostream& os, const VerifierSettings& s)
{
    os << "VerificationSettings"
       << "(\n  plot_results=" << s.plot_results
       << ",\n  maximum_parameter_depth=" << s.maximum_parameter_depth
       << ",\n  use_param_midpoints_for_proving=" << s.use_param_midpoints_for_proving
       << ",\n  use_param_midpoints_for_disproving=" << s.use_param_midpoints_for_disproving
       << ",\n  enable_backward_refinement_for_safety_proving=" << s.enable_backward_refinement_for_safety_proving
       << "\n)\n";
    return os;
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
