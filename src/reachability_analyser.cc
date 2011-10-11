/*****************************************************************************
 *            reachability_analyser.cc
 *
 *  Copyright  2006-11  Alberto Casagrande, Pieter Collins, Davide Bresolin, 
 *                      Luca Geretti
 *
 *****************************************************************************/

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
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "exceptions.h"

#include "numeric.h"

#include "vector.h"
#include "matrix.h"

#include "box.h"
#include "list_set.h"
#include "denotable_set.h"

#include "orbit.h"

#include "hybrid_time.h"
#include "hybrid_automaton.h"

#include "evolver_interface.h"

#include "reachability_analyser.h"
#include "reachability_restriction.h"
#include "logging.h"

#include "graphics.h"

#include "workers.h"


namespace Ariadne {

HybridReachabilityAnalyser::
~HybridReachabilityAnalyser()
{
}


HybridReachabilityAnalyser::
HybridReachabilityAnalyser(
		const HybridAutomatonInterface& system,
		const ReachabilityRestriction& restriction)
	: _settings(new SettingsType(system,restriction.bounding_box()))
	, _restriction(new ReachabilityRestriction(restriction))
	, _system(system.clone())
{
    this->charcode = "a";
}


HybridReachabilityAnalyser::
HybridReachabilityAnalyser(
		const HybridAutomatonInterface& system,
		const HybridBoxes& domain,
		int accuracy)
	: _system(system.clone())
{
	bool EQUAL_GRID_FOR_ALL_LOCATIONS = false;

    HybridFloatVector hmad = getHybridMidpointAbsoluteDerivatives(system,domain);
    HybridGrid grid = getHybridGrid(hmad,domain,EQUAL_GRID_FOR_ALL_LOCATIONS);

	_restriction.reset(new ReachabilityRestriction(domain,grid,accuracy));

	_settings.reset(new SettingsType(system,domain));

    this->charcode = "a";
}

int
HybridReachabilityAnalyser::
_accuracy() const
{
	return _restriction->accuracy();
}

const HybridGrid&
HybridReachabilityAnalyser::
_grid() const
{
	return _restriction->grid();
}


HybridBoxes
HybridReachabilityAnalyser::
_domain() const
{
	return _restriction->bounding_box();
}


HybridReachabilityAnalyser::EvolverPtrType
HybridReachabilityAnalyser::
_get_tuned_evolver(
        const SystemType& sys,
        uint ADD_TAB_OFFSET) const
{
    EvolverPtrType evolver(new ImageSetHybridEvolver(sys));

    evolver->verbosity = this->verbosity - ADD_TAB_OFFSET;
    evolver->tab_offset = this->tab_offset + ADD_TAB_OFFSET;

    HybridFloatVector hmad = getHybridMidpointAbsoluteDerivatives(sys,_domain());

    uint time_limit_for_result = _settings->time_limit_for_result - (time(NULL) - _start_time);

    evolver->tune_settings(_domain(),_grid(),_accuracy(),time_limit_for_result);

    return evolver;
}


void
HybridReachabilityAnalyser::
tune_settings(
        const Set<Identifier>& locked_params_ids,
        const HybridConstraintSet& constraint_set,
        uint time_limit_for_result)
{
    ARIADNE_LOG(1, "Tuning settings for analysis...");

    _settings->locked_parameters_ids = locked_params_ids;
    _settings->constraint_set = constraint_set;
    _settings->time_limit_for_result = time_limit_for_result;
}


void
HybridReachabilityAnalyser::
_plot_reach(
		const SetApproximationType& reach,
		string plot_dirpath,
		string name_prefix) const
{
	char mgd_char[10];
	sprintf(mgd_char,"%i",_accuracy());
	name_prefix.append(mgd_char);
	plot(plot_dirpath,name_prefix,reach);
}


std::pair<HybridDenotableSet,HybridDenotableSet>
HybridReachabilityAnalyser::
_upper_reach_evolve(
		const SystemType& sys,
        const HybridDenotableSet& initial_set,
        const HybridTime& time,
        bool enable_premature_termination_on_blocking_event,
        ContinuousEvolutionDirection direction) const
{
    const unsigned EVOLVER_TAB_OFFSET = 4;

	std::pair<HDS,HDS> result;
	HDS& reach = result.first;
	HDS& evolve = result.second;

	ARIADNE_LOG(3,"Evolving and discretising...");

	ARIADNE_LOG(4,"Initial size = " << initial_set.size());

	list<EnclosureType> initial_enclosures = to_enclosures(initial_set);

	const EvolverPtrType& evolver = _get_tuned_evolver(sys,EVOLVER_TAB_OFFSET);
	UpperReachEvolveWorker worker(evolver,initial_enclosures,time,
			_grid(),_accuracy(),enable_premature_termination_on_blocking_event,direction);

	result = worker.get_result();

    ARIADNE_LOG(4,"Reach size = " << reach.size());
    ARIADNE_LOG(4,"Final size = " << evolve.size());

    _check_timeout();

    return result;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_evolve(
        const HybridBoundedConstraintSet& initial_set,
        const TimeType& time) const
{
	return lower_reach_evolve(initial_set,time).second;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_reach(
        const HybridBoundedConstraintSet& initial_set,
        const TimeType& time) const
{
	return lower_reach_evolve(initial_set,time).first;
}


std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
lower_reach_evolve(
        const HybridBoundedConstraintSet& initial_set,
        const TimeType& time) const
{
	_start_time = std::time(NULL);

    const unsigned EVOLVER_TAB_OFFSET = 3;

    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach_evolve(...)");

    const HybridGrid& grid = _grid();

    HDS reach(grid); HDS evolve(grid);

    list<EnclosureType> initial_enclosures = _enclosures_from_discretised_initial_set_midpoints(initial_set);

    const EvolverPtrType& evolver = _get_tuned_evolver(*_system,EVOLVER_TAB_OFFSET);
	ARIADNE_LOG(3,"Computing evolution...");
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        Orbit<EnclosureType> orbit = evolver->orbit(*encl_it,time,LOWER_SEMANTICS);

        HDS local_reach = outer_approximation(orbit.reach(),grid,_accuracy());
        HDS local_evolve = outer_approximation(orbit.final(),grid,_accuracy());

        HybridFloatVector local_epsilon = get_epsilon(orbit.reach(),grid,_accuracy());

        HybridBoxes widened_outer_domain = widen(_restriction->outer_domain_box(),local_epsilon);
        if (!widened_outer_domain.superset(local_reach.bounding_box()))
			throw DomainException("The lower reach is not inside the over-approximated domain. Check the domain used for variables.");

        _restriction->apply_to(local_reach);

        reach.adjoin(local_reach);
        evolve.adjoin(local_evolve);

        _check_timeout();
    }

    ARIADNE_LOG(3,"Reach size = " << reach.size());
    ARIADNE_LOG(3,"Final size = " << evolve.size());

    return make_pair(reach,evolve);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_evolve(
        const HybridBoundedConstraintSet& initial_set,
        const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_evolve(set,time)");
 
	HDS reach, evolve;
	make_lpair(reach,evolve) = upper_reach_evolve(initial_set, time); // Runs the upper_reach_evolve routine on its behalf

    return evolve;
}



HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_reach(
        const HybridBoundedConstraintSet& initial_set,
        const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach(set,time)");

	HDS reach, evolve;
	make_lpair(reach,evolve) = upper_reach_evolve(initial_set, time);

    return reach;
}


std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
upper_reach_evolve(
        const HybridBoundedConstraintSet& initial_set,
        const TimeType& time) const
{
	_start_time = std::time(NULL);

    const unsigned EVOLVER_TAB_OFFSET = 4;

    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach_evolve(system,set,time)");
    ARIADNE_LOG(3,"initial_set="<<initial_set);

    HybridGrid grid = _grid();

    HDS found(grid),evolve(grid),reach(grid);

    Float real_time=time.continuous_time();
    uint discrete_steps=time.discrete_time();
    Float lock_to_grid_time=_settings->lock_to_grid_time;
    uint time_steps=uint(real_time/lock_to_grid_time);
    Float remaining_time=real_time-time_steps*lock_to_grid_time;
    if(time_steps == 0) {
        time_steps=1;
        remaining_time=0.0;
        lock_to_grid_time=real_time;
    }
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remaining_time(remaining_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time);
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time);

    ARIADNE_LOG(3,"Computing initial evolution...");
    HDS initial = _restriction->outer_intersection_with(initial_set);
    ARIADNE_LOG(4,"initial.size()="<<initial.size());

    const EvolverPtrType& evolver = _get_tuned_evolver(*_system,EVOLVER_TAB_OFFSET);
    std::list<EnclosureType> initial_enclosures = to_enclosures(initial);
    ARIADNE_LOG(4,"Starting from " << initial_enclosures.size() << " enclosures.");
    for (std::list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); ++encl_it) {
        Orbit<EnclosureType> orbit = evolver->orbit(*encl_it,hybrid_lock_to_grid_time,UPPER_SEMANTICS);
        reach.adjoin(outer_approximation(orbit.reach(),grid,_accuracy()));
        evolve.adjoin(outer_approximation(orbit.final(),grid,_accuracy()));
        _check_timeout();
    }
    ARIADNE_LOG(4,"Reach size ="<<reach.size());
    ARIADNE_LOG(4,"Final size ="<<evolve.size());

    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...");
        make_lpair(found,evolve) = _upper_reach_evolve(*_system,evolve,hybrid_lock_to_grid_time);
        reach.adjoin(found);
    }

    ARIADNE_LOG(5,"remaining_time="<<remaining_time);
    if(!evolve.empty() && remaining_time > 0) {
        ARIADNE_LOG(3,"computing evolution for the remaining time...");
        make_lpair(found,evolve) = _upper_reach_evolve(*_system,evolve,hybrid_remaining_time);
        reach.adjoin(found);
    }

    _restriction->apply_to(reach);
    _restriction->apply_to(evolve);

    ARIADNE_LOG(4,"Reach size = " << reach.size());
    ARIADNE_LOG(4,"Final size = " << evolve.size());

    return std::make_pair(reach,evolve);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
outer_chain_reach(
		const HybridBoundedConstraintSet& initial_set,
		ContinuousEvolutionDirection direction) const
{
	SetApproximationType initial = _restriction->outer_intersection_with(initial_set);

	return outer_chain_reach(initial,direction);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
outer_chain_reach(
		const SetApproximationType& initial,
		ContinuousEvolutionDirection direction) const
{
	_start_time = time(NULL);

    ARIADNE_LOG(1,"Performing outer chain reachability...");

	SetApproximationType reach;

	RealParameterSet original_parameters = _system->nonsingleton_parameters();

	ARIADNE_LOG(2,"Splitting the nonsingleton nonlocked parameter set...");

	std::list<RealParameterSet> split_parameter_set_list = _getSplitParameterSetList();

	uint i = 0;
	for (std::list<RealParameterSet>::const_iterator set_it = split_parameter_set_list.begin(); set_it != split_parameter_set_list.end(); ++set_it) {
		ARIADNE_LOG(1,"Split parameter set #" << ++i << "/" << split_parameter_set_list.size() << " : " << *set_it);

		_system->substitute_all(*set_it);

		SetApproximationType local_reach = _outer_chain_reach_splitted(*_system,initial,direction);

		reach.adjoin(local_reach);
	}

	_system->substitute_all(original_parameters);

	if (reach.empty())
		throw ReachEmptySizeException("The outer chain reachability is empty: check the initial set.");

	return reach;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
_outer_chain_reach_splitted(
		const SystemType& sys,
		const SetApproximationType& initial,
		ContinuousEvolutionDirection direction) const
{
    const Float& lock_to_grid_time = _settings->lock_to_grid_time;
    const int& lock_to_grid_steps = _settings->lock_to_grid_steps;

    HybridGrid grid = _grid();
    SetApproximationType new_final(grid), new_reach(grid), reach(grid), final(grid);

    SetApproximationType current_initial = initial;

    ARIADNE_LOG(2,"Computing recurrent " << (direction == DIRECTION_FORWARD ? "forward" : "backward") << " evolution...");
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);

    if (current_initial.empty()) {
    	ARIADNE_LOG(3,"Empty initial set, will skip calculation.");
    }

    uint i=0;
    while (!current_initial.empty())
	{
    	ARIADNE_LOG(2,"Iteration " << i++);

        static const bool ignore_activations = true;
		make_lpair(new_reach,new_final) = _upper_reach_evolve(sys,current_initial,
				hybrid_lock_to_grid_time,ignore_activations,direction);

	    _check_timeout();

        ARIADNE_LOG(3,"Removing the previously reached and final sets...");

        new_final.remove(final);
		new_reach.remove(reach);
	    ARIADNE_LOG(4,"Reach size after removal = "<<new_reach.size());
	    ARIADNE_LOG(4,"Final size after removal = "<<new_final.size());

	    ARIADNE_LOG(3,"Restricting the resulting reached and final sets...");

	    _restriction->apply_to(new_final);
	    _restriction->apply_to(new_reach);
		ARIADNE_LOG(4,"Reach size after restricting = "<<new_reach.size());
		ARIADNE_LOG(4,"Final size after restricting = "<<new_final.size());

        reach.adjoin(new_reach);
		final.adjoin(new_final);

		ARIADNE_LOG(3,"Determining the new initial set from the jump sets...");

       	current_initial = (direction == DIRECTION_FORWARD ?
       			_restriction->forward_jump_set(new_reach,sys) : _restriction->backward_jump_set(new_reach,sys));

        current_initial.adjoin(new_final);

        _check_timeout();
    }

	ARIADNE_LOG(2,"Found a total of " << reach.size() << " reached cells.");

    return reach;
}


std::pair<HybridDenotableSet,HybridFloatVector>
HybridReachabilityAnalyser::
lower_chain_reach_and_epsilon(const HybridBoundedConstraintSet& initial_set) const
{
	_start_time = time(NULL);

	HybridGrid grid = _grid();
	SetApproximationType reach(grid);
	HybridSpace state_space = _system->state_space();

    ARIADNE_LOG(1,"Performing lower chain reach with epsilon...");

	HybridFloatVector epsilon;
	for (HybridSpace::const_iterator hs_it = state_space.begin(); hs_it != state_space.end(); ++hs_it)
		epsilon.insert(std::pair<DiscreteLocation,Vector<Float> >(hs_it->first,Vector<Float>(hs_it->second)));

	RealParameterSet original_parameters = _system->nonsingleton_parameters();

	ARIADNE_LOG(2,"Splitting the nonsingleton nonlocked parameter set...");

	std::list<RealParameterSet> split_parameter_set_list = _getSplitParameterSetList();
	std::list<RealParameterSet> split_midpoint_set_list = getMidpointsSet(split_parameter_set_list);

	try {

		uint i = 0;
		for (std::list<RealParameterSet>::const_iterator set_it = split_midpoint_set_list.begin();
														set_it != split_midpoint_set_list.end();
														++set_it) {
			ARIADNE_LOG(1,"Split parameters set #" << ++i << "/" << split_midpoint_set_list.size() << " : " << *set_it);

			_system->substitute_all(*set_it);

			SetApproximationType local_reach(grid);
			HybridFloatVector local_epsilon;
			make_lpair<SetApproximationType,HybridFloatVector>(local_reach,local_epsilon) =
					_lower_chain_reach_and_epsilon(*_system,initial_set);

			reach.adjoin(local_reach);
			epsilon = max_elementwise(epsilon,local_epsilon);
		}

	} catch (ReachUnsatisfiesConstraintException& ex) {
		_system->substitute_all(original_parameters);
		throw ex;
	}

    _system->substitute_all(original_parameters);

	return std::pair<SetApproximationType,HybridFloatVector>(reach,epsilon);
}


std::pair<HybridDenotableSet,HybridFloatVector>
HybridReachabilityAnalyser::
_lower_chain_reach_and_epsilon(
		const SystemType& system,
		const HybridBoundedConstraintSet& initial_set) const
{
    const unsigned EVOLVER_TAB_OFFSET = 3;

	typedef std::list<EnclosureType> EL;
	typedef std::map<DiscreteLocation,uint> HUM;

	HybridGrid grid = _grid();
	TimeType lock_time(_settings->lock_to_grid_time,_settings->lock_to_grid_steps);

    HDS reach(grid);
	HybridSpace state_space = system.state_space();

	HybridFloatVector epsilon;
	for (HybridSpace::const_iterator hs_it = state_space.begin(); hs_it != state_space.end(); ++hs_it)
		epsilon.insert(std::pair<DiscreteLocation,Vector<Float> >(hs_it->first,Vector<Float>(hs_it->second)));

    EL initial_enclosures = _enclosures_from_discretised_initial_set_midpoints(initial_set);

    initial_enclosures = _restriction->filter(initial_enclosures);

    ARIADNE_LOG(2,"Computing recurrent evolution...");

    const EvolverPtrType& evolver = _get_tuned_evolver(system,EVOLVER_TAB_OFFSET);

    uint i=0;
    while (!initial_enclosures.empty()) {
		ARIADNE_LOG(2,"Iteration " << i++);
		ARIADNE_LOG(3,"Initial enclosures size = " << initial_enclosures.size());

		LowerReachEpsilonWorker worker(evolver,initial_enclosures,lock_time,grid,_accuracy());

		ARIADNE_LOG(3,"Evolving and discretising...");

		EL final_enclosures;
		std::pair<HUM,HUM> evolve_sizes;
		HDS local_reach;
		HybridFloatVector local_epsilon;
		make_ltuple<std::pair<HUM,HUM>,EL,HDS,HybridFloatVector>(evolve_sizes,final_enclosures,
				local_reach,local_epsilon) = worker.get_result();

		epsilon = max_elementwise(epsilon,local_epsilon);

        HybridBoxes widened_outer_domain = widen(_restriction->outer_domain_box(),local_epsilon);
        if (!widened_outer_domain.superset(local_reach.bounding_box()))
			throw DomainException("The lower reach is not inside the over-approximated domain. Check the domain used for variables.");

		_restriction->apply_to(local_reach);

		_lower_chain_reach_and_epsilon_constraint_check(system,local_reach,local_epsilon);

		ARIADNE_LOG(3,"Reach size before removal = " << local_reach.size());

		local_reach.remove(reach);

		ARIADNE_LOG(3,"Reach size after removal  = " << local_reach.size());
		if (local_reach.empty())
			break;

		reach.adjoin(local_reach);

		ARIADNE_LOG(3,"Final enclosures size = " << final_enclosures.size());

		HUM& adjoined_evolve_sizes = evolve_sizes.first;
		HUM& superposed_evolve_sizes = evolve_sizes.second;
		_filter_enclosures(final_enclosures,initial_enclosures,
				adjoined_evolve_sizes,superposed_evolve_sizes);

	    _check_timeout();
	}

	return std::pair<HDS,HybridFloatVector>(reach,epsilon);
}


void HybridReachabilityAnalyser::
_lower_chain_reach_and_epsilon_constraint_check(
		const SystemType& system,
		const HDS& reach,
		const HybridFloatVector& epsilon) const
{
	ARIADNE_LOG(3,"Checking constraint satisfaction...");

	HDS definitely_infeasible_reach = inner_difference(reach,_settings->constraint_set);

	// It is not necessary to check for eps-infeasible cells if no definitely infeasible cells are present
	if (!definitely_infeasible_reach.empty()) {

		// Identify the subspace of those locations where infeasible cells were found
		HybridSpace proj_space;
		for (HDS::locations_const_iterator loc_it = definitely_infeasible_reach.locations_begin();
				loc_it != definitely_infeasible_reach.locations_end(); ++loc_it) {
			if (!loc_it->second.empty())
				proj_space[loc_it->first] = loc_it->second.dimension();
		}

		if (_has_eps_definitely_infeasible_subset(reach,epsilon,proj_space)) {
			throw ReachUnsatisfiesConstraintException("The lower reached region partially does not satisfy the constraint.");
		} else {
			ARIADNE_LOG(3,"The reached set satisfies the constraint when epsilon is considered.");
		}

	} else {
		ARIADNE_LOG(3,"The reached set satisfies the constraint regardless of epsilon.");
	}
}


bool
HybridReachabilityAnalyser::
_has_eps_definitely_infeasible_subset(
		const HybridDenotableSet& reach,
		const HybridFloatVector& eps,
		const HybridSpace& space) const
{
	// Project the grid, constraint, residual reach and epsilon
	HybridGrid proj_grid;
	for (HybridSpace::const_iterator hs_it = space.begin(); hs_it != space.end(); ++hs_it) {
		proj_grid[hs_it->first] = _grid().find(hs_it->first)->second;
	}
	HybridConstraintSet proj_constraint_set;
	for (HybridSpace::const_iterator hs_it = space.begin(); hs_it != space.end(); ++hs_it) {
		proj_constraint_set[hs_it->first] = _settings->constraint_set.find(hs_it->first)->second;
	}
	HDS proj_reach(proj_grid);
	for (HybridSpace::const_iterator hs_it = space.begin(); hs_it != space.end(); ++hs_it) {
		proj_reach[hs_it->first] = reach.find(hs_it->first)->second;
	}
	HybridFloatVector proj_eps;
	for (HybridSpace::const_iterator hs_it = space.begin(); hs_it != space.end(); ++hs_it) {
		proj_eps[hs_it->first] = eps.find(hs_it->first)->second;
	}

	// On the projection, get the eps-enlarged constraint: this is given by the eps-enlarged codomain corresponding
	// to the image of the eps-enlarged possibly feasible set

	HybridDenotableSet proj_possibly_feasible = _restriction->outer_intersection_with(proj_constraint_set);

	HybridVectorFunction proj_constraint_functions = proj_constraint_set.functions();
	HybridBoxes proj_eps_constraint_codomain = eps_codomain(proj_possibly_feasible, proj_eps, proj_constraint_functions);
	HybridConstraintSet proj_eps_constraint(proj_constraint_functions,proj_eps_constraint_codomain);
	HybridDenotableSet eps_definitely_infeasible_proj_reach = inner_difference(proj_reach,proj_eps_constraint);
	
	return !eps_definitely_infeasible_proj_reach.empty();
}


void
HybridReachabilityAnalyser::
_filter_enclosures(
		std::list<EnclosureType>& final_enclosures,
		std::list<EnclosureType>& initial_enclosures,
		const std::map<DiscreteLocation,uint>& adjoined_evolve_sizes,
		const std::map<DiscreteLocation,uint>& superposed_evolve_sizes) const
{
	while (!final_enclosures.empty()) {
		EnclosureType encl = final_enclosures.front();
		final_enclosures.pop_front();

		const DiscreteLocation& loc = encl.location();

		/* If pruning is to be performed, push only a fraction of the final_enclosures into the initial_enclosures;
		 * otherwise, push indiscriminately.
		 */
		if (_settings->enable_lower_pruning) {
			Float coverage_ratio = (Float)adjoined_evolve_sizes.find(loc)->second/((Float)superposed_evolve_sizes.find(loc)->second);

			if (initial_enclosures.size() <= 2 || rand() < 2*coverage_ratio*RAND_MAX)
				initial_enclosures.push_back(encl);
		} else
			initial_enclosures.push_back(encl);
	}
}



std::list<RealParameterSet>
HybridReachabilityAnalyser::
_getSplitParameterSetList() const
{
    std::list<RealParameterSet> result_parameter_set_list;

    std::list<std::pair<Float,RealParameterSet> > working_scored_parameter_set_list;

    RealParameterSet initial_parameter_set = _system->nonsingleton_parameters();

    HybridFloatVector hmad = getHybridMidpointAbsoluteDerivatives(*_system,_domain());

    if (!initial_parameter_set.empty()) {

        Float initial_score = _getDerivativeWidthsScore(hmad);

        working_scored_parameter_set_list.push_back(std::pair<Float,RealParameterSet>(initial_score,initial_parameter_set));

        do {
            _updateSplitParameterSetLists(working_scored_parameter_set_list,result_parameter_set_list,
                    initial_score,hmad);
        } while (working_scored_parameter_set_list.size() > 0);

    } else {
        // We add the empty set, since not splitting or having no splittable parameters is the same
        result_parameter_set_list.push_back(initial_parameter_set);
    }

    ARIADNE_LOG(2,"Final parameter set list size: " << result_parameter_set_list.size());

	return result_parameter_set_list;
}


Float
HybridReachabilityAnalyser::
_getDerivativeWidthsScore(const HybridFloatVector& hmad) const
{
    Float result = 0;

    HybridBoxes domain = _domain();
    for (HybridBoxes::const_iterator domain_it = domain.begin(); domain_it != domain.end(); ++domain_it) {
        const DiscreteLocation& loc = domain_it->first;
        const Box& bbx = domain_it->second;

        // If the domain is empty, there is no need to calculate anything there (indeed, it would give an error)
        if (!bbx.empty()) {
			Vector<Float> der_widths = getDerivativeWidths(*_system,loc,bbx);
			Vector<Float> mad = hmad.find(loc)->second;

			for (unsigned i=0; i < bbx.size(); ++i) {
				if (mad[i] > 0) {
					result += der_widths[i]/mad[i];
				}
			}
        }
    }

    return result/domain.size();
}


std::pair<Float,Float>
HybridReachabilityAnalyser::
_getSplitDerivativeWidthsScores(
        const RealParameter& param,
        const HybridFloatVector& hmad) const
{
    const Interval& param_val = param.value();
    Float left_result = 0;
    Float right_result = 0;

    HybridBoxes domain = _domain();
    for (HybridBoxes::const_iterator domain_it = domain.begin(); domain_it != domain.end(); ++domain_it) {
        const DiscreteLocation& loc = domain_it->first;
        const Box& domain_box = domain_it->second;
        const Vector<Float>& mad = hmad.find(loc)->second;

        if (!domain_box.empty()) {
			Interval left_half_val = Interval(param_val.lower(),param_val.midpoint());
			_system->substitute(RealParameter(param.name(),left_half_val));
			Vector<Float> left_der_widths = getDerivativeWidths(*_system,loc,domain_box);
			Interval right_half_val = Interval(param_val.midpoint(),param_val.upper());
			_system->substitute(RealParameter(param.name(),right_half_val));
			Vector<Float> right_der_widths = getDerivativeWidths(*_system,loc,domain_box);
			_system->substitute(param);

			for (unsigned i=0; i < domain_box.size(); ++i) {
				if (mad[i] > 0) {
					left_result += left_der_widths[i]/mad[i];
					right_result += right_der_widths[i]/mad[i];
				}
			}
        }
    }

    return std::pair<Float,Float>(left_result/domain.size(),right_result/domain.size());
}


void
HybridReachabilityAnalyser::
_updateSplitParameterSetLists(
        std::list<std::pair<Float,RealParameterSet> >& working_scored_parameter_set_list,
        std::list<RealParameterSet>& result_parameter_set_list,
        const Float& initial_score,
        const HybridFloatVector& hmad) const
{
    std::pair<Float,RealParameterSet> working_scored_parameter_set = working_scored_parameter_set_list.back();
    working_scored_parameter_set_list.pop_back();

    const Float& previous_score = working_scored_parameter_set.first;
    const RealParameterSet& working_parameter_set = working_scored_parameter_set.second;

    RealParameterSet original_parameter_set = _system->nonsingleton_parameters();

    _system->substitute_all(working_parameter_set);

    Identifier best_param_id = "";
    Interval best_interval;
    std::pair<Float,Float> best_scores(std::numeric_limits<Float>::infinity(),std::numeric_limits<Float>::infinity());

    for (RealParameterSet::const_iterator param_it = working_parameter_set.begin();
            param_it != working_parameter_set.end(); ++param_it) {

        std::pair<Float,Float> local_scores = _getSplitDerivativeWidthsScores(*param_it,hmad);

        if (max(local_scores.first,local_scores.second) < max(best_scores.first,best_scores.second)) {
            best_scores = local_scores;
            best_param_id = param_it->name();
            best_interval = param_it->value();
        }
    }

    _system->substitute_all(original_parameter_set);

    Float ratio = (previous_score-max(best_scores.first,best_scores.second))/initial_score;

    if (ratio > _settings->splitting_parameters_target_ratio) {

        RealParameterSet new_parameters_left, new_parameters_right;
        for (RealParameterSet::const_iterator param_it = working_parameter_set.begin();
                param_it != working_parameter_set.end(); ++param_it) {
            if (param_it->name() != best_param_id) {
                new_parameters_left.insert(*param_it);
                new_parameters_right.insert(*param_it);
            } else {
                new_parameters_left.insert(RealParameter(best_param_id,Interval(best_interval.lower(),best_interval.midpoint())));
                new_parameters_right.insert(RealParameter(best_param_id,Interval(best_interval.midpoint(),best_interval.upper())));
            }
        }

        working_scored_parameter_set_list.push_front(std::pair<Float,RealParameterSet>(best_scores.first,new_parameters_left));
        working_scored_parameter_set_list.push_front(std::pair<Float,RealParameterSet>(best_scores.second,new_parameters_right));
    } else {
        result_parameter_set_list.push_back(working_parameter_set);
    }
}



HybridReachabilityAnalyserSettings::HybridReachabilityAnalyserSettings(
		const SystemType& sys,
		const HybridBoxes& domain)
    : time_limit_for_result(std::numeric_limits<uint>::max()),
      lock_to_grid_time(getLockToGridTime(sys,domain)),
      lock_to_grid_steps(1),
      constraint_set(sys.state_space()),
      splitting_parameters_target_ratio(0.05),
      enable_lower_pruning(true)
{
	HybridSpace sys_space = sys.state_space();
	ARIADNE_ASSERT_MSG(sys_space.size() == domain.size(),
			"The system and domain are of mismatched hybrid space sizes.");

	for (HybridSpace::const_iterator space_it = sys_space.begin(); space_it != sys_space.end(); ++space_it) {
		HybridBoxes::const_iterator domain_it = domain.find(space_it->first);
		ARIADNE_ASSERT_MSG(domain_it != domain.end(),
				"The domain does not have the " << space_it->first.name() << " location.");
		ARIADNE_ASSERT_MSG(domain_it->second.dimension() == space_it->second,
				"The domain and system space do not match dimensions in the " << space_it->first.name() << " location.");
	}
}


std::ostream&
operator<<(std::ostream& os, const HybridReachabilityAnalyserSettings& s)
{
    os << "DiscreteEvolutionSettings"
       << ",\n  time_limit_for_result=" << s.time_limit_for_result
       << ",\n  lock_to_grid_time=" << s.lock_to_grid_time
       << ",\n  lock_to_grid_steps=" << s.lock_to_grid_steps
       << ",\n  splitting_constants_target_ratio=" << s.splitting_parameters_target_ratio
       << ",\n  enable_lower_pruning=" << s.enable_lower_pruning
       << "\n)\n";
    return os;
}


HybridFloatVector min_cell_widths(
		const HybridGrid& grid,
		int maximum_grid_depth)
{
    const Float divider = 1 << maximum_grid_depth;

	HybridFloatVector result;

	for (HybridGrid::const_iterator grid_it = grid.begin(); grid_it != grid.end(); ++grid_it)
		result.insert(std::pair<DiscreteLocation,Vector<Float> >(grid_it->first,grid_it->second.lengths()/divider));

	return result;
}

void
remove_nonlocked_parameters(
        RealParameterSet& params,
        const Set<Identifier>& locked_params_ids)
{
    for (Set<Identifier>::const_iterator locked_parameter_it = locked_params_ids.begin();
                                         locked_parameter_it != locked_params_ids.end();
                                         ++locked_parameter_it) {
        for (RealParameterSet::iterator param_it = params.begin();
                                             param_it != params.end();
                                             ++param_it) {
            if (*locked_parameter_it == param_it->name()) {
                params.erase(param_it);
                break;
            }
        }
    }
}


HybridFloatVector
getHybridDerivativeWidths(
		const HybridReachabilityAnalyser::SystemType& system,
		const HybridBoxes& bounding_domain)
{
	const HybridSpace hspace = system.state_space();

	ARIADNE_ASSERT_MSG(bounding_domain.size() == hspace.size(), "The bounding domain must match the state space.");

	HybridFloatVector result;

	for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); hs_it++) {
		const DiscreteLocation& loc = hs_it->first;
		const Box& bx = bounding_domain.find(loc)->second;
		result.insert(pair<DiscreteLocation,Vector<Float> >(loc,getDerivativeWidths(system,loc,bx)));
	}

	return result;
}

Vector<Float>
getDerivativeWidths(
        const HybridReachabilityAnalyser::SystemType& system,
        const DiscreteLocation& loc,
        const Box& bx)
{
    const uint dim = bx.size();

    Vector<Interval> der = system.dynamic_function(loc)(bx);

    Vector<Float> result(dim);
    for (uint i=0;i<dim;i++)
        result[i] = der[i].width();

    return result;
}


std::list<RealParameterSet>
getMidpointsSet(const std::list<RealParameterSet>& intervals_set)
{
	std::list<RealParameterSet> result;

	for (std::list<RealParameterSet>::const_iterator set_it = intervals_set.begin(); set_it != intervals_set.end(); ++set_it) {
		RealParameterSet midpoints;
		for (RealParameterSet::const_iterator param_it = set_it->begin(); param_it != set_it->end(); ++param_it) {
			midpoints.insert(RealParameter(param_it->name(),param_it->value().midpoint()));
		}
		result.push_back(midpoints);
	}

	return result;
}


HybridGrid
getHybridGrid(
		const HybridFloatVector& hmad,
		const HybridBoxes& domain,
		bool equal_for_all_locations)
{
	return (equal_for_all_locations ? HybridGrid(HybridSpace(domain),getGrid(hmad,domain)) :
			getHybridGrid(hmad,domain));
}


HybridGrid
getHybridGrid(
		const HybridFloatVector& hmad,
		const HybridBoxes& domain)
{
	// Get the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	HybridGrid hg;

	// The lengths of the grid cell
	std::map<DiscreteLocation,Vector<Float> > hybridgridlengths;

	// Get the minimum domain length for each variable
	Vector<Float> minDomainLengths(css,std::numeric_limits<double>::infinity());
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		minDomainLengths = min_elementwise(minDomainLengths,domain.find(hfv_it->first)->second.widths());
	}

	// Initialize the minimum non-zero length for each variable as the minimum domain lengths
	Vector<Float> minNonZeroLengths = minDomainLengths;

	// Get the maximum derivative/domainlength ratio among variables and locations
	Float maxratio = 0.0;
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		for (uint i=0;i<css;i++) {
			maxratio = max(maxratio,hfv_it->second[i]/minDomainLengths[i]);
		}
	}

	// Update the minimum non zero lengths
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		for (uint i=0;i<css;i++) {
			if (hfv_it->second[i] > 0) {
				minNonZeroLengths[i] = min(minNonZeroLengths[i],hfv_it->second[i]/maxratio);
			}
		}
	}

	// Get the grid lengths from the derivatives
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
	{
		// Initialize the lengths
		Vector<Float> gridlengths(css);

		// Assign the lengths
		for (uint i=0;i<css;i++) {
			gridlengths[i] = (hfv_it->second[i] > 0) ? hfv_it->second[i]/maxratio : minNonZeroLengths[i];
		}

		// Add the pair to the hybrid lengths
		hybridgridlengths.insert(make_pair<DiscreteLocation,Vector<Float> >(hfv_it->first,gridlengths));
	}

	// Populate the grid, centered on the centre of the domain
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		const DiscreteLocation& loc = hfv_it->first;
		hg[loc] = Grid(domain.find(loc)->second.centre(),hybridgridlengths[loc]);
	}

	return hg;
}


Grid
getGrid(
		const HybridFloatVector& hmad,
		const HybridBoxes& domain)
{
	// Get the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	// The lengths of the grid cell
	Vector<Float> gridlengths(css,std::numeric_limits<double>::infinity());

	// Get the minimum domain length for each variable
	Vector<Float> minDomainLengths(css,std::numeric_limits<double>::infinity());
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		minDomainLengths = min_elementwise(minDomainLengths,domain.find(hfv_it->first)->second.widths());
	}

	// Initialize the minimum non-zero length for each variable as the minimum domain lengths
	Vector<Float> minNonZeroLengths = minDomainLengths;

	// Get the maximum derivative/domainlength ratio among variables and locations
	Float maxratio = 0.0;
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		for (uint i=0;i<css;i++) {
			maxratio = max(maxratio,hfv_it->second[i]/minDomainLengths[i]);
		}
	}

	// Update the minimum non zero lengths
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		for (uint i=0;i<css;i++) {
			if (hfv_it->second[i] > 0) {
				minNonZeroLengths[i] = min(minNonZeroLengths[i],hfv_it->second[i]/maxratio);
			}
		}
	}

	// Get the grid lengths from the derivatives
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		// Assign the lengths
		for (uint i=0;i<css;i++) {
			gridlengths[i] = min(gridlengths[i],(hfv_it->second[i] > 0) ? hfv_it->second[i]/maxratio : minNonZeroLengths[i]);
		}
	}

	return Grid(gridlengths);
}


Float
getLockToGridTime(
		const HybridReachabilityAnalyser::SystemType& sys,
		const HybridBoxes& domain)
{
	Float result = 0;

	const HybridSpace hspace = sys.state_space();
	for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); hs_it++) {

		const DiscreteLocation& loc = hs_it->first;
		const uint dim = hs_it->second;
		const Box& loc_domain = domain.find(loc)->second;

		// The domain can be empty, when constructed from an empty reachability
		if (!loc_domain.empty()) {

			// Gets the first order derivatives in respect to the dynamic of the location, applied to the corresponding domain
			Vector<Interval> der_bbx = sys.dynamic_function(loc)(loc_domain);

			for (uint i=0;i<dim;i++) {
				Float maxAbsDer = abs(der_bbx[i]).upper();
				if (maxAbsDer > 0)
					result = max(result,loc_domain[i].width()/maxAbsDer);
			}
		}
	}

	return result;
}


HybridFloatVector
getHybridMaximumAbsoluteDerivatives(
		const HybridReachabilityAnalyser::SystemType& sys,
		const HybridDenotableSetPtr& outer_approximation,
		const HybridBoxes& domain_constraint)
{
	HybridFloatVector result;

	const HybridSpace hspace = sys.state_space();
	for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); hs_it++) {

		const DiscreteLocation& loc = hs_it->first;
		const uint dim = hs_it->second;

		Vector<Interval> der_bbx;

		// Insert the corresponding hmad pair, initialized with zero maximum absolute derivatives
		result.insert(pair<DiscreteLocation,Vector<Float> >(loc,Vector<Float>(dim)));

		// If the reached region for the location exists and is not empty, check its cells, otherwise use the whole domain
		if (outer_approximation && outer_approximation->has_location(loc) && !(*outer_approximation)[loc].empty()) {

			DenotableSetType reach = (*outer_approximation)[loc];
			for (DenotableSetType::const_iterator cells_it = reach.begin(); cells_it != reach.end(); cells_it++) {
				// Gets the derivative bounds
				der_bbx = sys.dynamic_function(loc)(cells_it->box());

				// For each variable, sets the maximum value
				for (uint i=0;i<dim;i++)
					result[loc][i] = max(result[loc][i], abs(der_bbx[i]).upper());
			}
		} else {
			// Gets the first order derivatives in respect to the dynamic of the mode, applied to the domain of the corresponding location
			der_bbx = sys.dynamic_function(loc)(domain_constraint.find(loc)->second);

			// Gets the maximum absolute derivatives
			for (uint i=0;i<dim;i++)
				result[loc][i] = abs(der_bbx[i]).upper();
		}
	}

	return result;
}


list<LocalisedEnclosureType>
HybridReachabilityAnalyser::
_enclosures_from_discretised_initial_set_midpoints(const HybridBoundedConstraintSet initial_set) const
{
	list<EnclosureType> result;

	HybridDenotableSet discretised_initial_set = _restriction->outer_intersection_with(initial_set);

	for(HybridDenotableSet::const_iterator cell_it=discretised_initial_set.begin();
		 cell_it!=discretised_initial_set.end(); ++cell_it) {

		const DiscreteLocation& loc = cell_it->first;
		const Box& cell_bx = cell_it->second.box();
		const BoundedConstraintSet loc_initial_set = initial_set.find(loc)->second;

		Box cell_domain_intersection = intersection(cell_bx,loc_initial_set.domain());
		Box midpoint_box(cell_domain_intersection.centre());

		if (loc_initial_set.unconstrained() || definitely(loc_initial_set.covers(midpoint_box)))
			result.push_back(EnclosureType(loc,ContinuousEnclosureType(midpoint_box)));
	}

	return result;
}


void
HybridReachabilityAnalyser::
_check_timeout() const
{
	time_t current_time = time(NULL);

	if (current_time - _start_time > _settings->time_limit_for_result)
		throw TimeoutException();
}


std::list<LocalisedEnclosureType>
to_enclosures(const HybridDenotableSet& reach)
{
	std::list<LocalisedEnclosureType> result;

	for (HybridDenotableSet::const_iterator cell_it = reach.begin(); cell_it != reach.end(); ++cell_it)
		result.push_back(LocalisedEnclosureType(cell_it->first,cell_it->second));

	return result;
}

HybridFloatVector
get_epsilon(
		const ListSet<LocalisedEnclosureType>& enclosures,
		const HybridGrid& grid,
		int accuracy) {

	HybridFloatVector result;
	HybridSpace state_space = grid.state_space();
	HybridFloatVector grid_lengths = grid.lengths();
	Float accuracy_divider = (1 << accuracy);
	for (HybridSpace::const_iterator hs_it = state_space.begin(); hs_it != state_space.end(); ++hs_it)
		result.insert(std::pair<DiscreteLocation,Vector<Float> >(hs_it->first,Vector<Float>(hs_it->second)));

	for (ListSet<LocalisedEnclosureType>::const_iterator encl_it = enclosures.begin(); encl_it != enclosures.end(); ++encl_it) {
		result[encl_it->first] = max_elementwise(result[encl_it->first],encl_it->second.bounding_box().widths());
	}
	for (HybridSpace::const_iterator hs_it = state_space.begin(); hs_it != state_space.end(); ++hs_it)
		result[hs_it->first] += grid_lengths[hs_it->first]/accuracy_divider;

	return result;
}


} // namespace Ariadne
