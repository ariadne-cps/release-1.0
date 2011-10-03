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
	, free_cores(0)
{
    this->charcode = "a";
}


HybridReachabilityAnalyser::
HybridReachabilityAnalyser(
		const HybridAutomatonInterface& system,
		const HybridBoxes& domain,
		int accuracy)
	: _system(system.clone())
	, free_cores(0)
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
        const HybridAutomatonInterface& sys,
        unsigned ADD_TAB_OFFSET,
        Semantics semantics) const
{
    EvolverPtrType evolver(new ImageSetHybridEvolver(sys));

    evolver->verbosity = this->verbosity - ADD_TAB_OFFSET;
    evolver->tab_offset = this->tab_offset + ADD_TAB_OFFSET;

    HybridFloatVector hmad = getHybridMidpointAbsoluteDerivatives(sys,_domain());

    evolver->tune_settings(_grid(),hmad,_accuracy(),this->free_cores,semantics);

    return evolver;
}


void
HybridReachabilityAnalyser::
tune_settings(
        const Set<Identifier>& locked_params_ids,
        const HybridConstraintSet& constraint_set,
        unsigned free_cores,
        bool enable_lower_reach_restriction_check,
        Semantics semantics)
{
    ARIADNE_LOG(1, "Tuning settings for analysis...");

    this->free_cores = free_cores;

    _settings->locked_parameters_ids = locked_params_ids;
    _settings->constraint_set = constraint_set;
    _settings->enable_lower_reach_restriction_check = enable_lower_reach_restriction_check;
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


HybridDenotableSet
HybridReachabilityAnalyser::
initial_cells_set(const HybridImageSet& initial_enclosure_set) const
{
    SetApproximationType result(_grid());

    result.adjoin_outer_approximation(initial_enclosure_set,_accuracy());

    _restriction->apply_to(result);

    return result;
}


HybridDenotableSet
HybridReachabilityAnalyser::
initial_cells_set(const HybridConstraintSet& initial_constraint_set) const
{
    return _restriction->possibly_infeasible_projection(initial_constraint_set);
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

	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;
	ARIADNE_ASSERT_MSG(concurrency>0 && concurrency <= boost::thread::hardware_concurrency(),"Error: concurrency must be positive and less than the maximum allowed.");

	list<EnclosureType> initial_enclosures = to_enclosures(initial_set);

	const EvolverPtrType& evolver = _get_tuned_evolver(sys,EVOLVER_TAB_OFFSET,UPPER_SEMANTICS);
	UpperReachEvolveWorker worker(evolver,initial_enclosures,time,
			_grid(),_accuracy(),enable_premature_termination_on_blocking_event,direction,concurrency);
	result = worker.get_result();

    ARIADNE_LOG(4,"Reach size = " << reach.size());
    ARIADNE_LOG(4,"Final size = " << evolve.size());

    return result;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_evolve(
        const HybridImageSet& initial_set,
        const TimeType& time) const
{
	return lower_reach_evolve(initial_set,time).second;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_reach(
        const HybridImageSet& initial_set,
        const TimeType& time) const
{
	return lower_reach_evolve(initial_set,time).first;
}


std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
lower_reach_evolve(
        const HybridImageSet& initial_set,
        const TimeType& time) const
{
    const unsigned EVOLVER_TAB_OFFSET = 3;

    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach_evolve(...)");

    const HybridGrid& grid = _grid();

    HDS reach(grid); HDS evolve(grid);

    list<EnclosureType> initial_enclosures = enclosures_from_split_initial_set_midpoints(initial_set,min_cell_widths(grid,_accuracy()));

    const EvolverPtrType& evolver = _get_tuned_evolver(*_system,EVOLVER_TAB_OFFSET,LOWER_SEMANTICS);
	ARIADNE_LOG(3,"Computing evolution...");
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        Orbit<EnclosureType> orbit = evolver->orbit(*encl_it,time,LOWER_SEMANTICS);

        HDS local_reach = outer_approximation(orbit.reach(),grid,_accuracy());
        HDS local_evolve = outer_approximation(orbit.final(),grid,_accuracy());

        if (_settings->enable_lower_reach_restriction_check) {
        	ARIADNE_ASSERT_MSG(!_restriction->restricts(local_reach), "The lower reach is not inside the reachability restriction. Check the bounding domain used.");
        }

        _restriction->apply_to(local_reach);

        reach.adjoin(local_reach);
        evolve.adjoin(local_evolve);
    }

    reach.recombine();
    evolve.recombine();

    ARIADNE_LOG(3,"Reach size = " << reach.size());
    ARIADNE_LOG(3,"Final size = " << evolve.size());

    return make_pair(reach,evolve);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_evolve(
        const HybridImageSet& initial_set,
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
        const HybridImageSet& initial_set,
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
        const HybridImageSet& initial_set,
        const TimeType& time) const
{
    const unsigned EVOLVER_TAB_OFFSET = 4;

    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach_evolve(system,set,time)");
    ARIADNE_LOG(3,"initial_set="<<initial_set);

    HybridGrid grid = _grid();

    HDS initial(grid),found(grid),evolve(grid),reach(grid);

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
    initial.adjoin_outer_approximation(initial_set,_accuracy());
    ARIADNE_LOG(4,"initial.size()="<<initial.size());

    const EvolverPtrType& evolver = _get_tuned_evolver(*_system,EVOLVER_TAB_OFFSET,UPPER_SEMANTICS);
    std::list<EnclosureType> initial_enclosures = to_enclosures(initial);
    ARIADNE_LOG(4,"Starting from " << initial_enclosures.size() << " enclosures.");
    for (std::list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); ++encl_it) {
        Orbit<EnclosureType> orbit = evolver->orbit(*encl_it,hybrid_lock_to_grid_time,UPPER_SEMANTICS);
        reach.adjoin(outer_approximation(orbit.reach(),grid,_accuracy()));
        evolve.adjoin(outer_approximation(orbit.final(),grid,_accuracy()));
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

    reach.recombine();
    evolve.recombine();

    ARIADNE_LOG(4,"Reach size = " << reach.size());
    ARIADNE_LOG(4,"Final size = " << evolve.size());

    return std::make_pair(reach,evolve);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
outer_chain_reach(
		const HybridImageSet& initial_set,
		ContinuousEvolutionDirection direction) const
{
	SetApproximationType initial(_grid());
    initial.adjoin_outer_approximation(initial_set,_accuracy());

	return outer_chain_reach(initial,direction);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
outer_chain_reach(
		const SetApproximationType& initial,
		ContinuousEvolutionDirection direction) const
{
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

	ARIADNE_ASSERT_MSG(!reach.empty(),"The outer chain reachability of " << _system->name() << " is empty: check the initial set.");

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

        ARIADNE_LOG(3,"Removing the previously reached and final sets...");

        new_final.remove(final);
		new_reach.remove(reach);
	    ARIADNE_LOG(4,"Reach size after removal = "<<new_reach.size());
	    ARIADNE_LOG(4,"Final size after removal = "<<new_final.size());

	    ARIADNE_LOG(3,"Restricting the reach and final sets...");

	    _restriction->apply_to(new_final);
	    _restriction->apply_to(new_reach);
		ARIADNE_LOG(4,"Reach size after restricting = "<<new_reach.size());
		ARIADNE_LOG(4,"Final size after restricting = "<<new_final.size());

		ARIADNE_LOG(3,"Determining the new initial set from the jump sets...");

		new_reach.mince(_accuracy());
		new_final.mince(_accuracy());
		ARIADNE_LOG(4,"Reach size after mincing = "<<new_reach.size());
		ARIADNE_LOG(4,"Final size after mincing = "<<new_final.size());

       	current_initial = (direction == DIRECTION_FORWARD ?
       			_restriction->forward_jump_set(new_reach,sys) : _restriction->backward_jump_set(new_reach,sys));

        current_initial.adjoin(new_final);
        current_initial.recombine();

        reach.adjoin(new_reach);
		final.adjoin(new_final);

		reach.recombine();
		final.recombine();
    }

	ARIADNE_LOG(2,"Found a total of " << reach.size() << " reached cells.");

    return reach;
}


std::pair<HybridDenotableSet,HybridFloatVector>
HybridReachabilityAnalyser::
_lower_chain_reach_and_epsilon(
		const SystemType& system,
		const HybridImageSet& initial_set) const
{
    const unsigned EVOLVER_TAB_OFFSET = 3;

	typedef std::list<EnclosureType> EL;
	typedef std::map<DiscreteLocation,uint> HUM;

	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;
	ARIADNE_ASSERT_MSG(concurrency>0 && concurrency <= boost::thread::hardware_concurrency(),"Error: concurrency must be positive and less than the maximum allowed.");

	HybridGrid grid = _grid();
	TimeType lock_time(_settings->lock_to_grid_time,_settings->lock_to_grid_steps);

    HDS reach(grid);
	HybridSpace state_space = system.state_space();

	HybridFloatVector epsilon;
	for (HybridSpace::const_iterator hs_it = state_space.begin(); hs_it != state_space.end(); ++hs_it)
		epsilon.insert(std::pair<DiscreteLocation,Vector<Float> >(hs_it->first,Vector<Float>(hs_it->second)));

    EL initial_enclosures = enclosures_from_split_initial_set_midpoints(initial_set,min_cell_widths(grid,_accuracy()));

    initial_enclosures = _restriction->filter(initial_enclosures);

    ARIADNE_LOG(2,"Computing recurrent evolution...");

    const EvolverPtrType& evolver = _get_tuned_evolver(system,EVOLVER_TAB_OFFSET,LOWER_SEMANTICS);

    uint i=0;
    while (!initial_enclosures.empty()) {
		ARIADNE_LOG(2,"Iteration " << i++);
		ARIADNE_LOG(3,"Initial enclosures size = " << initial_enclosures.size());

		LowerReachEpsilonWorker worker(evolver,initial_enclosures,lock_time,grid,_accuracy(),concurrency);

		ARIADNE_LOG(3,"Evolving and discretising...");

		EL final_enclosures;
		std::pair<HUM,HUM> evolve_sizes;
		HDS local_reach;
		HybridFloatVector local_epsilon;
		make_ltuple<std::pair<HUM,HUM>,EL,HDS,HybridFloatVector>(evolve_sizes,final_enclosures,
				local_reach,local_epsilon) = worker.get_result();

		epsilon = max_elementwise(epsilon,local_epsilon);

	    if (_settings->enable_lower_reach_restriction_check) {
	    	ARIADNE_ASSERT_MSG(!_restriction->restricts(local_reach),
	    			"The lower reach is not inside the reachability restriction. Check the bounding domain used.");
	    }
		_restriction->apply_to(local_reach);

		if (!_settings->constraint_set.empty())
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

	HDS definitely_infeasible_reach = reach;
	HDS possibly_feasible_reach = possibly_overlapping_subset(reach,_settings->constraint_set);
	definitely_infeasible_reach.remove(possibly_feasible_reach);

	// It is not necessary to check for eps-infeasible cells if no infeasible cells are present
	if (!definitely_infeasible_reach.empty()) {

		// Identify the subspace of those locations where infeasible cells were found
		HybridSpace proj_space;
		for (HDS::locations_const_iterator loc_it = definitely_infeasible_reach.locations_begin();
				loc_it != definitely_infeasible_reach.locations_end(); ++loc_it) {
			if (!loc_it->second.empty())
				proj_space[loc_it->first] = loc_it->second.dimension();
		}

		if (_has_definitely_infeasible_subset(reach,epsilon,proj_space)) {
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
_has_definitely_infeasible_subset(
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
	HybridDenotableSet proj_possibly_feasible = _restriction->possibly_feasible_projection(proj_constraint_set);
	HybridVectorFunction proj_constraint_functions = proj_constraint_set.functions();
	HybridBoxes proj_eps_constraint_codomain = eps_codomain(proj_possibly_feasible, proj_eps, proj_constraint_functions);
	HybridConstraintSet proj_eps_constraint(proj_constraint_functions,proj_eps_constraint_codomain);

	HybridDenotableSet proj_eps_possibly_feasible_reach = possibly_overlapping_subset(proj_reach,proj_eps_constraint);
	proj_reach.remove(proj_eps_possibly_feasible_reach);

	return !proj_reach.empty();
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


std::pair<HybridDenotableSet,HybridFloatVector>
HybridReachabilityAnalyser::
lower_chain_reach_and_epsilon(const HybridImageSet& initial_set) const
{
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

    reach.recombine();

	return std::pair<SetApproximationType,HybridFloatVector>(reach,epsilon);
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

        ARIADNE_LOG(2,"Evaluating split sets...");

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

        Vector<Float> der_widths = getDerivativeWidths(*_system,loc,bbx);
        Vector<Float> mad = hmad.find(loc)->second;

        for (unsigned i=0; i < bbx.size(); ++i) {
            if (mad[i] > 0) {
                result += der_widths[i]/mad[i];
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
        const Box& cell_bx = domain_it->second;
        const Vector<Float>& mad = hmad.find(loc)->second;

        Interval left_half_val = Interval(param_val.lower(),param_val.midpoint());
        _system->substitute(RealParameter(param.name(),left_half_val));
        Vector<Float> left_der_widths = getDerivativeWidths(*_system,loc,cell_bx);
        Interval right_half_val = Interval(param_val.midpoint(),param_val.upper());
        _system->substitute(RealParameter(param.name(),right_half_val));
        Vector<Float> right_der_widths = getDerivativeWidths(*_system,loc,cell_bx);
        _system->substitute(param);

        ARIADNE_LOG(5,"Box " << *domain_it << ": " << left_der_widths << ", " << right_der_widths);

        for (unsigned i=0; i < cell_bx.size(); ++i) {
            if (mad[i] > 0) {
                left_result += left_der_widths[i]/mad[i];
                right_result += right_der_widths[i]/mad[i];
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
    : lock_to_grid_time(getLockToGridTime(sys,domain)),
      lock_to_grid_steps(1),
      constraint_set(),
      splitting_parameters_target_ratio(0.05),
      enable_lower_reach_restriction_check(false),
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


list<EnclosureType> enclosures_of_domains_midpoints(
		const list<HybridBox>& domains,
		const HybridImageSet& img_set)
{
	list<EnclosureType> result;

	for (list<HybridBox>::const_iterator box_it = domains.begin(); box_it != domains.end(); ++box_it) {
		const DiscreteLocation& loc = box_it->first;
		const Box& domain = box_it->second;
		const VectorFunction& func = img_set.find(loc)->second.function();
		Vector<Interval> domain_centre_box(domain.centre());
		result.push_back(EnclosureType(loc,ContinuousEnclosureType(func,domain_centre_box)));
	}

	return result;
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

		// Gets the first order derivatives in respect to the dynamic of the location, applied to the corresponding domain
		Vector<Interval> der_bbx = sys.dynamic_function(loc)(loc_domain);

		for (uint i=0;i<dim;i++) {
			Float maxAbsDer = abs(der_bbx[i]).upper();
			if (maxAbsDer > 0)
				result = max(result,loc_domain[i].width()/maxAbsDer);
		}
	}

	return result;
}

HybridFloatVector
getHybridMidpointAbsoluteDerivatives(
		const HybridReachabilityAnalyser::SystemType& sys,
		const HybridBoxes& bounding_domain)
{
	HybridFloatVector result;

	const HybridSpace hspace = sys.state_space();
	for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); hs_it++) {

		const DiscreteLocation& loc = hs_it->first;
		const uint dim = hs_it->second;

		// Gets the first order derivatives in respect to the dynamic of the mode, applied to the domain of the corresponding location
		Vector<Interval> der_bbx = sys.dynamic_function(loc)(bounding_domain.find(loc)->second);

		Vector<Float> local_result(dim);
		for (uint i=0;i<dim;i++)
			local_result[i] = abs(der_bbx[i]).midpoint();

		result.insert(make_pair(loc,local_result));
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


list<EnclosureType>
enclosures_from_split_initial_set_midpoints(
		const HybridImageSet initial_set,
		const HybridFloatVector max_cell_widths)
{
	list<EnclosureType> result;

	for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
		 loc_iter!=initial_set.locations_end(); ++loc_iter) {
		const DiscreteLocation& loc = loc_iter->first;
		const ImageSet& img_set = loc_iter->second;
		const VectorFunction& img_func = img_set.function();
		Vector<Float> local_max_cell_widths = max_cell_widths.find(loc)->second;

		list<Box> temporaries;
		temporaries.push_back(img_set.domain());

		/* While there are boxes:
		 * a) pick one, pop it and check it against the max_cell_widths
		 * 	 i) if larger, split it and put the couple into the temporaries
		 *   ii) if not, push the image resulting from the centre of such domain
		 */
		while (!temporaries.empty()) {
			Box domain = temporaries.back();
			temporaries.pop_back();
			Box codomain = img_func.evaluate(domain);

			bool hasBeenSplit = false;
			for (uint dim=0; dim < codomain.size(); ++dim) {
				if (codomain[dim].width() > local_max_cell_widths[dim]) {
					std::pair<Box,Box> split_boxes = domain.split(dim);
					temporaries.push_back(split_boxes.first);
					temporaries.push_back(split_boxes.second);
					hasBeenSplit = true;
					break;
				}
			}
			if (!hasBeenSplit)
				result.push_back(EnclosureType(loc,ContinuousEnclosureType(img_func,domain.centre())));
		}
	}

	return result;
}

std::list<EnclosureType>
to_enclosures(const HybridDenotableSet& reach)
{
	std::list<EnclosureType> result;

	for (HybridDenotableSet::const_iterator cell_it = reach.begin(); cell_it != reach.end(); ++cell_it)
		result.push_back(EnclosureType(cell_it->first,cell_it->second));

	return result;
}

std::list<EnclosureType>
restrict_enclosures(
		const std::list<EnclosureType> enclosures,
		const HybridDenotableSet& restriction)
{
	std::list<EnclosureType> result;

	for (std::list<EnclosureType>::const_iterator encl_it = enclosures.begin(); encl_it != enclosures.end(); ++encl_it) {
		HybridBox hbox(encl_it->first,encl_it->second.bounding_box());
		if (possibly(restriction.overlaps(hbox)))
			result.push_back(*encl_it);
	}

	return result;
}


} // namespace Ariadne
