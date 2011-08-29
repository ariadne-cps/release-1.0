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
#include "grid_set.h"

#include "orbit.h"

#include "hybrid_time.h"
#include "hybrid_automaton.h"

#include "evolver_interface.h"
#include "taylor_calculus.h"

#include "reachability_analyser.h"
#include "logging.h"

#include "graphics.h"

#include "workers.h"


namespace Ariadne {

HybridReachabilityAnalyser::
~HybridReachabilityAnalyser()
{
}

HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const HybridAutomatonInterface& system)
	: _settings(new SettingsType(system))
	, _system(system.clone())
	, free_cores(0)
{
    this->charcode = "a";
}


HybridReachabilityAnalyser::EvolverPtrType
HybridReachabilityAnalyser::
_get_tuned_evolver(
        const HybridAutomatonInterface& sys,
        int accuracy,
        unsigned ADD_TAB_OFFSET,
        Semantics semantics) const
{
    EvolverPtrType evolver(new ImageSetHybridEvolver(sys));

    evolver->verbosity = this->verbosity - ADD_TAB_OFFSET;
    evolver->tab_offset = this->tab_offset + ADD_TAB_OFFSET;

    HybridFloatVector hmad = getHybridMaximumAbsoluteDerivatives(sys,
            _settings->reachability_restriction,_settings->domain_bounds);

    evolver->tune_settings(_settings->grid,hmad,accuracy,this->free_cores,semantics);

    return evolver;
}


void
HybridReachabilityAnalyser::
tune_settings(
        const HybridBoxes& domain,
        const Set<Identifier>& locked_params_ids,
        const HybridGridTreeSetPtr& outer_approximation,
        const HybridGridTreeSetPtr& reachability_restriction,
        const HybridConstraintSet& constraint_set,
        bool EQUAL_GRID_FOR_ALL_LOCATIONS,
        int accuracy,
        unsigned free_cores,
        Semantics semantics)
{
    ARIADNE_LOG(1, "Tuning settings for analysis...");

    this->free_cores = free_cores;

    _settings->maximum_grid_depth = accuracy;
    _settings->reachability_restriction = reachability_restriction;

    _settings->domain_bounds = domain;
    ARIADNE_LOG(2, "Domain: " << domain);
    _settings->constraint_set = constraint_set;
    ARIADNE_LOG(2, "Constraint set: " << constraint_set);
    _settings->lock_to_grid_time = getLockToGridTime(*_system,domain);
    ARIADNE_LOG(2, "Lock to grid time: " << _settings->lock_to_grid_time);
    _settings->locked_parameters_ids = locked_params_ids;
    ARIADNE_LOG(2, "Locked parameters IDs: " << locked_params_ids);

    ARIADNE_LOG(2, "Derivatives evaluation source: " << (outer_approximation ? "Outer approximation" : "Domain box"));
    HybridFloatVector hmad = getHybridMaximumAbsoluteDerivatives(*_system,outer_approximation,domain);
    ARIADNE_LOG(2, "System parameters: " << dynamic_cast<HybridAutomaton&>(*_system).parameters());
    ARIADNE_LOG(2, "Derivatives bounds: " << hmad);
    _settings->grid = HybridGrid(getHybridGrid(hmad,domain,EQUAL_GRID_FOR_ALL_LOCATIONS));
    ARIADNE_LOG(2, "Grid lengths: " << _settings->grid.lengths());

    // Checking for grid and reachability restriction grid equality is required for upper semantics only
    if (reachability_restriction && semantics == UPPER_SEMANTICS)
        ARIADNE_ASSERT_MSG(_settings->grid == reachability_restriction->grid(),
            "For upper semantics, the analyser grid and the reachability restriction grid must match.");
}


void
HybridReachabilityAnalyser::
_plot_reach(
		const HybridGridTreeSet& reach,
		string plot_dirpath,
		string name_prefix) const
{
	char mgd_char[10];
	sprintf(mgd_char,"%i",_settings->maximum_grid_depth);
	name_prefix.append(mgd_char);
	plot(plot_dirpath,name_prefix,reach);
}


HybridGridTreeSet
HybridReachabilityAnalyser::
initial_cells_set(const HybridImageSet& initial_enclosure_set) const
{
    const int& accuracy = _settings->maximum_grid_depth;

    HybridGridTreeSet result(_settings->grid);

    result.adjoin_outer_approximation(initial_enclosure_set,accuracy);
    result.mince(accuracy);
    if (_settings->reachability_restriction) {
        result.restrict(*_settings->reachability_restriction);
    }

    return result;
}


HybridGridTreeSet
HybridReachabilityAnalyser::
initial_cells_set(const HybridConstraintSet& initial_constraint_set) const
{
    ARIADNE_ASSERT_MSG(_settings->reachability_restriction,
            "A reachability restriction is required for computing the initial cell set from a constraint set.");

    const int& accuracy = _settings->maximum_grid_depth;

    HybridGridTreeSet result(_settings->grid);

    result = *_settings->reachability_restriction;
    result.mince(accuracy);
    result.remove(definitely_covered_cells(result,initial_constraint_set));

    return result;
}


std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridReachabilityAnalyser::
_upper_reach_evolve(
		const SystemType& sys,
        const HybridGridTreeSet& initial_set,
        const HybridTime& time,
        const int accuracy) const
{
	std::list<EnclosureType> initial_enclosures = cells_to_smallest_enclosures(initial_set,accuracy);

	return _upper_reach_evolve(sys,initial_enclosures,time,false,DIRECTION_FORWARD,accuracy);
}

std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridReachabilityAnalyser::
_upper_reach_evolve(
		const SystemType& sys,
        const list<EnclosureType>& initial_enclosures,
        const HybridTime& time,
        bool enable_premature_termination_on_blocking_event,
        ContinuousEvolutionDirection direction,
        int accuracy) const
{
    const unsigned EVOLVER_TAB_OFFSET = 4;

	std::pair<GTS,GTS> result;
	GTS& reach = result.first;
	GTS& evolve = result.second;

	ARIADNE_LOG(4,"Evolving and discretising...");

	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;
	ARIADNE_ASSERT_MSG(concurrency>0 && concurrency <= boost::thread::hardware_concurrency(),"Error: concurrency must be positive and less than the maximum allowed.");

	const HybridGrid& grid = _settings->grid;

	const EvolverPtrType& evolver = _get_tuned_evolver(sys,accuracy,EVOLVER_TAB_OFFSET,UPPER_SEMANTICS);
	UpperReachEvolveWorker worker(evolver,initial_enclosures,time,
			grid,accuracy,enable_premature_termination_on_blocking_event,direction,concurrency);
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

    HybridGrid grid = _settings->grid;
    const int accuracy = _settings->maximum_grid_depth;

    GTS reach(grid); GTS evolve(grid);
    // No reachability restriction involved
    GTS reachability_restriction(grid);

    list<EnclosureType> initial_enclosures = enclosures_from_split_domain_midpoints(initial_set,min_cell_widths(grid,accuracy));

    const EvolverPtrType& evolver = _get_tuned_evolver(*_system,accuracy,EVOLVER_TAB_OFFSET,LOWER_SEMANTICS);
	ARIADNE_LOG(3,"Computing evolution...");
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        Orbit<EnclosureType> orbit = evolver->orbit(*encl_it,time,LOWER_SEMANTICS);
        reach.adjoin(outer_approximation(orbit.reach(),grid,accuracy));
        evolve.adjoin(outer_approximation(orbit.final(),grid,accuracy));
    }

    if (_settings->reachability_restriction) {
        reach.restrict(*_settings->reachability_restriction);
        evolve.restrict(*_settings->reachability_restriction);
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
 
	GTS reach, evolve;
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

	GTS reach, evolve;
	make_lpair(reach,evolve) = upper_reach_evolve(initial_set, time);

    return reach;
}


std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
upper_reach_evolve(
        const HybridImageSet& initial_set,
        const TimeType& time) const
{
    const unsigned EVOLVER_TAB_OFFSET = 3;

    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach_evolve(system,set,time)");
    ARIADNE_LOG(3,"initial_set="<<initial_set);

    HybridGrid grid = _settings->grid;
    const int accuracy = _settings->maximum_grid_depth;

    GTS initial(grid),found(grid),evolve(grid),reach(grid);

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
    initial.adjoin_outer_approximation(initial_set,accuracy);

    const EvolverPtrType& evolver = _get_tuned_evolver(*_system,accuracy,EVOLVER_TAB_OFFSET,UPPER_SEMANTICS);
    std::list<EnclosureType> initial_enclosures = cells_to_smallest_enclosures(initial,accuracy);
    for (std::list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); ++encl_it) {
        Orbit<EnclosureType> orbit = evolver->orbit(*encl_it,hybrid_lock_to_grid_time,UPPER_SEMANTICS);
        reach.adjoin(outer_approximation(orbit.reach(),grid,accuracy));
        evolve.adjoin(outer_approximation(orbit.final(),grid,accuracy));
    }

    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...");
        make_lpair(found,evolve) = _upper_reach_evolve(*_system,evolve,hybrid_lock_to_grid_time,accuracy);
        ARIADNE_LOG(4,"found.size()="<<found.size());
        ARIADNE_LOG(4,"evolve.size()="<<evolve.size());
        reach.adjoin(found);
        ARIADNE_LOG(3,"found "<<found.size()<<" cells.");
    }

    ARIADNE_LOG(5,"remaining_time="<<remaining_time);
    if(!evolve.empty() && remaining_time > 0) {
        ARIADNE_LOG(3,"computing evolution for the remaining time...");
        make_lpair(found,evolve) = _upper_reach_evolve(*_system,evolve,hybrid_remaining_time,accuracy);
        reach.adjoin(found);
    }

    if (_settings->reachability_restriction) {
        reach.restrict(*_settings->reachability_restriction);
        evolve.restrict(*_settings->reachability_restriction);
    }

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
    ARIADNE_LOG(1,"Performing outer chain reachability from an enclosure set...");

	HybridGridTreeSet initial(_settings->grid);
    initial.adjoin_outer_approximation(initial_set,_settings->maximum_grid_depth);
    std::list<EnclosureType> initial_enclosures = cells_to_smallest_enclosures(initial,_settings->maximum_grid_depth);

	return _outer_chain_reach(initial_enclosures,direction);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
outer_chain_reach(
		const HybridGridTreeSet& initial,
		ContinuousEvolutionDirection direction) const
{
    ARIADNE_LOG(1,"Performing outer chain reachability from a grid set...");

    std::list<EnclosureType> initial_enclosures = cells_to_smallest_enclosures(initial,_settings->maximum_grid_depth);

	return _outer_chain_reach(initial_enclosures,direction);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
_outer_chain_reach(const std::list<EnclosureType>& initial_enclosures,
		ContinuousEvolutionDirection direction) const
{
	HybridGridTreeSet reach;

	RealParameterSet original_parameters = _system->nonsingleton_parameters();

	ARIADNE_LOG(2,"Splitting the nonsingleton nonlocked parameter set...");

	std::list<RealParameterSet> split_parameter_set_list = _getSplitParameterSetList();

	try {
		uint i = 0;
		for (std::list<RealParameterSet>::const_iterator set_it = split_parameter_set_list.begin(); set_it != split_parameter_set_list.end(); ++set_it) {
			ARIADNE_LOG(1,"Split parameter set #" << ++i << "/" << split_parameter_set_list.size() << " : " << *set_it);

			_system->substitute_all(*set_it);

			HybridGridTreeSet local_reach = _outer_chain_reach_splitted(*_system,initial_enclosures,direction);

			reach.adjoin(local_reach);
		}
	} catch (ReachOutOfDomainException ex) {
		_system->substitute_all(original_parameters);
		throw ex;
	} catch (ReachUnsatisfiesConstraintException ex) {
	    _system->substitute_all(original_parameters);
		throw ex;
	}

	_system->substitute_all(original_parameters);

	ARIADNE_ASSERT_MSG(!reach.empty(),"The outer chain reachability of " << _system->name() << " is empty: check the initial set.");

	return reach;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
_outer_chain_reach_splitted(
		const SystemType& sys,
		const std::list<EnclosureType>& initial_enclosures,
		ContinuousEvolutionDirection direction) const
{
    bool no_restriction = !_settings->reachability_restriction;

	ARIADNE_ASSERT_MSG(!(direction == DIRECTION_BACKWARD && no_restriction),
			"The reachability restriction must be set if backward reachability is to be performed.");

    const Float& lock_to_grid_time = _settings->lock_to_grid_time;
    const int& lock_to_grid_steps = _settings->lock_to_grid_steps;
    const int& maximum_grid_depth = _settings->maximum_grid_depth;

    HybridGrid grid = _settings->grid;
    HybridGridTreeSet new_final(grid), new_reach(grid), reach(grid), final(grid);

    list<EnclosureType> working_enclosures = initial_enclosures;

    ARIADNE_LOG(2,"Computing recurrent " << (direction == DIRECTION_FORWARD ? "forward" : "backward") << " evolution...");
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);

    uint i=0;
    while (!working_enclosures.empty())
	{
    	ARIADNE_LOG(2,"Iteration " << i++);

        ARIADNE_LOG(3,"Initial enclosures size = " << working_enclosures.size());

        static const bool ignore_activations = true;
        make_lpair(new_reach,new_final) = _upper_reach_evolve(sys,working_enclosures,
        		hybrid_lock_to_grid_time,ignore_activations,direction,maximum_grid_depth);

        new_final.remove(final);
		new_reach.remove(reach);
	    ARIADNE_LOG(3,"Reach size after removal = "<<new_reach.size());
	    ARIADNE_LOG(3,"Final size after removal = "<<new_final.size());

	    if (!no_restriction) {
	    	new_final.restrict(*_settings->reachability_restriction);
	    	new_reach.restrict(*_settings->reachability_restriction);
		    ARIADNE_LOG(3,"Reach size after restricting = "<<new_reach.size());
		    ARIADNE_LOG(3,"Final size after restricting = "<<new_final.size());
	    }

		new_reach.mince(maximum_grid_depth);
		new_final.mince(maximum_grid_depth);
		ARIADNE_LOG(3,"Reach size after mincing = "<<new_reach.size());
		ARIADNE_LOG(3,"Final size after mincing = "<<new_final.size());

        working_enclosures.clear();

        if (direction == DIRECTION_FORWARD)
        	_outer_chain_reach_forward_pushTargetCells(sys,new_reach,working_enclosures,no_restriction);
        else
        	_outer_chain_reach_backward_pushSourceCells(sys,new_reach,working_enclosures);

		_outer_chain_reach_pushLocalFinalCells(new_final,working_enclosures,no_restriction);

        reach.adjoin(new_reach);
		final.adjoin(new_final);

		reach.recombine();
		final.recombine();
    }

	ARIADNE_LOG(2,"Found a total of " << reach.size() << " reached cells.");

    return reach;
}


void
HybridReachabilityAnalyser::
_outer_chain_reach_forward_pushTargetCells(
		const SystemType& system,
		const HybridGridTreeSet& reachCells,
		std::list<EnclosureType>& result_enclosures,
		bool use_domain_checking) const
{
	for (HybridGridTreeSet::const_iterator cell_it = reachCells.begin(); cell_it != reachCells.end(); cell_it++)
	{
		const DiscreteLocation& loc = cell_it->first;
		const Box& bx = cell_it->second.box();
		const Box& domain = _settings->domain_bounds[loc];

		ARIADNE_LOG(4,"Checking box "<< bx <<" in location " << loc.name());

		if (use_domain_checking && bx.disjoint(domain))
			throw ReachOutOfDomainException("a reach enclosure is outside the domain");

		if (!_outer_chain_reach_isOutsideInvariants(system,loc,bx)) {
			_outer_chain_reach_forward_pushTargetEnclosures(system,loc,ContinuousEnclosureType(bx),
					_settings->grid,result_enclosures,use_domain_checking);
		}
	}
}


void
HybridReachabilityAnalyser::
_outer_chain_reach_backward_pushSourceCells(
		const SystemType& system,
		const HybridGridTreeSet& targetCells,
		std::list<EnclosureType>& result_enclosures) const
{
    HybridGridTreeSet sourceCellsOverapprox = *_settings->reachability_restriction;

	// Adopt the minimum granularity when checking the cells
	sourceCellsOverapprox.mince(_settings->maximum_grid_depth);

	for (HybridGridTreeSet::const_iterator cell_it = sourceCellsOverapprox.begin(); cell_it != sourceCellsOverapprox.end(); ++cell_it) {
		const DiscreteLocation& loc = cell_it->first;
		const Box& bx = cell_it->second.box();

		ARIADNE_LOG(4,"Checking box "<< bx <<" in location " << loc.name());

		if (!_outer_chain_reach_isOutsideInvariants(system,loc,bx))
			_outer_chain_reach_backward_pushSourceEnclosures(system,loc,ContinuousEnclosureType(bx),
			        targetCells,_settings->grid,result_enclosures);
	}
}

bool
HybridReachabilityAnalyser::
_outer_chain_reach_isOutsideInvariants(
		const SystemType& system,
		const DiscreteLocation& location,
		const Box& bx) const
{
	ARIADNE_LOG(5,"Checking invariants...");

	boost::shared_ptr<CalculusInterface<TaylorModel> > tc(new TaylorCalculus());

	Set<DiscreteEvent> events = system.events(location);
	for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {
		const DiscreteEvent event = *event_it;
		EventKind kind = system.event_kind(location,event);
		if (kind == INVARIANT) {
			const ScalarFunction& activation = system.invariant_function(location,event);
			tribool is_active = tc->active(VectorFunction(1,activation),bx);
			if (definitely(is_active)) {
				ARIADNE_LOG(5,"Invariant '" << event.name() << "' is definitely active: transitions will not be checked.");
				return true;
			}
		}
	}

	return false;
}

void
HybridReachabilityAnalyser::
_outer_chain_reach_forward_pushTargetEnclosures(
		const SystemType& system,
		const DiscreteLocation& sourceLocation,
		const ContinuousEnclosureType& sourceEnclosure,
		const HybridGrid& grid,
		std::list<EnclosureType>& result_enclosures,
		bool use_domain_checking) const
{
    long numCellDivisions = (1<<_settings->maximum_grid_depth);

	ARIADNE_LOG(5,"Checking transitions...");

	boost::shared_ptr<CalculusInterface<TaylorModel> > tc(new TaylorCalculus());

	Set<DiscreteEvent> events = system.events(sourceLocation);
	for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {
		const DiscreteEvent event = *event_it;
		EventKind kind = system.event_kind(sourceLocation,event);
		if (kind == URGENT || kind == PERMISSIVE) {
			ARIADNE_LOG(6,"Target: "<<system.target(sourceLocation,event)<<", Kind: " << kind);

			if (_is_transition_feasible(system.guard_function(sourceLocation,event),kind,
					system.dynamic_function(sourceLocation),sourceEnclosure,UPPER_SEMANTICS)) {
				const DiscreteLocation& target_loc = system.target(sourceLocation,event);
				const ContinuousEnclosureType target_encl = tc->reset_step(
						system.reset_function(sourceLocation,event),sourceEnclosure);
				const Box& target_bounding = _settings->domain_bounds[target_loc];
				const Vector<Float> minTargetCellWidths = grid[target_loc].lengths()/numCellDivisions;

				pushSplitTargetEnclosures(result_enclosures,target_loc,target_encl,minTargetCellWidths,target_bounding,use_domain_checking);
			}
		}
	}
}

void
HybridReachabilityAnalyser::
_outer_chain_reach_backward_pushSourceEnclosures(
		const SystemType& system,
		const DiscreteLocation& sourceLocation,
		const ContinuousEnclosureType& sourceEnclosure,
		const HybridGridTreeSet& targetCells,
		const HybridGrid& grid,
		std::list<EnclosureType>& result_enclosures) const
{
	ARIADNE_LOG(5,"Checking transitions...");

	boost::shared_ptr<CalculusInterface<TaylorModel> > tc(new TaylorCalculus());

	Set<DiscreteEvent> events = system.events(sourceLocation);
	// In order to add a source hybrid box, just one possibly overlapping target enclosure suffices
	for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {
		const DiscreteEvent event = *event_it;
		EventKind kind = system.event_kind(sourceLocation,event);
		if (kind == URGENT || kind == PERMISSIVE) {
			ARIADNE_LOG(6,"Target: "<<system.target(sourceLocation,event)<<", Kind: " << kind);

			if (_is_transition_feasible(system.guard_function(sourceLocation,event),kind,
					system.dynamic_function(sourceLocation),sourceEnclosure,UPPER_SEMANTICS)) {
				const DiscreteLocation& target_loc = system.target(sourceLocation,event);
				const ContinuousEnclosureType target_encl = tc->reset_step(
						system.reset_function(sourceLocation,event),sourceEnclosure);
				const HybridBox targetHBox(target_loc,target_encl.bounding_box());

				if (possibly(targetCells.overlaps(targetHBox))) {
					ARIADNE_LOG(6,"Target enclosure overlaps with target cells: adding the source enclosure.");
					result_enclosures.push_back(EnclosureType(sourceLocation,sourceEnclosure));
					break;
				}
			}
		}
	}
}

bool
HybridReachabilityAnalyser::
_is_transition_feasible(
		const ScalarFunction& activation,
		EventKind event_kind,
		const VectorFunction& dynamic,
		const ContinuousEnclosureType& source,
		Semantics semantics) const
{
	bool result = false;

	const bool is_urgent = (event_kind == URGENT);

	boost::shared_ptr<CalculusInterface<TaylorModel> > tc(new TaylorCalculus());

	tribool is_guard_active = tc->active(VectorFunction(1,activation),source);

	ARIADNE_LOG(6,"Guard activity: " << is_guard_active);

	/*
	 * a) If the guard is definitely active and the transition is urgent, then we are definitely outside the related invariant
	 * b) If the transition is not urgent, it suffices to have a possibly active guard: we then must perform the transition
	 * c) If the transition is urgent and the guard is only possibly active, we check the crossing:
	 *    i) If it is negative, then no transition is possible
	 *	 ii) If it is possibly positive, then we must take the transition
	 */

	if (definitely(is_guard_active) && is_urgent) {
		ARIADNE_LOG(6,"Definitely active and urgent: the set is outside the implicit invariant of the transition, infeasible.");
		result = false;
	} else if (possibly(is_guard_active) && !is_urgent) {
		ARIADNE_LOG(6,"Possibly active and permissive: feasible.");
		result = true;
	} else if (possibly(is_guard_active) && is_urgent) {
		ARIADNE_LOG(6,"Possibly active and urgent: checking whether the crossing is nonnegative...");
		tribool positive_crossing = positively_crossing(source.bounding_box(),dynamic,activation);
		if (possibly(positive_crossing)) {
			ARIADNE_LOG(6,"Possibly positive: feasible.");
			result = true;
		} else {
			ARIADNE_LOG(6,"Negative: infeasible.");
			result = false;
		}
	} else {
		ARIADNE_LOG(6,"Inactive: infeasible.");
		result = false;
	}

	return result;
}

void
HybridReachabilityAnalyser::
_outer_chain_reach_pushLocalFinalCells(
		const HybridGridTreeSet& finalCells,
		std::list<EnclosureType>& result_enclosures,
		bool use_domain_checking) const
{
	for (GTS::const_iterator cell_it = finalCells.begin(); cell_it != finalCells.end(); ++cell_it) {
		const DiscreteLocation& loc = cell_it->first;
		const Box& domain = _settings->domain_bounds[loc];
		const Box& bx = cell_it->second.box();

		if (use_domain_checking && bx.disjoint(domain)) {
			ARIADNE_LOG(4,"Discarding enclosure " << bx << " from final cell outside the domain in location " << loc.name() <<".");
			throw ReachOutOfDomainException("a final cell is outside the domain");
		} else {
			result_enclosures.push_back(EnclosureType(cell_it->first,ContinuousEnclosureType(cell_it->second.box())));
		}
	}
}


std::pair<HybridGridTreeSet,HybridFloatVector>
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

	HybridGrid grid = _settings->grid;
	TimeType lock_time(_settings->lock_to_grid_time,_settings->lock_to_grid_steps);
	const int accuracy = _settings->maximum_grid_depth;

	bool have_restriction = _settings->reachability_restriction;

    HybridGridTreeSet reach(grid);
	HybridSpace state_space = system.state_space();

	HybridFloatVector epsilon;
	for (HybridSpace::const_iterator hs_it = state_space.begin(); hs_it != state_space.end(); ++hs_it)
		epsilon.insert(std::pair<DiscreteLocation,Vector<Float> >(hs_it->first,Vector<Float>(hs_it->second)));

    EL initial_enclosures = enclosures_from_split_domain_midpoints(initial_set,min_cell_widths(grid,accuracy));

    if (have_restriction)
    	initial_enclosures = restrict_enclosures(initial_enclosures,*_settings->reachability_restriction);

    ARIADNE_LOG(2,"Computing recurrent evolution...");

    const EvolverPtrType& evolver = _get_tuned_evolver(system,accuracy,EVOLVER_TAB_OFFSET,LOWER_SEMANTICS);

    uint i=0;
    while (!initial_enclosures.empty()) {
		ARIADNE_LOG(2,"Iteration " << i++);
		EL final_enclosures;
		std::pair<HUM,HUM> evolve_sizes;

		HUM& adjoined_evolve_sizes = evolve_sizes.first;
		HUM& superposed_evolve_sizes = evolve_sizes.second;

		GTS local_reach;
		HybridFloatVector local_epsilon;


		ARIADNE_LOG(3,"Initial enclosures size = " << initial_enclosures.size());

		LowerReachEpsilonWorker worker(evolver,initial_enclosures,lock_time,grid,accuracy,concurrency);

		ARIADNE_LOG(3,"Evolving and discretising...");

		make_ltuple<std::pair<HUM,HUM>,EL,GTS,HybridFloatVector>(evolve_sizes,final_enclosures,
				local_reach,local_epsilon) = worker.get_result();

		epsilon = max_elementwise(epsilon,local_epsilon);

		if (!_settings->constraint_set.empty()) {
			HybridGridTreeSet local_reachability_restriction;
			if (!have_restriction)
				local_reachability_restriction.adjoin_outer_approximation(_settings->domain_bounds,accuracy);
			else
			    local_reachability_restriction = *_settings->reachability_restriction;

			HybridGridTreeSet possibly_feasible_cells = Ariadne::possibly_feasible_cells(local_reach,
			        _settings->constraint_set,local_epsilon,local_reachability_restriction,accuracy);

			if (possibly_feasible_cells.size() < local_reach.size()) {
				throw ReachUnsatisfiesConstraintException("The lower reached region partially does not satisfy the constraint.");
			}
		}

		ARIADNE_LOG(3,"Reach size before removal = " << local_reach.size());

		local_reach.remove(reach);
		ARIADNE_LOG(3,"Reach size after removal  = " << local_reach.size());
		if (local_reach.empty())
			break;

		reach.adjoin(local_reach);

		ARIADNE_LOG(3,"Final enclosures size = " << final_enclosures.size());

		bool use_domain_checking = !have_restriction;
		_filter_enclosures(final_enclosures,initial_enclosures,
				adjoined_evolve_sizes,superposed_evolve_sizes,use_domain_checking);
	}

	return std::pair<HybridGridTreeSet,HybridFloatVector>(reach,epsilon);
}


void
HybridReachabilityAnalyser::
_filter_enclosures(
		std::list<EnclosureType>& final_enclosures,
		std::list<EnclosureType>& initial_enclosures,
		const std::map<DiscreteLocation,uint>& adjoined_evolve_sizes,
		const std::map<DiscreteLocation,uint>& superposed_evolve_sizes,
		bool use_domain_checking) const
{
	while (!final_enclosures.empty()) {
		EnclosureType encl = final_enclosures.front();
		final_enclosures.pop_front();

		const DiscreteLocation& loc = encl.location();
		const Box& encl_box = encl.continuous_state_set().bounding_box();

		if (use_domain_checking && encl_box.disjoint(_settings->domain_bounds[loc])) {
			ARIADNE_FAIL_MSG("Found an enclosure in location " << loc.name() << " with bounding box " << encl.continuous_state_set().bounding_box() <<
					" lying outside the domain in lower semantics: the domain is incorrect.");
		}

		/* If pruning is to be performed, push only a fraction of the final_enclosures into the initial_enclosures;
		 * otherwise, push indiscriminately.
		 */
		if (_settings->enable_lower_pruning) {
			Float coverage_ratio = (Float)adjoined_evolve_sizes.find(loc)->second/(Float)superposed_evolve_sizes.find(loc)->second;

			if (initial_enclosures.size() <= 2 || rand() < 2*coverage_ratio*RAND_MAX)
				initial_enclosures.push_back(encl);
		} else
			initial_enclosures.push_back(encl);
	}
}


std::pair<HybridGridTreeSet,HybridFloatVector>
HybridReachabilityAnalyser::
lower_chain_reach_and_epsilon(const HybridImageSet& initial_set) const
{
	HybridGrid grid = _settings->grid;
	HybridGridTreeSet reach(grid);
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

			HybridGridTreeSet local_reach(grid);
			HybridFloatVector local_epsilon;
			make_lpair<HybridGridTreeSet,HybridFloatVector>(local_reach,local_epsilon) =
					_lower_chain_reach_and_epsilon(*_system,initial_set);

			reach.adjoin(local_reach);
			epsilon = max_elementwise(epsilon,local_epsilon);

			ARIADNE_LOG(2,"Epsilon: " << local_epsilon);
		}

	} catch (ReachUnsatisfiesConstraintException ex) {
		_system->substitute_all(original_parameters);
		throw ex;
	}

    _system->substitute_all(original_parameters);

	return std::pair<HybridGridTreeSet,HybridFloatVector>(reach,epsilon);
}


std::list<RealParameterSet>
HybridReachabilityAnalyser::
_getSplitParameterSetList() const
{
    std::list<RealParameterSet> result_parameter_set_list;

    std::list<std::pair<Float,RealParameterSet> > working_scored_parameter_set_list;

    const int& accuracy = _settings->maximum_grid_depth;

    RealParameterSet initial_parameter_set = _system->nonsingleton_parameters();
    remove_nonlocked_parameters(initial_parameter_set,_settings->locked_parameters_ids);

    HybridFloatVector hmad = getHybridMaximumAbsoluteDerivatives(*_system,
            _settings->reachability_restriction,_settings->domain_bounds);

    ARIADNE_LOG(2,"Hybrid Maximum Absolute Derivatives: " << hmad);

    if (!initial_parameter_set.empty()) {

        HybridGridTreeSet discretised_domain;
        if (_settings->reachability_restriction) {
            ARIADNE_LOG(2,"Creating the discretised domain from the reachability restriction...");
            discretised_domain.adjoin(*_settings->reachability_restriction);
        } else {
            ARIADNE_LOG(2,"Creating the discretised domain from the bounding domain...");
            discretised_domain.adjoin_outer_approximation(_settings->domain_bounds,accuracy);
        }
        /* We recombine since large cells lie in the center of the region (at least when the tree set represents a
         * reachability restriction), and in that region we usually are not interested in having very fine evolution.
         * This solution accomodates a more reasonable computation time for the split calculation. */
        discretised_domain.recombine();

        Float initial_score = _getDerivativeWidthsScore(hmad,discretised_domain);

        ARIADNE_LOG(2,"Discretised domain size: " << discretised_domain.size());
        ARIADNE_LOG(2,"Evaluating split sets...");

        working_scored_parameter_set_list.push_back(std::pair<Float,RealParameterSet>(initial_score,initial_parameter_set));

        do {
            _updateSplitParameterSetLists(working_scored_parameter_set_list,result_parameter_set_list,
                    initial_score,hmad,discretised_domain);

            ARIADNE_LOG(3,"Working scored parameter set list size = " << working_scored_parameter_set_list.size());
            ARIADNE_LOG(3,"Result parameter set list size = " << result_parameter_set_list.size());
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
_getDerivativeWidthsScore(
        const HybridFloatVector& hmad,
        const HybridGridTreeSet& discretised_domain) const
{
    Float result = 0;

    // Reads: \sum_{cells} \sum_{dim} der_widths(cell_loc,dim)/hmad(cell_loc,dim), under hmad(cell_loc,dim) > 0

    for (HybridGridTreeSet::const_iterator cell_it = discretised_domain.begin();
            cell_it != discretised_domain.end(); ++cell_it) {
        const DiscreteLocation& loc = cell_it->box().first;
        const Box& cell_bx = cell_it->box().second;

        Vector<Float> der_widths = getDerivativeWidths(*_system,loc,cell_bx);
        Vector<Float> mad = hmad.find(loc)->second;

        for (unsigned i=0; i < cell_bx.size(); ++i) {
            if (mad[i] > 0) {
                result += der_widths[i]/mad[i];
            }
        }
    }

    return result/discretised_domain.size();
}


std::pair<Float,Float>
HybridReachabilityAnalyser::
_getSplitDerivativeWidthsScores(
        const RealParameter& param,
        const HybridFloatVector& hmad,
        const HybridGridTreeSet& discretised_domain) const
{
    const Interval& param_val = param.value();

    Float left_result = 0;
    Float right_result = 0;

    for (HybridGridTreeSet::const_iterator cell_it = discretised_domain.begin();
            cell_it != discretised_domain.end(); ++cell_it) {
        const DiscreteLocation& loc = cell_it->box().first;
        const Box& cell_bx = cell_it->box().second;
        const Vector<Float>& mad = hmad.find(loc)->second;

        Interval left_half_val = Interval(param_val.lower(),param_val.midpoint());
        _system->substitute(RealParameter(param.name(),left_half_val));
        Vector<Float> left_der_widths = getDerivativeWidths(*_system,loc,cell_bx);
        Interval right_half_val = Interval(param_val.midpoint(),param_val.upper());
        _system->substitute(RealParameter(param.name(),right_half_val));
        Vector<Float> right_der_widths = getDerivativeWidths(*_system,loc,cell_bx);
        _system->substitute(param);

        ARIADNE_LOG(5,"Cell " << cell_it->box() << ": " << left_der_widths << ", " << right_der_widths);

        for (unsigned i=0; i < cell_bx.size(); ++i) {
            if (mad[i] > 0) {
                left_result += left_der_widths[i]/mad[i];
                right_result += right_der_widths[i]/mad[i];
            }
        }
    }

    return std::pair<Float,Float>(left_result/discretised_domain.size(),right_result/discretised_domain.size());
}


void
HybridReachabilityAnalyser::
_updateSplitParameterSetLists(
        std::list<std::pair<Float,RealParameterSet> >& working_scored_parameter_set_list,
        std::list<RealParameterSet>& result_parameter_set_list,
        const Float& initial_score,
        const HybridFloatVector& hmad,
        const HybridGridTreeSet& discretised_domain) const
{
    std::pair<Float,RealParameterSet> working_scored_parameter_set = working_scored_parameter_set_list.back();
    working_scored_parameter_set_list.pop_back();

    const Float& previous_score = working_scored_parameter_set.first;
    const RealParameterSet& working_parameter_set = working_scored_parameter_set.second;

    ARIADNE_LOG(3,"Previous score: " << previous_score);
    ARIADNE_LOG(3,"Working parameter set: " << working_parameter_set);

    RealParameterSet original_parameter_set = _system->nonsingleton_parameters();

    _system->substitute_all(working_parameter_set);

    Identifier best_param_id = "";
    Interval best_interval;
    std::pair<Float,Float> best_scores(std::numeric_limits<Float>::infinity(),std::numeric_limits<Float>::infinity());

    for (RealParameterSet::const_iterator param_it = working_parameter_set.begin();
            param_it != working_parameter_set.end(); ++param_it) {

        ARIADNE_LOG(4,"Parameter " << *param_it);

        std::pair<Float,Float> local_scores = _getSplitDerivativeWidthsScores(*param_it,hmad,discretised_domain);

        ARIADNE_LOG(4,"Scores: " << local_scores);

        if (max(local_scores.first,local_scores.second) < max(best_scores.first,best_scores.second)) {
            best_scores = local_scores;
            best_param_id = param_it->name();
            best_interval = param_it->value();
        }
    }

    _system->substitute_all(original_parameter_set);

    ARIADNE_LOG(3,"Best result: " << best_param_id << " " << best_interval << " with scores " << best_scores);

    Float ratio = (previous_score-max(best_scores.first,best_scores.second))/initial_score;
    ARIADNE_LOG(3,"Current ratio: " << ratio);

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

        ARIADNE_LOG(3,"New parameters left: " << new_parameters_left);
        ARIADNE_LOG(3,"New parameters right: " << new_parameters_right);

        working_scored_parameter_set_list.push_front(std::pair<Float,RealParameterSet>(best_scores.first,new_parameters_left));
        working_scored_parameter_set_list.push_front(std::pair<Float,RealParameterSet>(best_scores.second,new_parameters_right));
    } else {
        result_parameter_set_list.push_back(working_parameter_set);
    }

    //getchar();
}


HybridReachabilityAnalyserSettings::HybridReachabilityAnalyserSettings(const SystemType& sys)
    : lock_to_grid_time(1.0),
      lock_to_grid_steps(1),
      maximum_grid_depth(6),
      domain_bounds(unbounded_hybrid_boxes(sys.state_space())),
      constraint_set(),
      reachability_restriction(),
      grid(HybridGrid(sys.state_space())),
      splitting_parameters_target_ratio(0.05),
      enable_lower_pruning(true)
{
}


std::ostream&
operator<<(std::ostream& os, const HybridReachabilityAnalyserSettings& s)
{
    os << "DiscreteEvolutionSettings"
       << "(\n  lock_to_grid_steps=" << s.lock_to_grid_steps
       << ",\n  lock_to_grid_time=" << s.lock_to_grid_time
       << ",\n  maximum_grid_depth=" << s.maximum_grid_depth
       << ",\n  bounding_domain=" << s.domain_bounds
       << ",\n  grid=" << s.grid
       << ",\n  splitting_constants_target_ratio=" << s.splitting_parameters_target_ratio
       << ",\n  enable_lower_pruning=" << s.enable_lower_pruning
       << "\n)\n";
    return os;
}


HybridFloatVector min_cell_widths(
		const HybridGrid& grid,
		int maximum_grid_depth)
{
	HybridFloatVector result;

	for (HybridGrid::const_iterator grid_it = grid.begin(); grid_it != grid.end(); ++grid_it)
		result.insert(std::pair<DiscreteLocation,Vector<Float> >(grid_it->first,grid_it->second.lengths()/(1 << maximum_grid_depth)));

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


list<EnclosureType>
enclosures_from_split_domain_midpoints(
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


void
pushSplitTargetEnclosures(
		std::list<EnclosureType>& initial_enclosures,
		const DiscreteLocation& target_loc,
		const ContinuousEnclosureType& target_encl,
		const Vector<Float>& minTargetCellWidths,
		const Box& target_domain_constraint,
		bool use_domain_checking)
{
	const Box& target_encl_box = target_encl.bounding_box();

	if (use_domain_checking && target_encl_box.disjoint(target_domain_constraint)) {
		throw ReachOutOfDomainException("A split target enclosure is out of the domain.");
	}

	bool hasSplit = false;
	for (uint i=0; i < minTargetCellWidths.size(); i++) {
		// If the enclosure has width larger than that of the minimum cell, split on that dimension
		if (target_encl_box[i].width() > minTargetCellWidths[i]) {
			hasSplit = true;
			std::pair<ContinuousEnclosureType,ContinuousEnclosureType> split_sets = target_encl.split(i);
			// Call recursively on the two enclosures
			pushSplitTargetEnclosures(initial_enclosures,target_loc,split_sets.first,minTargetCellWidths,target_domain_constraint,use_domain_checking);
			pushSplitTargetEnclosures(initial_enclosures,target_loc,split_sets.second,minTargetCellWidths,target_domain_constraint,use_domain_checking);
		}
	}

	// If we could not split, we put the hybrid enclosure into the new initial_enclosures
	if (!hasSplit)
		initial_enclosures.push_back(EnclosureType(target_loc,target_encl));
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
getHybridMaximumAbsoluteDerivatives(
		const HybridReachabilityAnalyser::SystemType& sys,
		const HybridGridTreeSetPtr& outer_approximation,
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

			// Get the GridTreeSet
			GridTreeSet reach = (*outer_approximation)[loc];
			// For each of its hybrid cells
			for (GridTreeSet::const_iterator cells_it = reach.begin(); cells_it != reach.end(); cells_it++) {
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

	// Returns
	return result;
}


bool
new_reachability_restriction_equals(
		const HybridGridTreeSet& new_restriction,
		const HybridGridTreeSet& old_restriction)
{
	HybridGridTreeSet old_restriction_copy = old_restriction;
	old_restriction_copy.remove(new_restriction);

	return old_restriction_copy.empty();
}

std::list<EnclosureType>
cells_to_smallest_enclosures(
		HybridGridTreeSet cells,
		int maximum_grid_depth)
{
    cells.mince(maximum_grid_depth);

    std::list<EnclosureType> enclosures;
    for (HybridGridTreeSet::const_iterator cell_it = cells.begin(); cell_it != cells.end(); ++cell_it) {
        EnclosureType encl(cell_it->first,cell_it->second.box());
        enclosures.push_back(encl);
    }

    return enclosures;
}


std::list<EnclosureType>
restrict_enclosures(
		const std::list<EnclosureType> enclosures,
		const HybridGridTreeSet& restriction)
{
	std::list<EnclosureType> result;

	for (std::list<EnclosureType>::const_iterator encl_it = enclosures.begin(); encl_it != enclosures.end(); ++encl_it) {
		HybridBox hbox(encl_it->first,encl_it->second.bounding_box());
		if (possibly(restriction.overlaps(hbox)))
			result.push_back(*encl_it);
	}

	return result;
}


HybridGridTreeSet
possibly_feasible_cells(
		const HybridGridTreeSet& reach,
		const HybridConstraintSet& constraint,
		const HybridFloatVector eps,
		HybridGridTreeSet reachability_restriction,
		int accuracy)
{
	reachability_restriction.mince(accuracy);
	HybridGridTreeSet feasible_reachability_restriction = possibly_overlapping_cells(reachability_restriction,constraint);
	HybridVectorFunction constraint_functions = constraint.functions();
	HybridBoxes eps_constraint_codomain = eps_codomain(feasible_reachability_restriction, eps, constraint_functions);

	HybridConstraintSet eps_constraint(constraint_functions,eps_constraint_codomain);

	return possibly_overlapping_cells(reach,eps_constraint);
}


} // namespace Ariadne
