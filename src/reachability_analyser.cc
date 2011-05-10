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

#include "settings.h"
#include "evolver_interface.h"
#include "taylor_calculus.h"

#include "discretiser.h"
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
HybridReachabilityAnalyser(const HybridDiscretiser<HybridEvolver::ContinuousEnclosureType>& discretiser)
    : _settings(new EvolutionSettingsType())
    , _discretiser(discretiser.clone())
	, free_cores(0)
{
}

const CalculusInterface<TaylorModel>&
HybridReachabilityAnalyser::
_getCalculusInterface() const
{
	ImageSetHybridEvolver& evolver = dynamic_cast<ImageSetHybridEvolver&>(*this->_discretiser->evolver());
	return evolver.getCalculusInterface();
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
HybridReachabilityAnalyser::_upper_reach(const HybridAutomaton& sys,
                                         const HybridGridTreeSet& set,
                                         const HybridTime& time,
                                         const int accuracy) const
{
    HybridGrid grid=grid_for(sys,*_settings);
    HybridGridTreeSet result(grid);
    HybridGridTreeSet cells=set;
    cells.mince(accuracy);
    for(HybridGridTreeSet::const_iterator iter=cells.begin(); iter!=cells.end(); ++iter) {
        EnclosureType enclosure=_discretiser->enclosure(*iter);
        result.adjoin(_discretiser->reach(sys,enclosure,time,grid,accuracy,UPPER_SEMANTICS));
    }
    return result;
}


std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridReachabilityAnalyser::_upper_reach_evolve(const HybridAutomaton& sys,
                                                const HybridGridTreeSet& set,
                                                const HybridTime& time,
                                                const int accuracy) const
{
    ARIADNE_LOG(4,"HybridReachabilityAnalyser::_upper_reach_evolve(...)\n");

    HybridGrid grid=grid_for(sys,*_settings);
	GTS initial_set = set;
	initial_set.mince(accuracy);

	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;

	ARIADNE_ASSERT_MSG(concurrency>0 && concurrency <= boost::thread::hardware_concurrency(),"Error: concurrency must be positive and less than the maximum allowed.");

	UpperReachEvolveWorker worker(_discretiser,sys,set,time,grid,accuracy,concurrency);

	std::pair<GTS,GTS> result = worker.get_result();

    ARIADNE_LOG(4,"Reach size = "<<result.first.size()<<"\n");
    ARIADNE_LOG(4,"Final size = "<<result.second.size()<<"\n");
    return result;
}

std::pair<HybridGridTreeSet,HybridGridTreeSet>
HybridReachabilityAnalyser::_upper_reach_evolve_continuous(const HybridAutomaton& sys,
                                                		   const list<EnclosureType>& initial_enclosures,
                                                		   const HybridTime& time,
                                                		   EvolutionDirection direction,
                                                		   int accuracy) const
{
	std::pair<GTS,GTS> result;
	GTS& reach = result.first;
	GTS& evolve = result.second;

	ARIADNE_LOG(4,"Evolving and discretising...\n");

	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;
	ARIADNE_ASSERT_MSG(concurrency>0 && concurrency <= boost::thread::hardware_concurrency(),"Error: concurrency must be positive and less than the maximum allowed.");

	HybridGrid grid=grid_for(sys,*_settings);

	UpperReachEvolveContinuousWorker worker(_discretiser,sys,initial_enclosures,time,direction,grid,accuracy,concurrency);
	result = worker.get_result();

    ARIADNE_LOG(4,"Reach size = "<<reach.size()<<"\n");
    ARIADNE_LOG(4,"Final size = "<<evolve.size()<<"\n");
    return result;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_evolve(const SystemType& system,
             const HybridImageSet& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_evolve(...)\n");

    HybridGrid grid=grid_for(system,*_settings);
    GTS evolve;

    list<EnclosureType> initial_enclosures = split_initial_set(initial_set,grid,_settings->maximum_grid_depth,LOWER_SEMANTICS);

	ARIADNE_LOG(3,"Computing evolution...\n");
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        GTS cell_final=_discretiser->evolve(system,*encl_it,time,grid,_settings->maximum_grid_depth,LOWER_SEMANTICS);
        evolve.adjoin(cell_final);
    }

    return evolve;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_reach(const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach(...)\n");

    HybridGrid grid=grid_for(system,*_settings);
    GTS reach(grid);

    list<EnclosureType> initial_enclosures = split_initial_set(initial_set,grid,_settings->maximum_grid_depth,LOWER_SEMANTICS);

    ARIADNE_LOG(3,"Evolving and discretising...\n");

    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        reach.adjoin(_discretiser->reach(system,*encl_it,time,grid,_settings->maximum_grid_depth,LOWER_SEMANTICS));
    }

	return reach;
}


std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
lower_reach_evolve(const SystemType& system,
                   const HybridImageSet& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach_evolve(...)\n");

    HybridGrid grid=grid_for(system,*_settings);
    GTS reach; GTS evolve;

    list<EnclosureType> initial_enclosures = split_initial_set(initial_set,grid,_settings->maximum_grid_depth,LOWER_SEMANTICS);

	ARIADNE_LOG(3,"Computing evolution...\n");
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        GTS cell_reach,cell_final;
        make_lpair(cell_reach,cell_final)=_discretiser->evolution(system,*encl_it,time,grid,_settings->maximum_grid_depth,LOWER_SEMANTICS);
        reach.adjoin(cell_reach);
        evolve.adjoin(cell_final);
    }

    return make_pair(reach,evolve);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_evolve(const SystemType& system,
             const HybridImageSet& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_evolve(system,set,time)\n");
 
	GTS reach, evolve;
	make_lpair(reach,evolve) = upper_reach_evolve(system, initial_set, time); // Runs the upper_reach_evolve routine on its behalf

    return evolve;
}



HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_reach(const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach(system,set,time)\n");

	GTS reach, evolve;
	make_lpair(reach,evolve) = upper_reach_evolve(system, initial_set, time);

    return reach;
}


std::pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
upper_reach_evolve(const SystemType& system,
                   const HybridImageSet& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach_evolve(system,set,time)\n");
    ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");

    HybridGrid grid=grid_for(system,*_settings);
    GTS found(grid),evolve(grid),reach(grid);
    int maximum_grid_depth = _settings->maximum_grid_depth;
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
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");

    list<EnclosureType> initial_enclosures = split_initial_set(initial_set,grid,maximum_grid_depth,UPPER_SEMANTICS);

	ARIADNE_LOG(3,"Computing initial evolution...\n");
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        GTS cell_reach,cell_final;
        make_lpair(cell_reach,cell_final)=_discretiser->evolution(system,*encl_it,hybrid_lock_to_grid_time,grid,maximum_grid_depth,UPPER_SEMANTICS);
        reach.adjoin(cell_reach);
        evolve.adjoin(cell_final);
    }

    for(uint i=1; i<time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        make_lpair(found,evolve) = _upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(4,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(4,"evolve.size()="<<evolve.size()<<"\n");
        reach.adjoin(found);
        ARIADNE_LOG(3,"found "<<found.size()<<" cells.\n");
    }

    ARIADNE_LOG(5,"remaining_time="<<remaining_time<<"\n");
    if(!evolve.empty() && remaining_time > 0) {
        ARIADNE_LOG(3,"computing evolution for the remaining time...\n");
        make_lpair(found,evolve) = _upper_reach_evolve(system,evolve,hybrid_remaining_time,maximum_grid_depth);
        reach.adjoin(found);
    }

    reach.recombine();
    ARIADNE_LOG(3,"reach="<<reach<<"\n");
    evolve.recombine();
    ARIADNE_LOG(3,"evolve="<<evolve<<"\n");
    return std::make_pair(reach,evolve);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
chain_reach(const SystemType& system,
            const HybridImageSet& initial_set) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::chain_reach(system,initial_set)\n");

    HybridBoxes domain = _settings->domain_bounds;
    Float transient_time = _settings->transient_time;
    int transient_steps = _settings->transient_steps;
    Float lock_to_grid_time = _settings->lock_to_grid_time;
    int lock_to_grid_steps = _settings->lock_to_grid_steps;
    int maximum_grid_depth = _settings->maximum_grid_depth;

    ARIADNE_LOG(4,"transient_time=("<<transient_time<<","<<transient_steps<<")\n");
    ARIADNE_LOG(4,"lock_to_grid_time=("<<lock_to_grid_time<<","<<lock_to_grid_steps<<")\n");

	// Checks consistency of the bounding domain in respect to the state space
	HybridSpace hspace = system.state_space();
	for (HybridSpace::locations_const_iterator hs_it = hspace.locations_begin(); hs_it != hspace.locations_end(); ++hs_it) {
		if (domain.find(hs_it->first) == domain.end()) {
			ARIADNE_FAIL_MSG("Error: the system state space and the bounding domain space do not match on the discrete space."); }		
		else if (hs_it->second != domain[hs_it->first].size()) {
			ARIADNE_FAIL_MSG("Error: the system state space and the bounding domain space do not match on the continuous space."); }}

    HybridGrid grid=grid_for(system,*_settings);
    GTS bounding(grid), evolve(grid), reach(grid), found(grid), intermediate(grid);

    bounding.adjoin_outer_approximation(domain,0u); 
    ARIADNE_LOG(4,"bounding_size(pre recombine)="<<bounding.size()<<"\n");
	bounding.recombine();
	HybridBoxes bounding_box = bounding.bounding_box(); // Used for the restriction check
    ARIADNE_LOG(4,"bounding_size(post recombine)="<<bounding.size()<<"\n");

    if(transient_time <= 0.0 || transient_steps <= 0) {
        transient_time = lock_to_grid_time;
        transient_steps = lock_to_grid_steps; }

    HybridTime hybrid_transient_time(transient_time, transient_steps);
    ARIADNE_LOG(3,"Computing first evolution step...\n");

    list<EnclosureType> initial_enclosures = split_initial_set(initial_set,grid,maximum_grid_depth,UPPER_SEMANTICS);

	ARIADNE_LOG(3,"Computing transient evolution...\n");
    for (list<EnclosureType>::const_iterator encl_it = initial_enclosures.begin(); encl_it != initial_enclosures.end(); encl_it++) {
        GTS cell_reach,cell_final;
        make_lpair(cell_reach,cell_final)=_discretiser->evolution(system,*encl_it,hybrid_transient_time,grid,maximum_grid_depth,UPPER_SEMANTICS);
        reach.adjoin(cell_reach);
        evolve.adjoin(cell_final);
    }

    evolve.restrict(bounding);

    ARIADNE_LOG(3,"found "<<reach.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
    
    ARIADNE_LOG(3,"Computing recurrent evolution...\n");
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);

	// While the final set has new cells in respect to the previous reach set, process them and increase the number of locks (starting from 1 due to the initial transient phase)
    for (uint i = 0; !evolve.empty(); i++) {
    	ARIADNE_LOG(3,"Iteration "<<i<<"\n");

		intermediate.adjoin(evolve);  

        make_lpair(found,evolve)=_upper_reach_evolve(system,evolve,hybrid_lock_to_grid_time,maximum_grid_depth);
        ARIADNE_LOG(4,"found.size()="<<found.size()<<"\n");
        ARIADNE_LOG(4,"evolve.size()="<<evolve.size()<<"\n");

        evolve.remove(intermediate); // Remove only the cells that are intermediate
        evolve.restrict(bounding);
        reach.adjoin(found);

        ARIADNE_LOG(4,"found "<<found.size()<<" cells, of which "<<evolve.size()<<" are new.\n");
    }
    reach.recombine();

    reach.restrict(bounding);

    return reach;
}

HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
chain_reach(
		const SystemType& system,
        const HybridImageSet& initial_set,
        const HybridBoxes& bounding_set) const
{
	_settings->domain_bounds = bounding_set;

	return chain_reach(system,initial_set);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
_outer_chain_reach_splitted(
		const SystemType& system,
		const std::list<EnclosureType>& initial_enclosures,
		EvolutionDirection direction,
		const HybridGridTreeSet& restriction) const
{
	ARIADNE_ASSERT_MSG(!(direction == DIRECTION_BACKWARD && restriction.empty()),
			"The restriction must be nonempty if backward reachability is to be performed.");

    const Float& lock_to_grid_time = _settings->lock_to_grid_time;
    const int& lock_to_grid_steps = _settings->lock_to_grid_steps;
    const int& maximum_grid_depth = _settings->maximum_grid_depth;

    HybridGrid grid=grid_for(system,*_settings);
    HybridGridTreeSet new_final(grid), new_reach(grid), reach(grid), final(grid);

    bool use_domain_checking = restriction.empty();

    list<EnclosureType> working_enclosures = initial_enclosures;

    ARIADNE_LOG(3,"Computing recurrent " << (direction == DIRECTION_FORWARD ? "forward" : "backward") << " evolution...\n");
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);

    uint i=0;
    while (!working_enclosures.empty())
	{
    	ARIADNE_LOG(3,"Iteration " << i++ << "\n");

        ARIADNE_LOG(4,"Initial enclosures size = " << working_enclosures.size() << "\n");

        make_lpair(new_reach,new_final) = _upper_reach_evolve_continuous(system,working_enclosures,hybrid_lock_to_grid_time,direction,maximum_grid_depth);

        new_final.remove(final);
		new_reach.remove(reach);
	    ARIADNE_LOG(4,"Reach size after removal = "<<new_reach.size()<<"\n");
	    ARIADNE_LOG(4,"Final size after removal = "<<new_final.size()<<"\n");

	    if (!use_domain_checking) {
	    	new_final.restrict(restriction);
	    	new_reach.restrict(restriction);
		    ARIADNE_LOG(4,"Reach size after restricting = "<<new_reach.size()<<"\n");
		    ARIADNE_LOG(4,"Final size after restricting = "<<new_final.size()<<"\n");
	    }

		new_reach.mince(maximum_grid_depth);
		new_final.mince(maximum_grid_depth);
		ARIADNE_LOG(4,"Reach size after mincing = "<<new_reach.size()<<"\n");
		ARIADNE_LOG(4,"Final size after mincing = "<<new_final.size()<<"\n");

        working_enclosures.clear();

        if (direction == DIRECTION_FORWARD)
        	_outer_chain_reach_forward_pushTargetCells(new_reach,system,working_enclosures,use_domain_checking);
        else
        	_outer_chain_reach_backward_pushSourceCells(new_reach,restriction,system,working_enclosures);

		_outer_chain_reach_pushLocalFinalCells(new_final,working_enclosures,use_domain_checking);

        reach.adjoin(new_reach);
		final.adjoin(new_final);

		reach.recombine();
		final.recombine();
    }

	ARIADNE_LOG(3,"Found a total of " << reach.size() << " reached cells.\n");

    return reach;
}


void
HybridReachabilityAnalyser::
_outer_chain_reach_forward_pushTargetCells(const HybridGridTreeSet& reachCells,
										   const SystemType& system,
										   std::list<EnclosureType>& result_enclosures,
										   bool use_domain_checking) const
{
	for (HybridGridTreeSet::const_iterator cell_it = reachCells.begin(); cell_it != reachCells.end(); cell_it++)
	{
		const DiscreteState& loc = cell_it->first;
		const Box& bx = cell_it->second.box();
		const Box& domain = _settings->domain_bounds[loc];

		ARIADNE_LOG(5,"Checking box "<< bx <<" in location " << loc.name() << "\n");

		if (use_domain_checking && bx.disjoint(domain))
			throw ReachOutOfDomainException("a reach enclosure is outside the domain");

		if (!_outer_chain_reach_isOutsideInvariants(bx,system.mode(loc).invariants())) {
			_outer_chain_reach_forward_pushTargetEnclosures(system.transitions(loc),ContinuousEnclosureType(bx),
					system.mode(loc).dynamic(),*_settings->grid,result_enclosures,use_domain_checking);
		}
	}
}


void
HybridReachabilityAnalyser::
_outer_chain_reach_backward_pushSourceCells(
		const HybridGridTreeSet& targetCells,
		HybridGridTreeSet sourceCellsOverapprox,
		const SystemType& system,
		std::list<EnclosureType>& result_enclosures) const
{
	// Adopt the minimum granularity when checking the cells
	sourceCellsOverapprox.mince(_settings->maximum_grid_depth);

	for (HybridGridTreeSet::const_iterator cell_it = sourceCellsOverapprox.begin(); cell_it != sourceCellsOverapprox.end(); ++cell_it) {
		const DiscreteState& loc = cell_it->first;
		const Box& bx = cell_it->second.box();

		ARIADNE_LOG(5,"Checking box "<< bx <<" in location " << loc.name() << "\n");

		if (!_outer_chain_reach_isOutsideInvariants(bx,system.mode(loc).invariants()))
			_outer_chain_reach_backward_pushSourceEnclosures(system.transitions(loc),loc,bx,targetCells,
					system.mode(loc).dynamic(),*_settings->grid,result_enclosures);
	}
}

bool
HybridReachabilityAnalyser::
_outer_chain_reach_isOutsideInvariants(const Box& bx,
									   const std::map<DiscreteEvent,VectorFunction>& invariants) const
{
	ARIADNE_LOG(6,"Checking invariants...\n");

	for (std::map<DiscreteEvent,VectorFunction>::const_iterator inv_it = invariants.begin(); inv_it != invariants.end(); ++inv_it) {
		const VectorFunction& activation = inv_it->second;
		tribool is_active = _getCalculusInterface().active(activation,bx);
		if (definitely(is_active)) {
			ARIADNE_LOG(6,"Invariant " << *inv_it << " is definitely active: transitions will not be checked.\n");
			return true;
		}
	}

	return false;
}

void
HybridReachabilityAnalyser::
_outer_chain_reach_forward_pushTargetEnclosures(
		const std::list<DiscreteTransition>& transitions,
		const ContinuousEnclosureType& source,
		const VectorFunction& dynamic,
		const HybridGrid& grid,
		std::list<EnclosureType>& result_enclosures,
		bool use_domain_checking) const
{
    long numCellDivisions = (1<<_settings->maximum_grid_depth);

	ARIADNE_LOG(6,"Checking transitions...\n");

	for (list<DiscreteTransition>::const_iterator trans_it = transitions.begin(); trans_it != transitions.end(); trans_it++)
	{
		ARIADNE_LOG(7,"Target: "<<trans_it->target()<<", Forced?: " << trans_it->forced() << "\n");

		if (_is_transition_feasible(*trans_it,dynamic,source)) {
			const DiscreteState& target_loc = trans_it->target();
			const ContinuousEnclosureType target_encl = _getCalculusInterface().reset_step(trans_it->reset(),source);
			const Box& target_bounding = _settings->domain_bounds[target_loc];
			const Vector<Float> minTargetCellWidths = grid[target_loc].lengths()/numCellDivisions;

			pushSplitTargetEnclosures(result_enclosures,target_loc,target_encl,minTargetCellWidths,target_bounding,use_domain_checking);
		}
	}
}

void
HybridReachabilityAnalyser::
_outer_chain_reach_backward_pushSourceEnclosures(
		const std::list<DiscreteTransition>& transitions,
		const DiscreteState& sourceLocation,
		const Box& sourceBox,
		const HybridGridTreeSet& targetCells,
		const VectorFunction& dynamic,
		const HybridGrid& grid,
		std::list<EnclosureType>& result_enclosures) const
{
	ARIADNE_LOG(6,"Checking transitions...\n");

	ContinuousEnclosureType sourceEncl(sourceBox);

	// In order to add a source hybrid box, just one possibly overlapping target enclosure is sufficient
	for (list<DiscreteTransition>::const_iterator trans_it = transitions.begin(); trans_it != transitions.end(); trans_it++)
	{
		ARIADNE_LOG(7,"Target: "<<trans_it->target()<<", Forced?: " << trans_it->forced() << "\n");

		if (_is_transition_feasible(*trans_it,dynamic,sourceEncl)) {
			const DiscreteState& target_loc = trans_it->target();
			const ContinuousEnclosureType target_encl = _getCalculusInterface().reset_step(trans_it->reset(),sourceEncl);
			const HybridBox targetHBox(target_loc,target_encl.bounding_box());

			if (possibly(targetCells.overlaps(targetHBox))) {
				ARIADNE_LOG(7,"Target enclosure overlaps with target cells: adding the source enclosure.");
				result_enclosures.push_back(EnclosureType(sourceLocation,sourceEncl));
				break;
			}
		}
	}
}

bool
HybridReachabilityAnalyser::
_is_transition_feasible(
		const DiscreteTransition& trans,
		const VectorFunction& dynamic,
		const ContinuousEnclosureType& source) const
{
	bool result = false;

	const VectorFunction& activation = trans.activation();
	const bool is_forced = trans.forced();

	tribool is_guard_active = _getCalculusInterface().active(activation,source);

	ARIADNE_LOG(8,"Guard activity: " << is_guard_active << "\n");

	/*
	 * a) If the guard is definitely active and the transition is forced, then we are definitely outside the related invariant
	 * b) If the transition is not forced, it suffices to have a possibly active guard: we then must perform the transition
	 * c) If the transition is forced and the guard is only possibly active, we check the crossing:
	 *    i) If it is negative, then no transition is possible
	 *	 ii) If it is possibly positive, then we must take the transition
	 */

	if (definitely(is_guard_active) && is_forced) {
		ARIADNE_LOG(8,"Definitely active and forced: the set is outside the implicit invariant of the transition, infeasible.\n");
		result = false;
	} else if (possibly(is_guard_active) && !is_forced) {
		ARIADNE_LOG(8,"Possibly active and permissive: feasible.\n");
		result = true;
	} else if (possibly(is_guard_active) && is_forced) {
		ARIADNE_LOG(8,"Possibly active and forced: checking whether the crossing is nonnegative...\n");
		tribool positive_crossing = positively_crossing(source.bounding_box(),dynamic,activation[0]);
		if (possibly(positive_crossing)) {
			ARIADNE_LOG(8,"Possibly positive: feasible.\n");
			result = true;
		} else {
			ARIADNE_LOG(8,"Negative: infeasible.\n");
			result = false;
		}
	} else {
		ARIADNE_LOG(8,"Inactive: infeasible.\n");
		result = false;
	}

	return result;
}

void
HybridReachabilityAnalyser::
_outer_chain_reach_pushLocalFinalCells(const HybridGridTreeSet& finalCells,
									   std::list<EnclosureType>& result_enclosures,
									   bool use_domain_checking) const
{
	for (GTS::const_iterator cell_it = finalCells.begin(); cell_it != finalCells.end(); ++cell_it) {
		const DiscreteState& loc = cell_it->first;
		const Box& domain = _settings->domain_bounds[loc];
		const Box& bx = cell_it->second.box();

		if (use_domain_checking && bx.disjoint(domain)) {
			ARIADNE_LOG(5,"Discarding enclosure " << bx << " from final cell outside the domain in location " << loc.name() <<".\n");
			throw ReachOutOfDomainException("a final cell is outside the domain");
		} else {
			result_enclosures.push_back(_discretiser->enclosure(*cell_it));
		}
	}
}


bool
HybridReachabilityAnalyser::
fb_refinement_check(
		SystemType& system,
		const HybridImageSet& initial_set,
		const HybridBoxes& target_bounds,
		const HybridGridTreeSet& reachability)
{
	bool result = false;

	HybridGridTreeSet forward_reach = reachability;
	HybridGridTreeSet backward_reach(*_settings->grid);

	int& maximum_grid_depth = _settings->maximum_grid_depth;
	int initial_maximum_grid_depth = maximum_grid_depth;

	time_t mytime;
	time(&mytime);
	string foldername = system.name()+"-fbr-png";
	mkdir(foldername.c_str(),0777);
	string timestring = asctime(localtime(&mytime));
	timestring.erase(std::remove(timestring.begin(), timestring.end(), '\n'), timestring.end());
	foldername = foldername+"/"+timestring;
	mkdir(foldername.c_str(),0777);
	string plot_dirpath = foldername;

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

		if (maximum_grid_depth >= _settings->highest_maximum_grid_depth) {
			result = false;
			break;
		}

		tribool backward_result = _fb_refinement_step_check(system,initial_set,target_bounds,forward_reach,
				backward_reach,DIRECTION_BACKWARD,plot_dirpath);

		if (!indeterminate(backward_result)) {
			result = backward_result;
			break;
		}

		if (maximum_grid_depth >= _settings->highest_maximum_grid_depth) {
			result = false;
			break;
		}

		tribool forward_result = _fb_refinement_step_check(system,initial_set,target_bounds,backward_reach,
				forward_reach,DIRECTION_FORWARD,plot_dirpath);

		if (!indeterminate(forward_result)) {
			result = forward_result;
			break;
		}
	}

	_settings->maximum_grid_depth = initial_maximum_grid_depth;

	return result;
}


tribool
HybridReachabilityAnalyser::
_fb_refinement_step_check(
		SystemType& system,
		const HybridImageSet& initial_set,
		const HybridBoxes& target_bounds,
		HybridGridTreeSet& old_restriction,
		HybridGridTreeSet& new_restriction,
		EvolutionDirection direction,
		string plot_dirpath)
{
	// TODO: Completely remove the plotting calls.
	bool plot = true;

	int& maximum_grid_depth = _settings->maximum_grid_depth;

	_fb_refine_settings();

	ARIADNE_LOG(1,"Depth " << maximum_grid_depth << "\n");
	old_restriction.mince(maximum_grid_depth);

	HybridGridTreeSet output(*_settings->grid);
	std::list<EnclosureType> starting_enclosures;
	for (HybridGridTreeSet::const_iterator cell_it = old_restriction.begin(); cell_it != old_restriction.end(); ++cell_it) {
		const DiscreteState& loc = cell_it->first;
		const Box& cell_box = cell_it->second.box();
		HybridBox hcell_box(loc,cell_box);

		bool add_enclosure = (direction == DIRECTION_FORWARD ?
				possibly(initial_set.overlaps(hcell_box)) : !Ariadne::covers(target_bounds,hcell_box));

		if (add_enclosure) {
			starting_enclosures.push_back(EnclosureType(loc,cell_box));
			if (plot) output.adjoin(*cell_it);
		}
	}

	ARIADNE_LOG(1,"Starting enclosures size: " << starting_enclosures.size() << "\n");

	if (plot) _plot_reach(output,plot_dirpath,(direction == DIRECTION_FORWARD ? "initial" : "final"));

	if (starting_enclosures.empty())
		return true;

	ARIADNE_LOG(1,"Performing " << (direction == DIRECTION_FORWARD ? "forward" : "backward") << " refinement...\n");

	new_restriction = _outer_chain_reach(system,starting_enclosures,direction,old_restriction);

	if (plot) _plot_reach(new_restriction,plot_dirpath,(direction == DIRECTION_FORWARD ? "forward" : "backward"));

	if (new_reachability_restriction_equals(new_restriction,old_restriction))
		return false;

	return indeterminate;
}

void
HybridReachabilityAnalyser::
_fb_refine_settings()
{
	_settings->maximum_grid_depth++;
	_discretiser->settings().maximum_enclosure_cell = getMaximumEnclosureCell(*_settings->grid,_settings->maximum_grid_depth);
	for (std::map<DiscreteState,Float>::iterator step_it = _discretiser->settings().hybrid_maximum_step_size.begin();
												 step_it != _discretiser->settings().hybrid_maximum_step_size.end();
												 ++step_it) {
		step_it->second /= 2;
	}
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
outer_chain_reach(
		SystemType& system,
		const HybridImageSet& initial_set,
		EvolutionDirection direction) const
{
	const std::list<EnclosureType> initial_enclosures = split_initial_set(initial_set,*_settings->grid,_settings->maximum_grid_depth,UPPER_SEMANTICS);

	return _outer_chain_reach(system,initial_enclosures,direction,_settings->reachability_restriction);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
_outer_chain_reach(
		SystemType& system,
		const std::list<EnclosureType>& initial_enclosures,
		EvolutionDirection direction,
		const HybridGridTreeSet& reachability_restriction) const
{
	HybridBoxes feasibility_bounds = unbounded_hybrid_boxes(system.state_space());

	// Uses dummy arguments
	return _outer_chain_reach(system,initial_enclosures,direction,reachability_restriction,false,feasibility_bounds,NOT_INSIDE_TARGET);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
outer_chain_reach(
		SystemType& system,
		const HybridImageSet& initial_set,
		EvolutionDirection direction,
		bool check_target_bounds_for_skipping,
		const HybridBoxes& target_bounds_for_skipping,
		OuterReachabilitySkippingPolicy skippingPolicy) const
{
	const std::list<EnclosureType> initial_enclosures = split_initial_set(initial_set,*_settings->grid,_settings->maximum_grid_depth,UPPER_SEMANTICS);

	return _outer_chain_reach(system,initial_enclosures,direction,_settings->reachability_restriction,check_target_bounds_for_skipping,target_bounds_for_skipping,skippingPolicy);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
_outer_chain_reach(
		SystemType& system,
		const std::list<EnclosureType>& initial_enclosures,
		EvolutionDirection direction,
		const HybridGridTreeSet& reachability_restriction,
		bool check_target_bounds_for_skipping,
		const HybridBoxes& target_bounds_for_skipping,
		OuterReachabilitySkippingPolicy skippingPolicy) const
{
	HybridGridTreeSet reach;

	RealConstantSet original_constants = system.nonsingleton_accessible_constants();

	std::list<RealConstantSet> split_intervals_set = _getSplitConstantsIntervalsSet(system,_settings->splitting_constants_target_ratio);

	try {
		uint i = 0;
		for (std::list<RealConstantSet>::const_iterator set_it = split_intervals_set.begin(); set_it != split_intervals_set.end(); ++set_it) {
			ARIADNE_LOG(2,"<Split constants set #" << ++i << " : " << *set_it << " >\n");

			system.substitute(*set_it);

			HybridGridTreeSet local_reach = _outer_chain_reach_splitted(system,initial_enclosures,direction,reachability_restriction);

			reach.adjoin(local_reach);

			if (check_target_bounds_for_skipping) {
				if (skippingPolicy == NOT_INSIDE_TARGET && definitely(!local_reach.subset(target_bounds_for_skipping))) {
					system.substitute(original_constants);
					throw ReachOutOfTargetException("The reached set is not inside the target region.");
				} else if (skippingPolicy == SUPERSET_OF_TARGET && definitely(superset(local_reach.bounding_box(),target_bounds_for_skipping))) {
					system.substitute(original_constants);
					throw ReachEnclosesTargetException("The reached set encloses the target region.");
				}
			}
		}
	} catch (ReachOutOfDomainException ex) {
		system.substitute(original_constants);
		throw ReachOutOfDomainException(ex.what());
	}

	system.substitute(original_constants);

	ARIADNE_ASSERT_MSG(!reach.empty(),"The outer chain reachability of " << system.name() << " is empty: check the initial set.");

	return reach;
}


std::pair<HybridReachabilityAnalyser::SetApproximationType,DisproveData>
HybridReachabilityAnalyser::
_lower_chain_reach(const SystemType& system,
				  const HybridImageSet& initial_set,
				  const HybridBoxes& feasibility_region,
				  bool terminate_as_soon_as_infeasible) const
{
	typedef std::list<EnclosureType> EL;
	typedef std::map<DiscreteState,uint> HUM;

	const uint concurrency = boost::thread::hardware_concurrency() - free_cores;
	ARIADNE_ASSERT_MSG(concurrency>0 && concurrency <= boost::thread::hardware_concurrency(),"Error: concurrency must be positive and less than the maximum allowed.");

	HybridGrid grid = *_settings->grid;
    GTS reach(grid);

	TimeType lock_time(_settings->lock_to_grid_time,_settings->lock_to_grid_steps);

	bool use_domain_checking = _settings->reachability_restriction.empty();

    DisproveData globalFalsInfo(system.state_space());

    EL initial_enclosures = split_initial_set(initial_set,grid,_settings->maximum_grid_depth,LOWER_SEMANTICS);

    ARIADNE_LOG(3,"Computing recurrent evolution...\n");

    uint i=0;
    while (!initial_enclosures.empty()) {
		ARIADNE_LOG(2,"Iteration " << i++ << "\n");
		EL final_enclosures;
		std::pair<HUM,HUM> evolve_sizes;

		HUM& adjoined_evolve_sizes = evolve_sizes.first;
		HUM& superposed_evolve_sizes = evolve_sizes.second;

		GTS new_reach;
		DisproveData newFalsInfo(system.state_space());

		ARIADNE_LOG(4,"Initial enclosures size = " << initial_enclosures.size() << "\n");

		LowerChainReachWorker worker(_discretiser,initial_enclosures,system,lock_time,grid,
									 _settings->maximum_grid_depth,concurrency, feasibility_region, terminate_as_soon_as_infeasible);

		ARIADNE_LOG(4,"Evolving and discretising...\n");

		make_ltuple<std::pair<HUM,HUM>,EL,GTS,DisproveData>(evolve_sizes,final_enclosures,
																new_reach,newFalsInfo) = worker.get_result();

		globalFalsInfo.updateWith(newFalsInfo);

		ARIADNE_LOG(4,"Reach size before removal = " << new_reach.size() << "\n");

		new_reach.remove(reach);
		ARIADNE_LOG(4,"Reach size after removal  = " << new_reach.size() << "\n");
		if (new_reach.empty())
			break;

		reach.adjoin(new_reach);

		ARIADNE_LOG(4,"Final enclosures size = " << final_enclosures.size() << "\n");

		_filter_enclosures(final_enclosures,initial_enclosures,adjoined_evolve_sizes,superposed_evolve_sizes,use_domain_checking);
	}

	return make_pair<GTS,DisproveData>(reach,globalFalsInfo);
}

void
HybridReachabilityAnalyser::
_filter_enclosures(
		std::list<EnclosureType>& final_enclosures,
		std::list<EnclosureType>& initial_enclosures,
		const std::map<DiscreteState,uint>& adjoined_evolve_sizes,
		const std::map<DiscreteState,uint>& superposed_evolve_sizes,
		bool use_domain_checking) const
{
	while (!final_enclosures.empty()) {
		EnclosureType encl = final_enclosures.front();
		final_enclosures.pop_front();

		const DiscreteState& loc = encl.location();
		const Box& encl_box = encl.continuous_state_set().bounding_box();

		if (use_domain_checking && encl_box.disjoint(_settings->domain_bounds[loc])) {
			ARIADNE_FAIL_MSG("Found an enclosure in location " << loc.name() << " with bounding box " << encl.continuous_state_set().bounding_box() <<
							 " lying outside the domain in lower semantics: the domain is incorrect.\n");
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

std::pair<HybridReachabilityAnalyser::SetApproximationType,DisproveData>
HybridReachabilityAnalyser::
lower_chain_reach(
		SystemType& system,
		const HybridImageSet& initial_set) const
{
	HybridBoxes safe_region = unbounded_hybrid_boxes(system.state_space());

	return lower_chain_reach(system,initial_set,safe_region,false);
}

std::pair<HybridReachabilityAnalyser::SetApproximationType,DisproveData>
HybridReachabilityAnalyser::
lower_chain_reach(
		SystemType& system,
		const HybridImageSet& initial_set,
		const HybridBoxes& safe_region,
		bool terminate_as_soon_as_disproved) const
{
	HybridGridTreeSet reach(*_settings->grid);
	DisproveData disproveData(system.state_space());

	RealConstantSet original_constants = system.nonsingleton_accessible_constants();

	std::list<RealConstantSet> split_intervals_set = _getSplitConstantsIntervalsSet(system,_settings->splitting_constants_target_ratio);
	std::list<RealConstantSet> split_midpoints_set = getSplitConstantsMidpointsSet(split_intervals_set);

	uint i = 0;
	// Progressively adds the results for each subsystem
	for (std::list<RealConstantSet>::const_iterator set_it = split_midpoints_set.begin();
												    set_it != split_midpoints_set.end();
												    ++set_it) {
		ARIADNE_LOG(2,"<Split constants set #" << ++i << " : " << *set_it << " >\n");

		system.substitute(*set_it);

		std::pair<HybridGridTreeSet,DisproveData> reachAndDisproveData =
				_lower_chain_reach(system,initial_set,safe_region,terminate_as_soon_as_disproved);

		reach.adjoin(reachAndDisproveData.first);
		disproveData.updateWith(reachAndDisproveData.second);

		ARIADNE_LOG(3,"Disprove data: " << reachAndDisproveData.second << "\n");

		if (terminate_as_soon_as_disproved && reachAndDisproveData.second.getIsDisproved())
			break;
	}
	system.substitute(original_constants);

	return std::pair<HybridGridTreeSet,DisproveData>(reach,disproveData);
}


void
HybridReachabilityAnalyser::
tuneEvolverSettings(
		const SystemType& system,
		const HybridFloatVector& hmad,
		uint maximum_grid_depth,
		Semantics semantics)
{
	HybridGrid grid = *_settings->grid;
	_discretiser->settings().maximum_enclosure_cell = getMaximumEnclosureCell(grid,maximum_grid_depth);
	ARIADNE_LOG(2, "Maximum enclosure cell: " << _discretiser->settings().maximum_enclosure_cell << "\n");
	_discretiser->settings().hybrid_maximum_step_size = getHybridMaximumStepSize(hmad,grid,maximum_grid_depth,semantics);
	ARIADNE_LOG(2, "Maximum step size: " << _discretiser->settings().hybrid_maximum_step_size << "\n");
}

std::list<RealConstantSet>
HybridReachabilityAnalyser::
_getSplitConstantsIntervalsSet(
		HybridAutomaton system,
		float tolerance) const
{
	const RealConstantSet& locked_constants = _settings->locked_constants;
	const HybridBoxes& domain = _settings->domain_bounds;

	const RealConstantSet original_constants = system.accessible_constants();

	const RealConstantIntMap split_factors = getSplitFactorsOfConstants(system,locked_constants,tolerance,domain);

	std::list<RealConstantSet> result;

	if (split_factors.empty()) {
		result.push_back(original_constants);
		return result;
	}

	// Creates a vector for all the interval splits (i.e. a jagged matrix)
	std::vector<std::vector<RealConstant> > split_intervals_set(split_factors.size());
	uint i=0;
	for (RealConstantIntMap::const_iterator factor_it = split_factors.begin();
											factor_it != split_factors.end();
											++factor_it)
		split_intervals_set[i++] = split(factor_it->first, factor_it->second);

	// Generates all the possible split combinations
	RealConstantSet initial_combination;
	std::vector<std::vector<RealConstant> >::iterator initial_col_it = split_intervals_set.begin();
	std::vector<RealConstant>::iterator initial_row_it = initial_col_it->begin();
	fillSplitSet(split_intervals_set,initial_col_it,initial_row_it,initial_combination,result);

	ARIADNE_LOG(3,"<Split factors: " << split_factors << ", size: " << result.size() << ">\n");

	return result;
}

list<HybridBasicSet<TaylorSet> >
split_initial_set(
		const HybridImageSet initial_set,
		const HybridGrid grid,
		int maximum_grid_depth,
		Semantics semantics)
{
	list<EnclosureType> result;

	// For each location, split the original enclosure into a list of enclosures
	// having widths lesser than the grid lengths at maximum depth
	for(HybridImageSet::locations_const_iterator loc_iter=initial_set.locations_begin();
		 loc_iter!=initial_set.locations_end(); ++loc_iter)
		{
		// Get the grid lengths at the maximum depth (i.e. the minimum cell widths)
		Vector<Float> mincell_widths = grid[loc_iter->first].lengths()/(1 << maximum_grid_depth);

		// Keep a list of the continuous enclosures to split
		list<ContinuousEnclosureType> continuous_enclosures;
		// Insert the original continuous enclosure
		continuous_enclosures.push_front(ContinuousEnclosureType(loc_iter->second));

		/* While there are continuous enclosures:
		 * a) pick one, pop it and check it against the mincell_lengths
		 * 	 i) if larger, split it and put the couple into the continuous_enclosures
		 *   ii) if not, put the enclosure (paired with the location) into the result
		 */
		while (!continuous_enclosures.empty()) {
			ContinuousEnclosureType enclosure = continuous_enclosures.back();
			continuous_enclosures.pop_back();

			bool hasBeenSplit = false;
			for (uint dim=0; dim < enclosure.dimension(); ++dim) {
				if (enclosure[dim].range().width() > mincell_widths[dim]) {
					pair<ContinuousEnclosureType,ContinuousEnclosureType> split_result = enclosure.split(dim);
					continuous_enclosures.push_front(split_result.first);
					continuous_enclosures.push_front(split_result.second);
					// Set the flag as true and skip the remaining checks
					hasBeenSplit = true;
					break;
				}
			}
			// If it has not been split, put the hybrid enclosure into the initial enclosures
			if (!hasBeenSplit)
			{
				// For the case of lower semantics, we want to reduce the enclosure to its centre in order to minimize the radius
				if (semantics == LOWER_SEMANTICS)
					enclosure = ContinuousEnclosureType(Vector<Interval>(enclosure.centre()));
				result.push_back(EnclosureType(loc_iter->first,enclosure));
			}
		}
	}

	// Returns
	return result;
}



RealConstantIntMap
getSplitFactorsOfConstants(
		SystemType& system,
		const RealConstantSet& locked_constants,
		const Float& targetRatioPerc,
		const HybridBoxes& bounding_domain)
{
	/*! Procedure:
	 * 1) Get the derivatives bounds for the system with constants set as their intervals
	 * 2) Get the derivatives bounds for the system with constants set as their midpoints
	 * 3) Get the maximum ratio from this bounds sets and set the target ratio as a fraction of this value
	 * 4) While the ratio is greater than the target ratio:
 	 *	 a) For each non-singleton accessible constant
 	 *	   i) Get the derivatives bounds for the system where the range of the chosen constant has been halved
 	 *	   ii) Update the best ratio
 	 *	 b) Increase the factor of the constant having the best ratio and update the system accordingly by halving the range of the constant
	 */

	// A lower threshold for the maximum ratio, in order to ignore cases where the ratios are too low: in that case,
	// splitting would have no effect and the procedure would not converge
	const float& maxRatioThreshold = 1e-8;

	RealConstantIntMap result;

	RealConstantSet working_constants = system.nonsingleton_accessible_constants();

	// Remove the locked constants
	for (RealConstantSet::const_iterator locked_constant_it = locked_constants.begin();
										 locked_constant_it != locked_constants.end();
										 ++locked_constant_it) {
		RealConstantSet::iterator original_constant_it = working_constants.find(*locked_constant_it);

		if (original_constant_it != working_constants.end())
			working_constants.erase(*original_constant_it);
	}

	if (working_constants.empty())
		return result;

	// Initializes the result and sets the system with the midpoints of the corresponding intervals
	for (RealConstantSet::const_iterator constant_it = working_constants.begin();
												 constant_it != working_constants.end();
												 ++constant_it) {
			result.insert(std::pair<RealConstant,int>(*constant_it,1));
			system.substitute(*constant_it,constant_it->value().midpoint());
	}

	// Gets the derivative widths corresponding to all accessible constants having midpoint value
	HybridFloatVector mid_der_widths = getDerivativeWidths(system, bounding_domain);
	// Restores the system to the original values
	system.substitute(working_constants);

	// While the ratio is sufficiently high, gets the best constant and substitutes half its interval into the system
	// If the maximum ratio is zero, then no constant affects the derivatives and splitting them would neither be necessary nor correct
	// given the current unfair implementation of _getBestConstantToSplit
	Float maxRatio = getMaxDerivativeWidthRatio(system, mid_der_widths, bounding_domain);
	if (maxRatio > maxRatioThreshold) {
		Float ratio = maxRatio;

		while (ratio > targetRatioPerc*maxRatio) {
			RealConstant bestConstant = getBestConstantToSplit(system, working_constants, mid_der_widths, bounding_domain);
			Interval originalInterval = bestConstant.value();
			Float quarterIntervalWidth = originalInterval.width()/4;
			Interval halvedInterval = Interval(originalInterval.midpoint()-quarterIntervalWidth,
											   originalInterval.midpoint()+quarterIntervalWidth);
			system.substitute(bestConstant,Real(halvedInterval));
			result[bestConstant]++;
			ratio = getMaxDerivativeWidthRatio(system, mid_der_widths, bounding_domain);
		}

		system.substitute(working_constants);
	}

	return result;
}

RealConstant
getBestConstantToSplit(
		SystemType& system,
		const RealConstantSet& working_constants,
		const HybridFloatVector& referenceWidths,
		const HybridBoxes& bounding_domain)
{
	RealConstant bestConstant = *working_constants.begin();
	Float bestLocalRatio = std::numeric_limits<Float>::infinity();

	for (RealConstantSet::const_iterator constant_it = working_constants.begin();
												 constant_it != working_constants.end();
												 ++constant_it) {
		// Modifies the system in order to have the range of the original given constant halved
		Real originalValue = system.accessible_constant_value(constant_it->name());
		Float quarterIntervalWidth = originalValue.width()/4;
		Interval halvedInterval = Interval(constant_it->value().midpoint()-quarterIntervalWidth,
										   constant_it->value().midpoint()+quarterIntervalWidth);
		system.substitute(*constant_it,Real(halvedInterval));

		Float localRatio = getMaxDerivativeWidthRatio(system,referenceWidths,bounding_domain);
		if (localRatio < bestLocalRatio) {
			bestLocalRatio = localRatio;
			bestConstant = RealConstant(constant_it->name(),originalValue);
		}

		// Restores the related constant to its original value
		system.substitute(*constant_it,originalValue);
	}

	return bestConstant;
}

void
pushSplitTargetEnclosures(
		std::list<EnclosureType>& initial_enclosures,
		const DiscreteState& target_loc,
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
getDerivativeWidths(
		const HybridAutomaton& system,
		const HybridBoxes& bounding_domain)
{
	ARIADNE_ASSERT_MSG(bounding_domain.size() == system.state_space().size(), "The bounding domain must be defined");

	HybridFloatVector result;

	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = system.state_space().locations_begin()->second;

	for (list<DiscreteMode>::const_iterator modes_it = system.modes().begin(); modes_it != system.modes().end(); modes_it++) {
		const DiscreteState& loc = modes_it->location();

		// Gets the first order derivatives in respect to the dynamic of the mode, applied to the domain of the corresponding location
		Vector<Interval> der = modes_it->dynamic()(bounding_domain.find(loc)->second);

		Vector<Float> der_widths(css);
		for (uint i=0;i<css;i++)
			der_widths[i] = der[i].width();

		result.insert(pair<DiscreteState,Vector<Float> >(loc,der_widths));
	}

	return result;
}

Float
getMaxDerivativeWidthRatio(
		const HybridAutomaton& system,
		const HybridFloatVector& referenceWidths,
		const HybridBoxes& bounding_domain)
{
	ARIADNE_ASSERT_MSG(bounding_domain.size() == system.state_space().size(), "The bounding domain must be defined.");

	Float result = 0;

	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = system.state_space().locations_begin()->second;

	// For each location and dimension, updates the result with the (w - wm)/wm derivative width ratio, excluding the undefined wm = 0 case
	for (list<DiscreteMode>::const_iterator modes_it = system.modes().begin(); modes_it != system.modes().end(); modes_it++) {
		const DiscreteState& loc = modes_it->location();

		Vector<Interval> der = modes_it->dynamic()(bounding_domain.find(loc)->second);

		for (uint i=0; i<css; ++i) {
			Float referenceWidth = referenceWidths.find(loc)->second[i];
			if (referenceWidth != 0)
				result = max(result,(der[i].width()-referenceWidth)/referenceWidth);
		}

	}

	return result;
}

/*! \brief Splits a RealConstant \a con into \a numParts parts.
 * \details Orders the subintervals by putting the second leftmost subinterval up to the rightmost, followed by the leftmost. */
std::vector<RealConstant>
split(
		const RealConstant& con,
		uint numParts)
{
	Interval bounds;
	Float lower, upper;

	std::vector<RealConstant> result(numParts,con);

	String name = con.name();
	Float intervalWidth = con.value().width();

	// Puts the first element
	lower = con.value().lower();
	upper = con.value().lower() + intervalWidth/numParts;
	result[numParts-1] = RealConstant(name,Interval(lower,upper));
	// Puts the last to the second element, in inverse order
	for (uint i=numParts;i>1;--i) {
		lower = con.value().lower() + intervalWidth*(i-1)/numParts;
		upper = con.value().lower() + intervalWidth*i/numParts;
		result[i-2] = RealConstant(name,Interval(lower,upper));
	}

	return result;
}

void fillSplitSet(
		const std::vector<std::vector<RealConstant> >& src,
		std::vector<std::vector<RealConstant> >::iterator col_it,
		std::vector<RealConstant>::iterator row_it,
		RealConstantSet s,
		std::list<RealConstantSet>& dest)
{
	if (col_it != src.end() && row_it != col_it->end()) {
		fillSplitSet(src,col_it,row_it+1,s,dest);
		s.insert(*row_it);
		col_it++;
		row_it = col_it->begin();

		if (col_it != src.end())
			fillSplitSet(src,col_it,row_it,s,dest);
		else
			dest.push_back(s);
	}
}

std::list<RealConstantSet>
getSplitConstantsMidpointsSet(const std::list<RealConstantSet>& intervals_set)
{
	std::list<RealConstantSet> result;

	for (std::list<RealConstantSet>::const_iterator set_it = intervals_set.begin(); set_it != intervals_set.end(); ++set_it) {
		RealConstantSet midpoints;
		for (RealConstantSet::const_iterator constant_it = set_it->begin(); constant_it != set_it->end(); ++constant_it) {
			midpoints.insert(RealConstant(constant_it->name(),constant_it->value().midpoint()));
		}
		result.push_back(midpoints);
	}

	return result;
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
	std::map<DiscreteState,Vector<Float> > hybridgridlengths;

	// Get the minimum domain length for each variable
	Vector<Float> minDomainLengths(css);
	for (uint i=0;i<css;i++) {
		minDomainLengths[i] = std::numeric_limits<double>::infinity();
	}
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		for (uint i=0;i<css;i++) {
			minDomainLengths[i] = min(minDomainLengths[i], domain.find(hfv_it->first)->second[i].width());
		}
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
		hybridgridlengths.insert(make_pair<DiscreteState,Vector<Float> >(hfv_it->first,gridlengths));
	}

	// Populate the grid, centered on the centre of the domain
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++) {
		const DiscreteState& loc = hfv_it->first;
		hg[loc] = Grid(domain.find(loc)->second.centre(),hybridgridlengths[loc]);
	}

	return hg;
}


Vector<Float>
getMaximumEnclosureCell(
		const HybridGrid& hgrid,
		int maximum_grid_depth)
{
	// Introduces a ratio (>1) in respect to the grid cell
	// NOTE: it is preferable to have the ratio slightly lesser than an integer multiple of the grid cell, so that
	// the overapproximation error due to discretization is minimized
	static const double RATIO = 1.9;

	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hgrid.locations_begin()->second.lengths().size();

	// Initializes the result
	Vector<Float> result(css);

	// For each location and dimension of the space
	for (HybridGrid::locations_const_iterator hg_it = hgrid.locations_begin(); hg_it != hgrid.locations_end(); hg_it++)
		for (uint i=0;i<css;i++)
			if (hg_it->second.lengths()[i] > result[i])
				result[i] = hg_it->second.lengths()[i];

	// Scales the cell in respect to the maximum grid depth
	for (uint i=0;i<css;i++)
		result[i] /= (1<<maximum_grid_depth);

	return RATIO*result;
}

Float
getLockToGridTime(
		const HybridAutomaton& system,
		const HybridBoxes& domain)
{
	// Get the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = system.state_space().locations_begin()->second;

	// The variable for the bounding box of the derivatives
	Vector<Interval> der;
    // The variable for the result
	Float result = 0;

	// For each mode
	for (list<DiscreteMode>::const_iterator modes_it = system.modes().begin(); modes_it != system.modes().end(); modes_it++)
	{
		// Gets the location
		const DiscreteState& loc = modes_it->location();

		// Gets the domain for this mode
		const Box& loc_domain = domain.find(loc)->second;

		// Gets the first order derivatives in respect to the dynamic of the mode, applied to the domain of the corresponding location
		der = modes_it->dynamic()(loc_domain);

		// Updates the lock time
		for (uint i=0;i<css;i++)
		{
			Float maxAbsDer = abs(der[i]).upper();
			if (maxAbsDer > 0)
				result = max(result,loc_domain[i].width()/maxAbsDer);
		}
	}

	return result;
}

HybridFloatVector
getHybridMaximumAbsoluteDerivatives(
		const HybridAutomaton& system,
		const HybridGridTreeSet& outer_approx_constraint,
		const HybridBoxes& domain_constraint)
{
	HybridFloatVector result;

	// Get the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = system.state_space().locations_begin()->second;

	// The variable for the bounding box of the derivatives
	Vector<Interval> der;

	// For each mode
	for (list<DiscreteMode>::const_iterator modes_it = system.modes().begin(); modes_it != system.modes().end(); modes_it++) {

		const DiscreteState& loc = modes_it->location();

		// Insert the corresponding hmad pair, initialized with zero maximum absolute derivatives
		result.insert(pair<DiscreteState,Vector<Float> >(loc,Vector<Float>(css)));

		// If the reached region for the location exists and is not empty, check its cells, otherwise use the whole domain
		if (outer_approx_constraint.has_location(loc) && !outer_approx_constraint[loc].empty()) {
			// Get the GridTreeSet
			GridTreeSet reach = outer_approx_constraint[loc];
			// For each of its hybrid cells
			for (GridTreeSet::const_iterator cells_it = reach.begin(); cells_it != reach.end(); cells_it++) {
				// Gets the derivative bounds
				der = modes_it->dynamic()(cells_it->box());

				// For each variable, sets the maximum value
				for (uint i=0;i<css;i++)
					result[loc][i] = max(result[loc][i], abs(der[i]).upper());
			}
		} else {
			// Gets the first order derivatives in respect to the dynamic of the mode, applied to the domain of the corresponding location
			der = modes_it->dynamic()(domain_constraint.find(loc)->second);

			// Gets the maximum absolute derivatives
			for (uint i=0;i<css;i++)
				result[loc][i] = abs(der[i]).upper();
		}
	}

	// Returns
	return result;
}


std::map<DiscreteState,Float>
getHybridMaximumStepSize(
		const HybridFloatVector& hmad,
		const HybridGrid& hgrid,
		int maximum_grid_depth,
		Semantics semantics)
{
	// We choose a coefficient for upper semantics such that an enclosure at maximum size is able to cross
	// urgent transitions in one step. For lower semantics we prefer to have a finer result.
	Float coefficient = (semantics == UPPER_SEMANTICS ? 2.0 : 1.0);

	// Gets the size of the continuous space (NOTE: taken as equal for all locations)
	const uint css = hmad.begin()->second.size();

	// Initialize the hybrid maximum step size
	std::map<DiscreteState,Float> hmss;

	// For each couple DiscreteState,Vector<Float>
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
	{
		// Initializes the maximum step size
		Float mss = 0.0;
		// For each dimension of the space, if the derivative is not zero,
		// evaluates the ratio between the minimum cell length and the derivative itself
		for (uint i=0;i<css;i++)
			if (hfv_it->second[i] > 0)
				mss = max(mss,hgrid[hfv_it->first].lengths()[i]/(1 << maximum_grid_depth)/hfv_it->second[i]);

		// Inserts the value (twice the value since the maximum enclosure is set as ~2 the grid cell)
		hmss.insert(std::pair<DiscreteState,Float>(hfv_it->first,coefficient*mss));
	}

	return hmss;
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


} // namespace Ariadne
