/*****************************************************************************************
 *            workers.h
 *
 *  Copyright  2010  Luca Geretti
 *
 *****************************************************************************************/

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

#include <boost/thread.hpp>
#include <boost/bind.hpp>

#include "ariadne.h"

/* Provides multi-threaded workers for various functions */

// Worker for the _upper_reach_evolve routine in reachability_analyser.cc
class UpperReachEvolveWorker
{
public:

	typedef HybridDenotableSet HDS;
    typedef HybridEvolver::EnclosureType EnclosureType;
	typedef HybridEvolver::ContinuousEnclosureType CE;
	typedef std::list<EnclosureType> EL;
	typedef ListSet<EnclosureType> ELS;
	typedef boost::shared_ptr<EvolverInterface<HybridAutomatonInterface,EnclosureType> > EvolverType;

	// Constructor
    UpperReachEvolveWorker(
    		const EvolverType& evolver,
    		list<EnclosureType>& initial_enclosures,
    		const HybridTime& time,
    		const HybridGrid& grid,
    		const int& accuracy,
    		bool ignore_activations,
    		ContinuousEvolutionDirection continuous_direction,
    		const uint& concurrency)
	: _evolver(evolver),
	  _initial_enclosures(initial_enclosures),
	  _time(time),
	  _grid(grid),
	  _accuracy(accuracy),
	  _ignore_activations(ignore_activations),
	  _continuous_direction(continuous_direction),
	  _concurrency(concurrency)
    {
		_reach = HDS(grid);
		_evolve = HDS(grid);
    }
 
    ~UpperReachEvolveWorker()
    {
    }

    std::pair<HDS,HDS> get_result() 
    {
    	_start();
    	_wait_completion();

		return make_pair<HDS,HDS>(_reach,_evolve);
    }
 
private:

	// A reference to the input variables
	const EvolverType& _evolver;
	list<EnclosureType>& _initial_enclosures;
	const HybridTime& _time;
	const HybridGrid& _grid;
	const int& _accuracy;
	const bool& _ignore_activations;
	const ContinuousEvolutionDirection& _continuous_direction;
	const uint& _concurrency;

	HDS _reach, _evolve;

    std::list<boost::shared_ptr<boost::thread> > _m_threads;
    boost::mutex _inp_mutex, _out_mutex;
 
	void _start()
	{
		for (uint i=0;i<_concurrency;i++)
	        _m_threads.push_back(boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&UpperReachEvolveWorker::_compute, this))));
	}

    void _compute() 
    {
        while (true) {
			_inp_mutex.lock();

			if (_initial_enclosures.empty()) {
				_inp_mutex.unlock();					
				break;
			} else {

				EnclosureType enclosure = _initial_enclosures.front();
				_initial_enclosures.pop_front();
				_inp_mutex.unlock();		

				// Get the enclosures from the initial enclosure, in a lock_time flight
				ELS current_reach_enclosures, current_evolve_enclosures;
				make_ltuple<ELS,ELS>(current_reach_enclosures,current_evolve_enclosures) =
										_evolver->reach_evolve(enclosure,_time,_ignore_activations,_continuous_direction,UPPER_SEMANTICS);

				_out_mutex.lock();

				// Get the discretisation
				HDS current_reach = outer_approximation(current_reach_enclosures,_grid,_accuracy);
				HDS current_evolve = outer_approximation(current_evolve_enclosures,_grid,_accuracy);

		        _reach.adjoin(current_reach);
        		_evolve.adjoin(current_evolve);

				_out_mutex.unlock();
			}
        }
    }       

    void _wait_completion() {
		for (std::list<boost::shared_ptr<boost::thread> >::iterator it = _m_threads.begin(); it != _m_threads.end(); it++)
			(*it)->join();
    }             
};

// Worker for the lower_reach_and_epsilon routine
class LowerReachEpsilonWorker
{
public:

	typedef HybridDenotableSet HDS;
    typedef HybridEvolver::EnclosureType EnclosureType;
	typedef HybridEvolver::ContinuousEnclosureType CE;
	typedef std::list<EnclosureType> EL;
	typedef ListSet<EnclosureType> ELS;
	typedef std::map<DiscreteLocation,uint> HUM;
	typedef boost::shared_ptr<EvolverInterface<HybridAutomatonInterface,EnclosureType> > EvolverType;

	// Constructor
    LowerReachEpsilonWorker(
    		const EvolverType& evolver,
			EL& initial_enclosures,
			const HybridTime& time,
			const HybridGrid& grid,
			const int& accuracy,
			const uint& concurrency)
	: _evolver(evolver),
	  _initial_enclosures(initial_enclosures),
	  _time(time),
	  _grid(grid),
	  _accuracy(accuracy),
	  _concurrency(concurrency)
    {
    	_reach = HDS(grid);
		_evolve_global = HDS(grid);

    	HybridSpace state_space = _evolver->system().state_space();
    	for (HybridSpace::const_iterator hs_it = state_space.begin(); hs_it != state_space.end(); ++hs_it)
    		_epsilon.insert(std::pair<DiscreteLocation,Vector<Float> >(hs_it->first,Vector<Float>(hs_it->second)));
    }
 
    ~LowerReachEpsilonWorker()
    {
    }

    // Create the threads and produce the required results
    tuple<std::pair<HUM,HUM>,EL,HDS,HybridFloatVector> get_result()
    {
    	_start();
    	_wait_completion();

		// Calculate the adjoined evolve sizes
		for (HDS::locations_const_iterator evolve_global_it = _evolve_global.locations_begin(); evolve_global_it != _evolve_global.locations_end(); evolve_global_it++)
			_adjoined_evolve_sizes[evolve_global_it->first] = evolve_global_it->second.size();

		return make_tuple<std::pair<HUM,HUM>,EL,HDS,HybridFloatVector>(make_pair<HUM,HUM>(_adjoined_evolve_sizes,_superposed_evolve_sizes),_final_enclosures,_reach,_epsilon);
    }
 
private:

	// A reference to the input variables
	const EvolverType& _evolver;
	EL& _initial_enclosures;
	const HybridTime& _time;
	const HybridGrid& _grid;
	const int& _accuracy;
	const uint& _concurrency;

	HDS _evolve_global;
	HUM _adjoined_evolve_sizes, _superposed_evolve_sizes;
	EL _final_enclosures;
	HDS _reach;
	HybridFloatVector _epsilon;

    std::list<boost::shared_ptr<boost::thread> > _m_threads;
    boost::mutex _inp_mutex, _out_mutex;
 
	void _start()
	{
		for (uint i=0;i<_concurrency;i++)
	        _m_threads.push_back(boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&LowerReachEpsilonWorker::_compute, this))));
	}

    void _compute() 
    {
        const Float accuracy_divider = 1 << _accuracy;

        while (true) {
			_inp_mutex.lock();

			if (_initial_enclosures.empty()) {
				_inp_mutex.unlock();					
				break;
			} else {
				HybridBasicSet<CE> current_initial_enclosure = _initial_enclosures.front();
				_initial_enclosures.pop_front();
				_inp_mutex.unlock();

				// Get the enclosures from the initial enclosure, in a lock_time flight
				ELS current_reach_enclosures, current_evolve_enclosures;
				make_ltuple<ELS,ELS>(current_reach_enclosures,current_evolve_enclosures) =
										_evolver->reach_evolve(current_initial_enclosure,_time,LOWER_SEMANTICS);

				_out_mutex.lock();

				// Get the discretisation
				HDS current_reach = outer_approximation(current_reach_enclosures,_grid,_accuracy);
				HDS current_evolve = outer_approximation(current_evolve_enclosures,_grid,_accuracy);

				_reach.adjoin(current_reach);
				_evolve_global.adjoin(current_evolve);

			    // Update the epsilon
				for (ELS::const_iterator encl_it = current_reach_enclosures.begin(); encl_it != current_reach_enclosures.end(); ++encl_it) {
					Vector<Float> encl_widths = encl_it->second.bounding_box().widths();
					Vector<Float> minimum_cell_lengths = _grid[encl_it->first].lengths()/accuracy_divider;
					_epsilon[encl_it->first] = max_elementwise(_epsilon[encl_it->first],minimum_cell_lengths);
				}

				// Add the number of cells of the current evolve to the superposed total at the end of the step, for each location
				for (HDS::locations_const_iterator evolve_it = current_evolve.locations_begin(); evolve_it != current_evolve.locations_end(); evolve_it++)
					_superposed_evolve_sizes[evolve_it->first] += current_evolve[evolve_it->first].size();
				// Add the current_enclosures to the final enclosures
				for (ELS::locations_const_iterator loc_it = current_evolve_enclosures.locations_begin(); loc_it != current_evolve_enclosures.locations_end(); loc_it++)
					for (ListSet<CE>::const_iterator list_it = loc_it->second.begin(); list_it != loc_it->second.end(); list_it++)
						_final_enclosures.push_back(EnclosureType(loc_it->first,*list_it));

				_out_mutex.unlock();
			}
        }
    }       

    void _wait_completion()
    {
		for (std::list<boost::shared_ptr<boost::thread> >::iterator it = _m_threads.begin(); it != _m_threads.end(); it++)
			(*it)->join();
    }             
};


