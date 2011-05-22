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

	typedef HybridGridTreeSet HGTS;
    typedef HybridEvolver::EnclosureType EnclosureType;
	typedef HybridEvolver::ContinuousEnclosureType CE;

	// Constructor
    UpperReachEvolveWorker(const boost::shared_ptr<HybridDiscretiser<CE> >& discretiser,
    					   const HybridAutomaton& sys,
    					   const HGTS& initial_set,
    					   const HybridTime& time,
    					   const HybridGrid& grid,
    					   const int& accuracy,
    					   const uint& concurrency)
	: _discretiser(discretiser),
	  _sys(sys), 
	  _initial_set(initial_set),
	  _time(time),
	  _grid(grid),
	  _accuracy(accuracy),
	  _concurrency(concurrency),
	  _cells_it(_initial_set.begin())
    {
		_reach = HGTS(grid);
		_evolve = HGTS(grid);
		_cells_it = _initial_set.begin();
    }
 
    ~UpperReachEvolveWorker()
    {
    }

    std::pair<HGTS,HGTS> get_result() 
    {
    	_start();
    	_wait_completion();

		return make_pair<HGTS,HGTS>(_reach,_evolve);
    }
 
private:

	// A reference to the input variables
	const boost::shared_ptr<HybridDiscretiser<CE> >& _discretiser;
	const HybridAutomaton& _sys;
	const HGTS& _initial_set;
	const HybridTime& _time;
	const HybridGrid& _grid;
	const int& _accuracy;
	const uint& _concurrency;
	// The sets in which the _compute function adjoins results
	HGTS _reach, _evolve;

	// The various threads
    std::list<boost::shared_ptr<boost::thread> > _m_threads;
	// The mutexes for synchronization
    boost::mutex _inp_mutex, _out_mutex;
	// The iterator for the cells
	HGTS::const_iterator _cells_it;

	friend class HybridDiscretiser<CE>;
 
	// Create and start the threads
	void _start()
	{
		for (uint i=0;i<_concurrency;i++)
	        _m_threads.push_back(boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&UpperReachEvolveWorker::_compute, this))));
	}

    // Compute the result
    void _compute() 
    {
		// Execution loop
        while (true)
        {
			_inp_mutex.lock();

			// If all cells have been picked
			if (_cells_it == _initial_set.end())
			{
				_inp_mutex.unlock();
				break;
			}
            else
			{
		        EnclosureType enclosure=_discretiser->enclosure(*_cells_it);
				++_cells_it;
				_inp_mutex.unlock();		

        		HGTS reach, evolve;
		        make_lpair(reach,evolve)=_discretiser->evolution(_sys,enclosure,_time,_grid,_accuracy,UPPER_SEMANTICS);
				
				_out_mutex.lock();
		        _reach.adjoin(reach);
        		_evolve.adjoin(evolve);
				_out_mutex.unlock();
			}
        }
    }       

	// Wait for completion of all threads
    void _wait_completion()
    {
		for (std::list<boost::shared_ptr<boost::thread> >::iterator it = _m_threads.begin(); it != _m_threads.end(); it++)
			(*it)->join();
    }             
};

// Worker for the _upper_reach_evolve_continuous routine in reachability_analyser.cc
class UpperReachEvolveContinuousWorker
{
public:

	typedef HybridGridTreeSet HGTS;
    typedef HybridEvolver::EnclosureType EnclosureType;
	typedef HybridEvolver::ContinuousEnclosureType CE;

	// Constructor
    UpperReachEvolveContinuousWorker(const boost::shared_ptr<HybridDiscretiser<CE> >& discretiser,
    								 const HybridAutomaton& sys,
    								 const list<EnclosureType>& initial_enclosures,
    								 const HybridTime& time,
    								 EvolutionDirection direction,
    								 const HybridGrid& grid,
    								 const int& accuracy,
    								 const uint& concurrency)
	: _discretiser(discretiser),
	  _sys(sys), 
	  _initial_enclosures(initial_enclosures),
	  _time(time),
	  _direction(direction),
	  _grid(grid),
	  _accuracy(accuracy),
	  _concurrency(concurrency),
	  _enclosures_it(_initial_enclosures.begin())
    {
		_reach = HGTS(grid);
		_evolve = HGTS(grid);
		_enclosures_it = _initial_enclosures.begin();
    }
 
    ~UpperReachEvolveContinuousWorker()
    {
    }

    std::pair<HGTS,HGTS> get_result() 
    {
    	EvolutionDirection saved_direction = _discretiser->evolver()->settings().direction;
    	_discretiser->evolver()->settings().direction = _direction;
    	_discretiser->evolver()->settings().enable_premature_termination_on_blocking_event = true;

    	_start();
    	_wait_completion();

    	_discretiser->evolver()->settings().direction = saved_direction;
    	_discretiser->evolver()->settings().enable_premature_termination_on_blocking_event = false;

		return make_pair<HGTS,HGTS>(_reach,_evolve);
    }
 
private:

	// A reference to the input variables
	const boost::shared_ptr<HybridDiscretiser<CE> >& _discretiser;
	const HybridAutomaton& _sys;
	const list<EnclosureType>& _initial_enclosures;
	const HybridTime& _time;
	const EvolutionDirection& _direction;
	const HybridGrid& _grid;
	const int& _accuracy;
	const uint& _concurrency;
	// The sets in which the _compute function adjoins results
	HGTS _reach, _evolve;

	// The various threads
    std::list<boost::shared_ptr<boost::thread> > _m_threads;
	// The mutexes for synchronization
    boost::mutex _inp_mutex, _out_mutex;
	// The iterator for the enclosures
	list<EnclosureType>::const_iterator _enclosures_it;

	friend class HybridDiscretiser<CE>;
 
	// Create and start the threads
	void _start()
	{
		for (uint i=0;i<_concurrency;i++)
	        _m_threads.push_back(boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&UpperReachEvolveContinuousWorker::_compute, this))));
	}

    // Compute the result
    void _compute() 
    {
		// Execution loop
        while (true)
        {
			_inp_mutex.lock();

			// If all enclosures have been picked
			if (_enclosures_it == _initial_enclosures.end())
			{
				_inp_mutex.unlock();					
				break;
			}
            else
			{
				EnclosureType enclosure = *_enclosures_it;
				++_enclosures_it;
				_inp_mutex.unlock();		

        		HGTS reach, evolve;
		        make_lpair(reach,evolve)=_discretiser->evolution(_sys,enclosure,_time,_grid,_accuracy,UPPER_SEMANTICS);

				_out_mutex.lock();
		        _reach.adjoin(reach);
        		_evolve.adjoin(evolve);
				_out_mutex.unlock();
			}
        }
    }       

	// Wait for completion of all threads
    void _wait_completion()
    {
		for (std::list<boost::shared_ptr<boost::thread> >::iterator it = _m_threads.begin(); it != _m_threads.end(); it++)
			(*it)->join();
    }             
};

// Worker for the lower_chain_reach routine
class LowerChainReachWorker
{
public:

	typedef HybridGridTreeSet HGTS;
    typedef HybridEvolver::EnclosureType EnclosureType;
	typedef HybridEvolver::ContinuousEnclosureType CE;
	typedef std::list<EnclosureType> EL;
	typedef ListSet<EnclosureType> ELS;
	typedef std::map<DiscreteState,uint> HUM;

	friend class HybridDiscretiser<CE>;

	// Constructor
    LowerChainReachWorker(const boost::shared_ptr<HybridDiscretiser<CE> >& discretiser,
					 EL& initial_enclosures,
					 const HybridAutomaton& sys, 
					 const HybridTime& time, 
					 const HybridGrid& grid,
					 const int& accuracy,
					 const uint& concurrency)
	: _discretiser(discretiser),
	  _initial_enclosures(initial_enclosures),
	  _sys(sys), 
	  _time(time),
	  _grid(grid),
	  _accuracy(accuracy),
	  _concurrency(concurrency),
	  _eps_lower_bounds(EpsilonLowerBounds(sys.state_space()))
    {
    	_reach = HGTS(grid);
		_evolve_global = HGTS(grid);
    }
 
    ~LowerChainReachWorker()
    {
    }

    // Create the threads and produce the required results
    tuple<std::pair<HUM,HUM>,EL,HGTS,EpsilonLowerBounds> get_result()
    {
    	_start();
    	_wait_completion();

		// Calculate the adjoined evolve sizes
		for (HGTS::locations_const_iterator evolve_global_it = _evolve_global.locations_begin(); evolve_global_it != _evolve_global.locations_end(); evolve_global_it++)
			_adjoined_evolve_sizes[evolve_global_it->first] = evolve_global_it->second.size();

		return make_tuple<std::pair<HUM,HUM>,EL,HGTS,EpsilonLowerBounds>(make_pair<HUM,HUM>(_adjoined_evolve_sizes,_superposed_evolve_sizes),_final_enclosures,_reach,_eps_lower_bounds);
    }
 
private:

	// A reference to the input variables
    const boost::shared_ptr<HybridDiscretiser<CE> >& _discretiser;
	EL& _initial_enclosures;
	const HybridAutomaton& _sys;
	const HybridTime& _time;
	const HybridGrid& _grid;
	const int& _accuracy;
	const uint& _concurrency;
	// The adjoined evolve common to all threads
	HGTS _evolve_global;
	// The sizes of the adjoined or superposed evolve
	HUM _adjoined_evolve_sizes, _superposed_evolve_sizes;
	// The final enclosures
	EL _final_enclosures;
	// The reached region
	HGTS _reach;
	// The epsilon bounds
	EpsilonLowerBounds _eps_lower_bounds;

	// The various threads
    std::list<boost::shared_ptr<boost::thread> > _m_threads;
	// The mutexes for synchronization
    boost::mutex _inp_mutex, _out_mutex;
 
	// Create and start the threads
	void _start()
	{
		for (uint i=0;i<_concurrency;i++)
	        _m_threads.push_back(boost::shared_ptr<boost::thread>(new boost::thread(boost::bind(&LowerChainReachWorker::_compute, this))));
	}

    // Compute the result
    void _compute() 
    {
		// Execution loop
        while (true)
        {
			_inp_mutex.lock();

			// If all initial enclosures have been picked
			if (_initial_enclosures.empty())
			{
				_inp_mutex.unlock();					
				break;
			}
            else
			{
				HybridBasicSet<CE> current_initial_enclosure = _initial_enclosures.front();
				_initial_enclosures.pop_front();
				_inp_mutex.unlock();

				HGTS current_reach, current_evolve;
				ELS current_reach_enclosures, current_evolve_enclosures;

				// Get the enclosures from the initial enclosure, in a lock_time flight
				make_ltuple<ELS,ELS>(current_reach_enclosures,current_evolve_enclosures) =
										_discretiser->evolver()->reach_evolve(_sys,current_initial_enclosure,_time,LOWER_SEMANTICS);

				// Get the discretisation
				current_reach = outer_approximation(current_reach_enclosures,_grid,_accuracy);
				current_evolve = outer_approximation(current_evolve_enclosures,_grid,_accuracy);

				_out_mutex.lock();

				_reach.adjoin(current_reach);
				_evolve_global.adjoin(current_evolve);

			    // Update the epsilon lower bounds
				for (ELS::const_iterator encl_it = current_reach_enclosures.begin(); encl_it != current_reach_enclosures.end(); ++encl_it) {
					Box encl_bounds = encl_it->second.bounding_box();
					_eps_lower_bounds.updateEpsilon(encl_it->first,encl_bounds.widths());
					_eps_lower_bounds.updateReachBounds(encl_it->first,encl_bounds);
				}

				// Add the number of cells of the current evolve to the superposed total at the end of the step, for each location
				for (HGTS::locations_const_iterator evolve_it = current_evolve.locations_begin(); evolve_it != current_evolve.locations_end(); evolve_it++)
					_superposed_evolve_sizes[evolve_it->first] += current_evolve[evolve_it->first].size();
				// Add the current_enclosures to the final enclosures
				for (ELS::locations_const_iterator loc_it = current_evolve_enclosures.locations_begin(); loc_it != current_evolve_enclosures.locations_end(); loc_it++)
					for (ListSet<CE>::const_iterator list_it = loc_it->second.begin(); list_it != loc_it->second.end(); list_it++)
						_final_enclosures.push_back(EnclosureType(loc_it->first,*list_it));

				_out_mutex.unlock();
			}
        }
    }       

	// Wait for completion of all threads
    void _wait_completion()
    {
		for (std::list<boost::shared_ptr<boost::thread> >::iterator it = _m_threads.begin(); it != _m_threads.end(); it++)
			(*it)->join();
    }             
};


