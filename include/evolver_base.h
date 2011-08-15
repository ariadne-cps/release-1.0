/***************************************************************************
 *            evolver_base.h
 *
 *  Copyright  2008  Pieter Collins
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
 
/*! \file evolver_base.h
 *  \brief Interface for computing a single time step of the evolution of a system.
 */

#ifndef ARIADNE_EVOLVER_BASE_H
#define ARIADNE_EVOLVER_BASE_H

#include "evolver_interface.h"
#include "list_set.h"

namespace Ariadne {
  
  
template<class SYS, class ES> class EvolverBase
    : public EvolverInterface<SYS,ES>
{
    typedef EvolverInterface<SYS,ES> Interface;
    typedef typename SYS::TimeType T;
    typedef ListSet<ES> ESL;
    typedef typename ESL::const_iterator ESLCI;

  public:

    virtual std::ostream& write(std::ostream& os) const {
        return os << "Evolver( ... )"; }

    EvolverBase(const SYS& system) :
    	_sys(system.clone())
	{
	}

    virtual const SYS& system() const { return *_sys; }

  public:
    Orbit<ES> orbit(const ES& initial_set, const T& time, Semantics semantics) const {
        Orbit<ES> orbit(initial_set);
        ESL final; ESL reachable; ESL intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,time,false,DIRECTION_FORWARD,semantics);
        orbit.adjoin_intermediate(intermediate); orbit.adjoin_reach(reachable); orbit.adjoin_final(final);
        return orbit;
    }
    //! \brief Compute an approximation to the evolution set under the given semantics. 
    ESL evolve(const ES& initial_set, const T& time, Semantics semantics) const {
        ESL final; ESL reachable; ESL intermediate; this->_evolution(final,reachable,intermediate,initial_set,time,false,DIRECTION_FORWARD,semantics); return final; }
    //! \brief Compute an approximation to the evolution set under the given semantics. 
    ESL reach(const ES& initial_set, const T& time, Semantics semantics) const {
        ESL final; ESL reachable; ESL intermediate; this->_evolution(final,reachable,intermediate,initial_set,time,false,DIRECTION_FORWARD,semantics); return reachable; }
    //! \brief Compute an approximation to the evolution set under the given semantics. 
    std::pair<ESL,ESL> reach_evolve(const ES& initial_set, const T& time, Semantics semantics) const {
        ESL final; ESL reachable; ESL intermediate; this->_evolution(final,reachable,intermediate,initial_set,time,false,DIRECTION_FORWARD,semantics); return std::make_pair(reachable,final); }
    //! \brief Compute an approximation to the evolution set under the given semantics.
    std::pair<ESL,ESL> reach_evolve(const ES& initial_set, const T& time, bool ignore_activations, ContinuousEvolutionDirection continuous_direction, Semantics semantics) const {
    	ESL final; ESL reachable; ESL intermediate; this->_evolution(final,reachable,intermediate,initial_set,time,ignore_activations,continuous_direction,semantics); return std::make_pair(reachable,final); }

  protected:
    virtual void _evolution(ESL& final, ESL& reachable, ESL& intermediate, const ES& initial, const T& time, bool ignore_activations, ContinuousEvolutionDirection direction, Semantics semantics) const = 0;
  protected:
    boost::shared_ptr<SYS> _sys;
};

inline
VectorFunction
get_directed_dynamic(const VectorFunction& dynamic, ContinuousEvolutionDirection direction)
{
	const ScalarFunction minus_one = ScalarFunction::constant(dynamic.result_size(),-1);
	return (direction == DIRECTION_FORWARD ? dynamic : dynamic*minus_one);
}
  
} // namespace Ariadne



#endif // ARIADNE_EVOLVER_BASE_H
