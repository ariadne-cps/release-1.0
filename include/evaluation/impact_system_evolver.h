/***************************************************************************
 *            impact_system_evolver.h
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
 
/*! \file impact_system_evolver.h
 *  \brief Methods for computing a single time step of the evolution of a vector_field.
 */

#ifndef ARIADNE_IMPACT_SYSTEM_EVOLVER_H
#define ARIADNE_IMPACT_SYSTEM_EVOLVER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

#include "base/tuple.h"

#include "evaluation/applicator_interface.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/satisfier_interface.h"
#include "evaluation/subdivider_interface.h"
#include "evaluation/reducer_interface.h"
#include "evaluation/evolver_interface.h"

#include "evaluation/evolver_base.h"

#include "evaluation/evolution_parameters.h"



namespace Ariadne {
  
    template<class R> class ImpactSystem;
    template<class R> class ApproximateTaylorModel;

    template<class Sys, class ES> class Evolver;
 
    /*! \ingroup Evolvers
     *  \brief A class for evolving a continuous-time dynamical system.
     */
    template<class ES> 
    class Evolver< ImpactSystem<typename ES::real_type>, ES>
      : public EvolverBase< ImpactSystem<typename ES::real_type>, ES>
    {
      typedef typename ES::real_type R;
      typedef typename traits<R>::approximate_arithmetic_type A;
      typedef typename traits<R>::interval_type I;
      typedef ImpactSystem<R> Sys;
      typedef typename Sys::time_type T;
      typedef Vector<I> IVec;
      typedef Box<R> Bx;
      typedef ListSet<ES> ESL;
      typedef ApproximateTaylorModel<R> ATM;
      typedef TimedSet<T,ES> TES;
      typedef ListSet<TES> TESL;
      typedef VectorField<R> VF;
      typedef FunctionInterface<R> FN;
     public:
      Evolver();
      Evolver(const EvolutionParameters<R>&);
      Evolver(const EvolutionParameters<R>&,const ApplicatorInterface<ES>&, const IntegratorInterface<ES>&);
      Evolver(const EvolutionParameters<R>&,const ApplicatorInterface<ES>&, const IntegratorInterface<ES>&, const SatisfierInterface<ES>&, const SubdividerInterface<ES>&, const ReducerInterface<ES>&);
      virtual Evolver<Sys,ES>* clone() const { return new Evolver<Sys,ES>(*this); }
     protected:
      virtual void _evolution(ESL& final, ESL& reachable, ESL& intermediate, const Sys& system, const ES& initial, const T& time, Semantics semantics, bool reach) const;
     public:
      using EvolverBase<Sys,ES>::evolve;
      using EvolverBase<Sys,ES>::reach;
      /*! \brief Compute an approximation to the evolution set under upper semantics. */
      ESL evolve(const Sys& system, const ES& initial_set, const T& time) const {
        ESL final; ESL reachable; ESL intermediate; this->_evolution(final,reachable,intermediate,system,initial_set,time,upper_semantics,false); return final; }
      /*! \brief Compute an approximation to the evolution set under upper semantics. */
      ESL reach(const Sys& system, const ES& initial_set, const T& time) const {
        ESL final; ESL reachable; ESL intermediate; this->_evolution(final,reachable,intermediate,system,initial_set,time,upper_semantics,true); return reachable; }
     public:
      void evolution(ESL& final_sets, ESL& reach_sets, ESL& intermediate_sets, const Sys& system, const ES& initial_set, const T& time) const {
        return this->_evolution(final_sets, reach_sets, intermediate_sets, system, initial_set, time, upper_semantics, true); }
     private:
      Interval<R> _grazing_time_interval(const ATM& flow_model, const ATM& guard_model, const ATM& initial_set_model, 
                                         const A& initial_time, const A& final_time) const;
      Interval<R> _crossing_time_interval(const ATM& flow_model, const ATM& guard_model, const ATM& initial_set_model, 
                                          const A& initial_time, const A& final_time) const;
      ATM _hitting_time_model(const ATM& flow_model, const ATM& guard_model, const ATM& initial_set_model, 
                              const A& initial_time, const A& final_time) const;

      ATM _reset_step(const ATM& map_model, const ATM& set_model) const;
      ATM _integration_step(const ATM& flow_model, const ATM& initial_set_model, const A& integration_time) const;
      ATM _integration_step(const ATM& flow_model, const ATM& initial_set_model, const ATM& integration_time_model) const;
      ATM _reachability_step(const ATM& flow_model, const ATM& initial_set_model, const A& initial_time, const A& final_time) const;
      ATM _reachability_step(const ATM& flow_model, const ATM& initial_set_model, const A& initial_time, const ATM& final_time_model) const;
      ATM _reachability_step(const ATM& flow_model, const ATM& initial_set_model, const ATM& initial_time_model, const ATM& final_time_model) const;
      ATM _reachability_time(const ATM& initial_time_model, const A& start_time, const A& finish_time) const;

      tribool _active(const ATM& guard_model, const ATM& _set_model) const;
      tribool _active(const FN& guard_function, const ATM& _set_model) const;
     private:
      // Helper functions for accessing parameters
      Rational maximum_step_size() const { 
        return this->_parameters->maximum_step_size(); }
      R maximum_enclosure_radius() const { 
        return this->_parameters->maximum_enclosure_radius(); }
     private:
      // Services provided by other classes
      std::pair<Rational,IVec> flow_bounds(const FN& vf, const IVec& bx) const {
        std::pair<Rational,Bx> bounds=this->_integrator->flow_bounds(VF(vf),Bx(bx),this->maximum_step_size()); 
        return std::make_pair(bounds.first,bounds.second.position_vectors()); }
      std::pair<Rational,IVec> flow_bounds(const VF& vf, const IVec& bx, const Rational& h) const {
        std::pair<Rational,Bx> bounds=this->_integrator->flow_bounds(VF(vf),Bx(bx),h); 
        return std::make_pair(bounds.first,bounds.second.position_vectors()); }
      R radius(const ES& s) const {
        return s.radius(); }
      IVec bounding_box(const ES& s) const {
        return s.bounding_box().position_vectors(); }
      ESL subdivide(const ES& s) const {
        return this->_subdivider->subdivide(s,this->maximum_enclosure_radius()); }
     private:
      // Helper functions for converting between sets and models.
      ApproximateTaylorModel<R> model(const ES& s, ushort d) const;
      ES set(const ApproximateTaylorModel<R>& s) const;
     private:
      // Helper functions for timed sets
      R radius(const TES& tes) const { return tes.set().radius(); }
      void adjoin_subdivision(TESL& tls, const TES& ts) const { 
        T t=ts.time();
        ESL subdivisions=this->_subdivider->subdivide(ts.set(),this->maximum_enclosure_radius());
        for(size_type i=0; i!=subdivisions.size(); ++i) {
          tls.adjoin(TES(t,subdivisions[i]));
        }
      }
      TES reduce(const TES& ts) const {
        return TES(ts.time(),this->_reducer->over_approximate(ts.set())); }
     private:
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
      boost::shared_ptr< IntegratorInterface<ES> >  _integrator;
      boost::shared_ptr< SubdividerInterface<ES> > _subdivider;
      boost::shared_ptr< ReducerInterface<ES> > _reducer;

      boost::shared_ptr< EvolutionProfiler >  _profiler;
     public:
      mutable uint verbosity;
    };



  
} // namespace Ariadne



#endif /* ARIADNE_IMPACT_SYSTEM_EVOLVER_H */