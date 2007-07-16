/***************************************************************************
 *            constraint_hybrid_evolver.code.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
#include <vector>

#include "../geometry/rectangle_expression.h"
#include "../geometry/set_interface.h"
#include "../geometry/hybrid_set.h"
#include "../geometry/timed_set.h"
#include "../system/map.h"
#include "../system/vector_field.h"
#include "../system/constraint_hybrid_automaton.h"
#include "../evaluation/applicator.h"
#include "../evaluation/integrator.h"

#include "../evaluation/lohner_integrator.h"

#include "../output/epsfstream.h"
#include "../output/logging.h"

#include "constraint_hybrid_evolver.h"


namespace {

using namespace Ariadne;

template<class V, class Iter> inline
void append(V& v, Iter begin, Iter end) 
{
  for(;begin!=end;++begin) {
    v.push_back(*begin);
  }
}

template<class V, class C> inline
void append(V& v, const C& c) 
{
  for(typename C::const_iterator iter=c.begin(); iter!=c.end(); ++iter) {
    v.push_back(*iter);
  }
}


} 



namespace Ariadne {
 
namespace Evaluation { static int& verbosity = hybrid_evolver_verbosity; }


template<class R>
Evaluation::ConstraintHybridEvolver<R>::~ConstraintHybridEvolver()
{
}


template<class R>
Evaluation::ConstraintHybridEvolver<R>::ConstraintHybridEvolver(Applicator<R>& a, Integrator<R>& i)
  : _applicator(&a), _integrator(dynamic_cast<LohnerIntegrator<R>*>(&i))
{
  if(!this->_integrator) {
    throw std::runtime_error("ConstraintHybridEvolver::ConstraintHybridEvolver(Applicator a, Integrator i): Invalid integrator");
  }
}



template<class R> inline
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type>
Evaluation::ConstraintHybridEvolver<R>::subdivide(const timed_set_type& timed_set)
{
  std::vector<timed_set_type> result;
  Geometry::ListSet<continuous_basic_set_type> subdivisions=timed_set.continuous_state_set().subdivide();
  for(typename std::vector<continuous_basic_set_type>::const_iterator iter=subdivisions.begin();
      iter!=subdivisions.end(); ++iter)
  {
    result.push_back(timed_set_type(timed_set.time(),timed_set.steps(),timed_set.discrete_state(),*iter));
  }
  return result;
}


template<class R> inline
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type>
Evaluation::ConstraintHybridEvolver<R>::working_sets(const Geometry::HybridListSet<continuous_basic_set_type>& initial_set)
{
  std::vector<timed_set_type> working_sets;
  for(typename Geometry::HybridListSet<continuous_basic_set_type>::const_iterator loc_iter=initial_set.begin();
      loc_iter!=initial_set.end(); ++loc_iter)
  {
    id_type loc_id=loc_iter->first;
    const Geometry::ListSet<continuous_basic_set_type>& list_set=*loc_iter->second;
    for(typename Geometry::ListSet<continuous_basic_set_type>::const_iterator bs_iter=list_set.begin();
        bs_iter!=list_set.end(); ++bs_iter)
    {
      const continuous_basic_set_type& bs=*bs_iter;
      working_sets.push_back(timed_set_type(0,0,loc_id,bs));
    }
  }
  return working_sets;
}





/*! \brief Compute the possible states reached by an unfored jump within time \a h.
 *
 * The generators are given by
 * \f[ D\Phi_2 \circ DF \circ D\Phi_1 G; \quad (D\Phi_i\circ DF \circ \dot{\Phi}_1 - \dot{\Phi}_2) (h/2) \f]
 */
template<class R>
typename Evaluation::ConstraintHybridEvolver<R>::continuous_basic_set_type
Evaluation::ConstraintHybridEvolver<R>::unforced_jump(const System::VectorFieldInterface<R>& dynamic1,
                                                      const System::VectorFieldInterface<R>& dynamic2,
                                                      const System::MapInterface<R>& reset,
                                                      const continuous_basic_set_type& initial_set,
                                                      const Geometry::ConstraintInterface<R>& activation,
                                                      const time_type& step_size)
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;;

  time_type required_step_size=step_size;

  Rectangle<R> bounding_box = initial_set.bounding_box();
  bounding_box = this->_integrator->estimate_flow_bounds(dynamic1,bounding_box,required_step_size);
  Point<I> mode1_bounds = bounding_box;
  Rectangle<R> mode1_bounding_box=bounding_box;
  bounding_box = this->_applicator->evaluate(reset,bounding_box);
  bounding_box = this->_integrator->estimate_flow_bounds(dynamic2,bounding_box,required_step_size);
  Point<I> mode2_bounds = bounding_box;
  Rectangle<R> mode2_bounding_box=bounding_box;
  assert(required_step_size==step_size);
  
  // Compute evolution assuming transition at time step_size/2
  time_type half_step_size=step_size/2;
  continuous_basic_set_type integral_set = initial_set;
  integral_set = this->_integrator->bounded_integration_step(dynamic1,integral_set,mode1_bounding_box,half_step_size);
  integral_set = this->_applicator->evaluate(reset,integral_set);
  integral_set = this->_integrator->bounded_integration_step(dynamic2,integral_set,mode2_bounding_box,half_step_size);
  
  // Compute extra generator associated to transition
  // FIXME: I don't think formula for jacobian is correct
  I hh=half_step_size;
  Vector<I> v=(dynamic2.jacobian(mode2_bounds)*(reset.jacobian(mode1_bounds)*dynamic1.image(mode1_bounds))-dynamic2.image(mode2_bounds))*hh;

  return continuous_basic_set_type(integral_set.centre(),integral_set.generators(),v);
}


/*! \brief Compute the possible states reached by an unfored jump within time \a h.
 *
 * The generators are given by
 * \f[ D\Phi_2 \circ DF \circ D\Phi_1 G; \quad (D\Phi_i\circ DF \circ \dot{\Phi}_1 - \dot{\Phi}_2) (h/2) \f]
 */
template<class R>
typename Evaluation::ConstraintHybridEvolver<R>::continuous_basic_set_type
Evaluation::ConstraintHybridEvolver<R>::forced_jump(const System::VectorFieldInterface<R>& dynamic1,
                                                    const System::VectorFieldInterface<R>& dynamic2,
                                                    const System::MapInterface<R>& reset,
                                                    const continuous_basic_set_type& initial_set,
                                                    const Geometry::ConstraintInterface<R>& guard,
                                                    const time_type& required_time_step)
{
  using namespace Numeric;
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;;

  ARIADNE_LOG(7,"forced_jump: initial_set="<<initial_set<<"\n");
  
  time_type step_size=required_time_step;

  Rectangle<R> mode1_bounding_box = this->_integrator->estimate_flow_bounds(dynamic1,initial_set.bounding_box(),step_size);

  R approximate_normal_flow_speed = inner_product( guard.gradient(mode1_bounding_box.centre()), dynamic1(mode1_bounding_box.centre()) ).centre();
  assert(approximate_normal_flow_speed!=0);

  time_type approximate_hitting_time = Interval<R>(-guard.value(initial_set.centre())/approximate_normal_flow_speed).centre();

    //ARIADNE_LOG(7,"approximate_hitting_time="<<approximate_hitting_time<<"\n");
  assert(approximate_hitting_time>=0);

  Rectangle<R> bounding_box = initial_set.bounding_box();
  bounding_box = this->_integrator->estimate_flow_bounds(dynamic1,bounding_box,step_size);
  Point<I> mode1_bounds = bounding_box;
  bounding_box = this->_applicator->evaluate(reset,bounding_box);
  bounding_box = this->_integrator->estimate_flow_bounds(dynamic2,bounding_box,step_size);
  Point<I> mode2_bounds = bounding_box;
  assert(required_time_step==step_size);
  
  // Compute time to jump as an interval affine map on the state

  // Reference centre and generators of zonotope
  const Point<I>& c=initial_set.centre();
  const Matrix<I>& G=initial_set.generators();
  
  Interval<R> cval=guard.value(c);
  Vector<I> cgrad=guard.gradient(mode1_bounds);
  
  // Direction of the vector field
  Vector<I> dir=dynamic1(mode1_bounds);

  // Normal direction of the vector field
  Interval<R> ndir=LinearAlgebra::inner_product(cgrad,dir);
  if(ndir.contains(0)) {
    throw std::runtime_error("Vector field is not normal to guard set");
  }

  // Time to jump for centre of zonoope
  Interval<R> tcentre=(cval-LinearAlgebra::inner_product(cgrad,c.position_vector()))/ndir;
  // Modification to jump time for generators of zonotope
  Vector<I> tgrad=(cgrad*G)/(-ndir);
  
  continuous_basic_set_type integral_set = initial_set;
  integral_set = this->_integrator->integrate(dynamic1,integral_set,approximate_hitting_time);
  integral_set = this->_applicator->evaluate(reset,integral_set);
  integral_set = this->_integrator->integrate(dynamic2,integral_set,required_time_step-approximate_hitting_time);
  
  // Compute extra generator associated to transition
  // FIXME: I don't think formula for jacobian is correct
  Vector<I> v(dynamic2.jacobian(mode2_bounds)*(reset.jacobian(mode1_bounds)*dynamic1.image(mode1_bounds))-dynamic2.image(mode2_bounds));
   

  return continuous_basic_set_type(integral_set.centre(),integral_set.generators()+LinearAlgebra::outer_product(v,tgrad));
}


template<class R>
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type>
Evaluation::ConstraintHybridEvolver<R>::lower_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                             const timed_set_type& initial_set,
                                                             time_type& step_size)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type> 
Evaluation::ConstraintHybridEvolver<R>::upper_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                             const timed_set_type& initial_set,
                                                             time_type& step_size)
{
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;;

      
  std::vector<timed_set_type> result;

  const id_type& discrete_state = initial_set.discrete_state();
  const continuous_basic_set_type& basic_set = initial_set.continuous_state_set();
  const time_type& time = initial_set.time();
  const discrete_time_type& steps = initial_set.steps();

  const mode_type& mode = automaton.mode(discrete_state);
  const System::VectorFieldInterface<R>& dynamic=mode.dynamic();

  const std::vector< constraint_const_pointer >& invariants=automaton.invariants(discrete_state);
  const std::map< id_type, constraint_const_pointer >& activations=automaton.activations(discrete_state);
  const std::map< id_type, constraint_const_pointer >& guards=automaton.guards(discrete_state);

  Rectangle<R> bounding_box=this->_integrator->estimate_flow_bounds(dynamic,basic_set.bounding_box(),step_size);
  continuous_basic_set_type reach_set = this->_integrator->bounded_reachability_step(dynamic,basic_set,bounding_box,step_size);
  continuous_basic_set_type integral_set = this->_integrator->bounded_integration_step(dynamic,basic_set,bounding_box,step_size);
  ARIADNE_LOG(4,"step_size="<<step_size<<", time="<<time<<"\n");
  ARIADNE_LOG(4,"reach_set="<<reach_set<<"\n");
  ARIADNE_LOG(4,"integral_set="<<integral_set<<"\n");

  if(step_size < this->_integrator->minimum_step_size()) {
    throw std::runtime_error("Minimum step size reached -- aborting.");
  }

  bool block = false;

  // Check for blocking by invariants
  ARIADNE_LOG(7,"invariants.size()="<<invariants.size()<<"\n");
  for(typename std::vector< constraint_const_pointer >::const_iterator constraint_ptr_iter = invariants.begin();
      constraint_ptr_iter != invariants.end(); ++constraint_ptr_iter)
  {
    const ConstraintInterface<R>& constraint=**constraint_ptr_iter;
    ARIADNE_LOG(7,"constraint_ptr="<<*constraint_ptr_iter<<"\n");
    ARIADNE_LOG(7,"constraint="<<constraint<<"\n");
    if(!satisfies(integral_set,constraint)) {
      ARIADNE_LOG(7,"  constraint not satisfied\n");
      block=true;
      break;
    }
  }

  // Check for blocking by guards (blocking if constraint is satisfied)
  ARIADNE_LOG(7,"guards.size()="<<guards.size()<<"\n");
  for(typename std::map< id_type, constraint_const_pointer >::const_iterator constraint_iter = guards.begin();
      constraint_iter != guards.end(); ++constraint_iter)
  {
    const ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(satisfies(integral_set,constraint)) {
      block=true;
      break;
    }
  }

  if(block==false) {
    result.push_back(timed_set_type(time+step_size,steps,discrete_state,integral_set));
  }


  // Check for activation 
  ARIADNE_LOG(7,"activations.size()="<<activations.size()<<"\n");
  for(typename std::map< id_type, constraint_const_pointer >::const_iterator constraint_iter = activations.begin();
      constraint_iter != activations.end(); ++constraint_iter)
  {
    id_type event_id = constraint_iter->first;
    const ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(!satisfies(reach_set,constraint)) {
    }
    else {
      const transition_type& transition = automaton.transition(event_id,discrete_state);
      const mode_type& destination = transition.destination();
      id_type destination_state = destination.id();
      const System::VectorFieldInterface<R>& destination_dynamic = destination.dynamic();
      const System::MapInterface<R>& reset = transition.reset();
      //continuous_basic_set_type jump_set = this->unforced_jump(dynamic,destination_dynamic,reset,basic_set,constraint,step_size);
      //result.push_back(timed_set_type(time+step_size,steps,destination_state,jump_set));
    }

  }

  // Check for guards 
  for(typename std::map< id_type, constraint_const_pointer >::const_iterator constraint_iter = guards.begin();
      constraint_iter != guards.end(); ++constraint_iter)
  {
    id_type event_id = constraint_iter->first;
    const ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(!satisfies(reach_set,constraint)) {
    }
    else {
      const transition_type& transition = automaton.transition(event_id,discrete_state);
      const mode_type& destination = transition.destination();
      id_type destination_state = destination.id();
      const System::VectorFieldInterface<R>& destination_dynamic = destination.dynamic();
      const System::MapInterface<R>& reset = transition.reset();
      //continuous_basic_set_type jump_set = this->forced_jump(dynamic,destination_dynamic,reset,basic_set,constraint,step_size);
      //result.push_back(timed_set_type(time+step_size,steps,destination_state,jump_set));
    }

  }

  return result;
}



template<class R>
std::vector<typename Evaluation::ConstraintHybridEvolver<R>::timed_set_type> 
Evaluation::ConstraintHybridEvolver<R>::upper_reachability_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                                const timed_set_type& initial_set,
                                                                time_type& step_size)
{
  using namespace LinearAlgebra;
  using namespace Geometry;
  using namespace System;;

      
  std::vector<timed_set_type> result;

  const id_type& discrete_state = initial_set.discrete_state();
  const continuous_basic_set_type& basic_set = initial_set.continuous_state_set();
  const time_type& time = initial_set.time();
  const discrete_time_type& steps = initial_set.steps();

  const mode_type& mode = automaton.mode(discrete_state);
  const System::VectorFieldInterface<R>& dynamic=mode.dynamic();

  const std::vector< constraint_const_pointer >& invariants=automaton.invariants(discrete_state);
  const std::map< id_type, constraint_const_pointer >& activations=automaton.activations(discrete_state);
  const std::map< id_type, constraint_const_pointer >& guards=automaton.guards(discrete_state);

  Rectangle<R> bounding_box=this->_integrator->estimate_flow_bounds(dynamic,basic_set.bounding_box(),step_size);
  continuous_basic_set_type reach_set = this->_integrator->bounded_reachability_step(dynamic,basic_set,bounding_box,step_size);
  continuous_basic_set_type integral_set = this->_integrator->bounded_integration_step(dynamic,basic_set,bounding_box,step_size);
  ARIADNE_LOG(4,"step_size="<<step_size<<", time="<<time<<"\n");
  ARIADNE_LOG(4,"reach_set="<<reach_set<<"\n");
  ARIADNE_LOG(4,"integral_set="<<integral_set<<"\n");

  if(step_size < this->_integrator->minimum_step_size()) {
    throw std::runtime_error("Minimum step size reached -- aborting.");
  }

  bool block = false;

  // Check for blocking by invariants
  ARIADNE_LOG(7,"invariants.size()="<<invariants.size()<<"\n");
  for(typename std::vector< constraint_const_pointer >::const_iterator constraint_ptr_iter = invariants.begin();
      constraint_ptr_iter != invariants.end(); ++constraint_ptr_iter)
  {
    const ConstraintInterface<R>& constraint=**constraint_ptr_iter;
    ARIADNE_LOG(7,"constraint_ptr="<<*constraint_ptr_iter<<"\n");
    ARIADNE_LOG(7,"constraint="<<constraint<<"\n");
    if(!satisfies(integral_set,constraint)) {
      ARIADNE_LOG(7,"  constraint not satisfied\n");
      block=true;
      break;
    }
  }

  // Check for blocking by guards (blocking if constraint is satisfied)
  ARIADNE_LOG(7,"guards.size()="<<guards.size()<<"\n");
  for(typename std::map< id_type, constraint_const_pointer >::const_iterator constraint_iter = guards.begin();
      constraint_iter != guards.end(); ++constraint_iter)
  {
    const ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(satisfies(integral_set,constraint)) {
      block=true;
      break;
    }
  }

  if(block==false) {
    result.push_back(timed_set_type(time+step_size,steps,discrete_state,integral_set));
  }


  // Check for activation 
  ARIADNE_LOG(7,"activations.size()="<<activations.size()<<"\n");
  for(typename std::map< id_type, constraint_const_pointer >::const_iterator constraint_iter = activations.begin();
      constraint_iter != activations.end(); ++constraint_iter)
  {
    id_type event_id = constraint_iter->first;
    const ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(!satisfies(reach_set,constraint)) {
    }
    else {
      const transition_type& transition = automaton.transition(event_id,discrete_state);
      const mode_type& destination = transition.destination();
      id_type destination_state = destination.id();
      const System::VectorFieldInterface<R>& destination_dynamic = destination.dynamic();
      const System::MapInterface<R>& reset = transition.reset();
      continuous_basic_set_type jump_set = this->unforced_jump(dynamic,destination_dynamic,reset,basic_set,constraint,step_size);
      result.push_back(timed_set_type(time+step_size,steps,destination_state,jump_set));
      assert(false); // FINISH
    }

  }

  // Check for guards 
  for(typename std::map< id_type, constraint_const_pointer >::const_iterator constraint_iter = guards.begin();
      constraint_iter != guards.end(); ++constraint_iter)
  {
    id_type event_id = constraint_iter->first;
    const ConstraintInterface<R>& constraint=*constraint_iter->second;
    if(!satisfies(reach_set,constraint)) {
    }
    else {
      const transition_type& transition = automaton.transition(event_id,discrete_state);
      const mode_type& destination = transition.destination();
      id_type destination_state = destination.id();
      const System::VectorFieldInterface<R>& destination_dynamic = destination.dynamic();
      const System::MapInterface<R>& reset = transition.reset();
      continuous_basic_set_type jump_set = this->forced_jump(dynamic,destination_dynamic,reset,basic_set,constraint,step_size);
      result.push_back(timed_set_type(time+step_size,steps,destination_state,jump_set));
    }

  }

  return result;
}



template<class R>
Geometry::HybridListSet<typename Evaluation::ConstraintHybridEvolver<R>::continuous_basic_set_type> 
Evaluation::ConstraintHybridEvolver<R>::upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                     const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                     time_type evolution_time,
                                                     size_type maximum_number_of_events)
{
  ARIADNE_LOG(2,"HybridEvolver::upper_evolve(HybridAutomaton automaton, ListSet initial_set, Time time, Integer maximum_number_of_events)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  using namespace Numeric;
  
  time_type maximum_step_size=this->_integrator->maximum_step_size();
  R maximum_basic_set_radius=this->_integrator->maximum_basic_set_radius();
  ARIADNE_LOG(7,"maximum_step_size="<<maximum_step_size<<", maximum_basic_set_radius="<<maximum_basic_set_radius<<"\n");
  
  time_type step_size=maximum_step_size;
  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets;
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  Geometry::HybridListSet< continuous_basic_set_type > final_set(automaton.locations());

  for(typename Geometry::HybridListSet<continuous_basic_set_type>::const_iterator loc_iter=initial_set.begin();
      loc_iter!=initial_set.end(); ++loc_iter)
  {
    id_type loc_id=loc_iter->first;
    const Geometry::ListSet<continuous_basic_set_type>& list_set=*loc_iter->second;
    for(typename Geometry::ListSet<continuous_basic_set_type>::const_iterator bs_iter=list_set.begin();
        bs_iter!=list_set.end(); ++bs_iter)
    {
      const continuous_basic_set_type& bs=*bs_iter;
      working_sets.push_back(timed_set_type(0,0,loc_id,bs));
    }
  }


  typedef typename std::vector< timed_set_type >::const_iterator list_set_const_iterator;
  while(!working_sets.empty()) {
    timed_set_type timed_set=working_sets.back(); working_sets.pop_back();
    
    step_size=Numeric::max(time_type(2*step_size),maximum_step_size);

    if(timed_set.time()==evolution_time) {
      final_set.adjoin(hybrid_basic_set_type(timed_set.discrete_state(),timed_set.continuous_state_set()));
    } else if(timed_set.continuous_state_set().radius()>maximum_basic_set_radius) {
      ::append(working_sets,this->subdivide(timed_set));
    } else {
      step_size=Numeric::max(step_size,time_type(evolution_time-timed_set.time()));
      ::append(working_sets,this->upper_evolution_step(automaton,timed_set,step_size));
    }
  }

  return final_set;
}


template<class R>
Geometry::HybridListSet<typename Evaluation::ConstraintHybridEvolver<R>::continuous_basic_set_type> 
Evaluation::ConstraintHybridEvolver<R>::upper_reach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                    time_type evolution_time, 
                                                    size_type maximum_number_of_events)
{
  ARIADNE_LOG(2,"HybridEvolver::upper_evolve(HybridAutomaton automaton, ListSet initial_set, Time time, Integer maximum_number_of_events)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  using namespace Numeric;
  
  time_type maximum_step_size=this->_integrator->maximum_step_size();
  R maximum_basic_set_radius=this->_integrator->maximum_basic_set_radius();
  ARIADNE_LOG(7,"maximum_step_size="<<maximum_step_size<<", maximum_basic_set_radius="<<maximum_basic_set_radius<<"\n");
  
  time_type step_size=maximum_step_size;
  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets;
  std::vector< timed_set_type > reached_sets;
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  Geometry::HybridListSet< continuous_basic_set_type > reach_set(automaton.locations());

  for(typename Geometry::HybridListSet<continuous_basic_set_type>::const_iterator loc_iter=initial_set.begin();
      loc_iter!=initial_set.end(); ++loc_iter)
  {
    id_type loc_id=loc_iter->first;
    const Geometry::ListSet<continuous_basic_set_type>& list_set=*loc_iter->second;
    for(typename Geometry::ListSet<continuous_basic_set_type>::const_iterator bs_iter=list_set.begin();
        bs_iter!=list_set.end(); ++bs_iter)
    {
      const continuous_basic_set_type& bs=*bs_iter;
      working_sets.push_back(timed_set_type(0,0,loc_id,bs));
    }
  }


  typedef typename std::vector< timed_set_type >::const_iterator list_set_const_iterator;
  while(!working_sets.empty()) {
    timed_set_type timed_set=working_sets.back(); working_sets.pop_back();
    
    step_size=Numeric::max(time_type(2*step_size),maximum_step_size);

    if(timed_set.time()==evolution_time) {
    } else if(timed_set.continuous_state_set().radius()>maximum_basic_set_radius) {
      ::append(working_sets,this->subdivide(timed_set));
    } else {
      step_size=Numeric::max(step_size,time_type(evolution_time-timed_set.time()));
      ::append(working_sets,this->upper_evolution_step(automaton,timed_set,step_size));
      ::append(reached_sets,this->upper_reachability_step(automaton,timed_set,step_size));
    }
  }

  for(typename std::vector<timed_set_type>::const_iterator reach_iter=reached_sets.begin();
      reach_iter!=reached_sets.end(); ++reach_iter)
  {
    reach_set.adjoin(hybrid_basic_set_type(reach_iter->discrete_state(),reach_iter->continuous_state_set()));
  }

  return reach_set;
}



      

} // namespace Ariadne
