/***************************************************************************
 *            hybrid_evolver-image.cc
 *
 *  Copyright  2008-10  Alberto Casagrande, Pieter Collins, Luca Geretti
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

#include "macros.h"
#include "array.h"
#include "tuple.h"
#include "stlio.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "taylor_set.h"
#include "taylor_function.h"
#include "taylor_model.h"
#include "orbit.h"
#include "taylor_calculus.h"
#include "settings.h"

#include "logging.h"

#include "hybrid_time.h"
#include "hybrid_automaton.h"
#include "hybrid_evolver-image.h"

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

void wait_for_keypress() {
    std::string str;
    getline(std::cin,str);
}


static const int BLOCKING_EVENT = -2;

class DegenerateCrossingException : public std::runtime_error {
  public:
    DegenerateCrossingException(const char* msg) : std::runtime_error(msg) { }
};

const DiscreteEvent ImageSetHybridEvolver::starting_event = -1;
const DiscreteEvent ImageSetHybridEvolver::finishing_event = -2;
const DiscreteEvent ImageSetHybridEvolver::blocking_event = -3;

typedef VectorFunction RealVectorFunction;
typedef ScalarFunction RealScalarFunction;

ImageSetHybridEvolver::ImageSetHybridEvolver(const SystemType& system)
    : EvolverBase(system),
      _settings(new EvolutionSettingsType()),
      _upper_toolbox(new TaylorCalculus()),
      _lower_toolbox(new TaylorCalculus())
{
	
}

ImageSetHybridEvolver::ImageSetHybridEvolver(
		const SystemType& system,
		const TaylorCalculus& upper_calculus,
		const TaylorCalculus& lower_calculus)
    : EvolverBase(system),
      _settings(new EvolutionSettingsType()),
      _upper_toolbox(new TaylorCalculus(upper_calculus)),
      _lower_toolbox(new TaylorCalculus(lower_calculus))
{
}

void
ImageSetHybridEvolver::
tune_settings(
		const HybridGrid& grid,
		const HybridFloatVector& hmad,
		AccuracyType accuracy,
		Semantics semantics)
{
	this->_settings->maximum_enclosure_cell = getMaximumEnclosureCell(grid,accuracy);
	ARIADNE_LOG(2, "Maximum enclosure cell: " << this->_settings->maximum_enclosure_cell << "\n");
	this->_settings->hybrid_maximum_step_size = getHybridMaximumStepSize(hmad,grid,accuracy,semantics);
	ARIADNE_LOG(2, "Maximum step size: " << this->_settings->hybrid_maximum_step_size << "\n");
}


void
ImageSetHybridEvolver::
_evolution(EnclosureListType& final_sets,
           EnclosureListType& reach_sets,
           EnclosureListType& intermediate_sets,
           const EnclosureType& initial_set,
           const TimeType& maximum_hybrid_time,
           bool ignore_activations,
           ContinuousEvolutionDirection direction,
           Semantics semantics) const
{
    ARIADNE_LOG(5,ARIADNE_PRETTY_FUNCTION<<"\n");
    ARIADNE_LOG(1,"Computing evolution up to "<<maximum_hybrid_time.continuous_time()<<" time units and "<<maximum_hybrid_time.discrete_time()<<" steps.\n");

	// The working sets, pushed back and popped front
    std::list< HybridTimedSetType > working_sets;
    _evolution_add_initialSet(working_sets,initial_set,semantics);

	// While there exists a working set, process it and increment the total
    uint i=0;
	while(!working_sets.empty()) {

		ARIADNE_LOG(2,"\n");
		ARIADNE_LOG(2,"Processed sets: " << i++ << ", remaining sets: " << working_sets.size() << "\n\n");

		// Get the least recent working set, pop it and update the corresponding size
		HybridTimedSetType current_set = working_sets.front();
		working_sets.pop_front();

		// Get the members of the current set
        DiscreteLocation initial_location=current_set.first;
        EventListType initial_events=current_set.second;
		SetModelType initial_set_model=current_set.third;
		TimeModelType initial_time_model=current_set.fourth;

		bool isEnclosureTooLarge = _is_enclosure_too_large(initial_set_model);
		bool subdivideOverTime = (initial_time_model.range().width() > this->_settings->hybrid_maximum_step_size[initial_location]/2);

		if(initial_time_model.range().lower()>=maximum_hybrid_time.continuous_time() ||
		   initial_events.size()>=uint(maximum_hybrid_time.discrete_time())) {
            ARIADNE_LOG(3,"Final time reached, adjoining result to final sets.\n");
            final_sets.adjoin(initial_location,this->getCalculusInterface(semantics).enclosure(initial_set_model));
        } else if (subdivideOverTime && this->_settings->enable_subdivisions) {
            ARIADNE_LOG(1,"WARNING: computed time range " << initial_time_model.range() << " width larger than half the maximum step size " << this->_settings->hybrid_maximum_step_size[initial_location] << ", subdividing over time.\n");
            _add_models_subdivisions_time(working_sets,initial_set_model,initial_time_model,initial_location,initial_events,semantics);
		} else if (semantics == UPPER_SEMANTICS && this->_settings->enable_subdivisions && isEnclosureTooLarge) {
            ARIADNE_LOG(1,"WARNING: computed set range " << initial_set_model.range() << " widths larger than maximum_enclosure_cell " << this->_settings->maximum_enclosure_cell << ", subdividing.\n");
            _add_models_subdivisions_autoselect(working_sets,initial_set_model,initial_time_model,initial_location,initial_events,semantics);
        } else if((semantics == LOWER_SEMANTICS || !this->_settings->enable_subdivisions) &&
                  this->_settings->enable_premature_termination_on_enclosure_size && isEnclosureTooLarge) {
            ARIADNE_LOG(1,"\n\nWARNING: Terminating evolution at time " << initial_time_model.value()
                        << " and set " << initial_set_model.centre() << " due to maximum enclosure bounds being exceeded.\n\n");
            if(semantics == UPPER_SEMANTICS)
                final_sets.adjoin(initial_location,this->getCalculusInterface(semantics).enclosure(initial_set_model));
        } else {
            this->_evolution_step(working_sets,reach_sets,intermediate_sets,current_set,maximum_hybrid_time,
            		ignore_activations,direction,semantics);
        }

		_logStepAtVerbosity1(working_sets,reach_sets,initial_events,initial_time_model,initial_set_model,initial_location);
    }
}


void
ImageSetHybridEvolver::
_evolution_step(std::list< HybridTimedSetType >& working_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const HybridTimedSetType& working_set,
                const TimeType& maximum_hybrid_time,
                bool ignore_activations,
                ContinuousEvolutionDirection direction,
                Semantics semantics) const
{
    // Use the following basic algorithm for computing an evolution step
    //
    // 1) Extract data about working set and location
    //
    // 2) Find all blocking events (invariants and urgent transitions) which are
    //    active at the starting time. If any events are definitely active,
    //    the corresponding reset occurs, and there is no continuous evolution.
    //    If there are possibly active events, these occur along with continuous
    //    evolution for upper semantics, and evolution terminates with lower
    //    semantics.
    //
    //    Input: Set: starting_set, map<Event,Expression> guards
    //    Output: map<Event,Tribool> initially_active
    //         (includes special "blocking event")
    //
    // 3) Compute the continuous evolution for a fixed step size h, over the
    //    time interval [-h,+h]
    //
    //    Input: Function: dynamic,
    //           Box: domain OR SetModel initial_set,
    //           Time: maximum_step_size
    //    Output: FunctionModel: flow_model
    //               OR SetModel: flow_set_model
    //            Box: flow_bounds
    //
    // 4) For each blocking event, compute the crossing time with the guard set.
    //    The computed crossing time may lie outside the flow time interval;
    //    any such crossings will be ignored. Non-transverse crossings may have
    //    large crossing time intervals.
    //
    //    If the transition is definitely not initially active and the crossing
    //    is transverse, then there are no problems. If the transition is
    //    possibly initially active, then the crossing may be in the "wrong"
    //    direction, i.e. the transition may become inactive. In this case, we
    //    have upper semantics (otherwise we would already have terminated the
    //    evolution) and the transition is considered inactive for evolution
    //    purposes.
    //
    //    If the crossing is not transverse, then lower evolution is blocking
    //    and upper evolution considers the transition as non-urgent. However,
    //    we should probably evolve close to the transition, and maybe split
    //    the evolution across the transversality boundary.
    //
    //    Input: map<Event,Expression>: resets,
    //           Box: flow_bounds OR FunctionModel: flow_model
    //             OR SetModel: flow_set_model
    //    Output: map<Event,ExpressionModel> crossing_times,
    //            set<Event> tangential_events
    //
    // 5) Compute the blocking time and blocking events. The blocking time is
    //    the minimum of the computed crossing times.
    //
    //    If there is a single  blocking event, the crossing is tranverse, and
    //    occurs between the starting and finishing times, then evolution
    //    proceeds according to this event.
    //
    //    If there are multiple blocking events, and the upper bound of the
    //    crossing time range is large, then we set the finishing time to
    //    just below the crossing time. This means that in the next step
    //    we may be better able to resolve the crossing.
    //
    //    If there are multiple blocking events and the upper bound of the
    //    blocking time range is small, then lower evolution terminates, and
    //    upper evolution proceeds according to the crossing time of each
    //    blocking event.
    //
    //    Input: map<Event,ExpressionModel>: crossing_times
    //    Output: set<Event> blocking_events, ExpressionModel: blocking_time
    //
    // 6) Compute the initial and final activation times of the non-blocking
    //    events. Tangential crossings are included in this computation,
    //    as they are treated as non-urgent.
    //
    //    Compute the maximum of the initial activation time and starting time,
    //    and the minimum of the final activation time and finishing/blocking
    //    time.
    //
    //    Input: map<Event,Expression>: guards,
    //           Box: domain OR FunctionModel: flow_model
    //               OR SetModel: flow_set_model
    //    Output: map<Event,(ExpressionModel,ExpressionModel)>: initial/final times
    //
    // 7) Apply the flows, guards and resets according to the computed
    //    event times.


    const double SMALL_RELATIVE_TIME=1./16;

    // Extract information about the working set
    DiscreteLocation location(0);
    IntegerType steps;
    EventListType events_history;
    SetModelType set_model;
    TimeModelType time_model;
    ARIADNE_LOG(9,"working_set = "<<working_set<<"\n");
    make_ltuple(location,events_history,set_model,time_model)=working_set;
    steps=events_history.size();

    // Extract information about the current location
    const RealVectorFunction dynamic=get_directed_dynamic(_sys->dynamic_function(location),direction);
    Set<DiscreteEvent> available_events = _sys->events(location);
    std::map<DiscreteEvent,RealScalarFunction> urgent_guards, permissive_guards, invariants;

    for (Set<DiscreteEvent>::const_iterator event_it = available_events.begin(); event_it != available_events.end(); ++event_it) {
    	EventKind kind = _sys->event_kind(location,*event_it);
    	switch (kind) {
    		case INVARIANT:
    			invariants.insert(std::pair<DiscreteEvent,RealScalarFunction>(*event_it,_sys->invariant_function(location,*event_it)));
    			urgent_guards.insert(std::pair<DiscreteEvent,RealScalarFunction>(*event_it,_sys->invariant_function(location,*event_it)));
    			break;
    		case URGENT:
    			urgent_guards.insert(std::pair<DiscreteEvent,RealScalarFunction>(*event_it,_sys->guard_function(location,*event_it)));
    			break;
    		case PERMISSIVE:
    			permissive_guards.insert(std::pair<DiscreteEvent,RealScalarFunction>(*event_it,_sys->guard_function(location,*event_it)));
    			break;
    		default:
    			ARIADNE_FAIL_MSG("Unhandled event kind.");
    	}
    }

    // Check to make sure dimensions are correct
    ARIADNE_ASSERT(set_model.argument_size()==time_model.argument_size());
    ARIADNE_ASSERT_MSG(set_model.result_size()==_sys->dimension(location),"set_model="<<set_model<<", location="<<location.name());

    _logEvolutionStepInitialState(events_history,time_model,location,set_model,dynamic,invariants,urgent_guards,permissive_guards);

    if (is_enclosure_to_be_discarded(set_model,urgent_guards,dynamic,semantics))
    	return;

    // Compute continuous evolution
    FlowSetModelType flow_set_model; BoxType flow_bounds; 
    Float time_step = this->_settings->hybrid_maximum_step_size[location];
    const Float maximum_time=maximum_hybrid_time.continuous_time();
    compute_flow_model(flow_set_model,flow_bounds,time_step,dynamic,set_model,time_model,maximum_time,semantics);

    ARIADNE_LOG(2,"flow_bounds = "<<flow_bounds<<"\n")
    ARIADNE_LOG(2,"time_step = "<<time_step<<"\n")
    ARIADNE_LOG(2,"flow_range = "<<flow_set_model.range()<<"\n");
    ARIADNE_LOG(2,"starting_set_range = "<<set_model.range()<<"\n");
    // Partial evaluation on flow set model to obtain final set must take scaled time equal to 1.0
    SetModelType finishing_set=partial_evaluate(flow_set_model.models(),set_model.argument_size(),1.0);
    ARIADNE_LOG(2,"finishing_set_range = "<<finishing_set.range()<<"\n")

    // Set special events and times; note that the time step is scaled to [0,1]
    TimeModelType zero_time_model = this->getCalculusInterface(semantics).time_model(0.0,Box(time_model.argument_size()));
    TimeModelType time_step_model = this->getCalculusInterface(semantics).time_model(1.0,Box(time_model.argument_size()));

    std::set<DiscreteEvent> blocking_events;
    TimeModelType blocking_time_model;
    std::set<DiscreteEvent> non_transverse_events;
    _compute_blocking_info(non_transverse_events,blocking_events,blocking_time_model,
    				  time_step_model,flow_set_model,urgent_guards,SMALL_RELATIVE_TIME,semantics);

    ActivationTimesType activation_times;
    if (!ignore_activations) {
    	_compute_activation_info(permissive_guards,activation_times,non_transverse_events,
    			flow_set_model,blocking_time_model,urgent_guards,invariants,semantics);
    }

    SetModelType reachable_set;
    _compute_and_adjoin_reachableSet(reach_sets,reachable_set,location,flow_set_model,zero_time_model,blocking_time_model,semantics);

    if(semantics!=LOWER_SEMANTICS || blocking_events.size()==1)
    	_computeEvolutionForEvents(working_sets,intermediate_sets,location,blocking_events,events_history,
    								activation_times,flow_set_model,time_model,blocking_time_model,time_step,ignore_activations,semantics);

}


ImageSetHybridEvolver::TimeModelType ImageSetHybridEvolver::
crossing_time(
		VectorFunction guard,
		const FlowSetModelType& flow_set_model,
		Semantics semantics) const
{
    try {
        TimeModelType crossing_time_model=this->getCalculusInterface(semantics).scaled_crossing_time(guard,flow_set_model);
        return crossing_time_model;
    }
    catch(DegenerateCrossingException e) {
        BoxType space_domain=project(flow_set_model.domain(),range(0,flow_set_model.argument_size()-1));
        Interval touching_time_interval=this->getCalculusInterface(semantics).scaled_touching_time_interval(guard,flow_set_model);
        TimeModelType touching_time_model=this->getCalculusInterface(semantics).time_model(touching_time_interval,space_domain);
        return touching_time_model;
    } // end non-transverse crossing
}


Interval jacobian2_range(const TaylorModel& tm);

Interval ImageSetHybridEvolver::
normal_derivative(VectorFunction guard, const FlowSetModelType& flow_set_model, const TimeModelType& crossing_time_model) const
{
    typedef TimeModelType GuardValueModelType;
    GuardValueModelType guard_flow_set_model=apply(guard,flow_set_model)[0];
    Interval normal_derivative=jacobian2_range(guard_flow_set_model);
    return normal_derivative;
}


void ImageSetHybridEvolver::
compute_initially_active_events(
		std::map<DiscreteEvent,tribool>& initially_active_events,
        const std::map<DiscreteEvent,RealScalarFunction>& urgent_guards,
        const ContinuousEnclosureType& initial_set,
        Semantics semantics) const
{
	ARIADNE_LOG(2,"Computing initially active events...\n");
    tribool blocking_event_initially_active=false;
    for(std::map<DiscreteEvent,RealScalarFunction>::const_iterator iter=urgent_guards.begin(); iter!=urgent_guards.end(); ++iter) {
        RealVectorFunction activation(1,iter->second);
        ARIADNE_LOG(3,"Evaluating urgent guard '" << iter->first.name() << "' with activation " << activation << ": ");
        tribool initially_active=this->getCalculusInterface(semantics).active(activation,initial_set);
        if(possibly(initially_active)) {
        	ARIADNE_LOG(3," possibly active\n");
            initially_active_events.insert(std::make_pair(iter->first,initially_active));
            blocking_event_initially_active = blocking_event_initially_active || initially_active;
        } else {
        	ARIADNE_LOG(3," inactive\n");
        }
    }
    initially_active_events.insert(std::make_pair(blocking_event,blocking_event_initially_active));
    return;
}

bool
ImageSetHybridEvolver::
has_nonnegative_crossing(
		const std::map<DiscreteEvent,RealScalarFunction>& urgent_guards,
		const RealVectorFunction dynamic,
		const Box& set_bounds,
		Semantics semantics) const
{
	ARIADNE_LOG(3,"Checking nonnegative crossings...\n");
    for(std::map<DiscreteEvent,RealScalarFunction>::const_iterator iter=urgent_guards.begin(); iter!=urgent_guards.end(); ++iter) {
        RealVectorFunction activation(1,iter->second);
        ARIADNE_LOG(3,"Guard: " << activation << "\n");
        tribool is_active = this->getCalculusInterface(semantics).active(activation,set_bounds);
        ARIADNE_LOG(3,"Active: " << is_active);
        tribool is_positively_crossing = positively_crossing(set_bounds,dynamic,activation[0]);
        ARIADNE_LOG(3,"; positively crossing: " << is_positively_crossing << "\n");
        if(possibly(is_active) && possibly(is_positively_crossing))
			return true;
    }
    return false;
}



tribool
positively_crossing(const Box& set_bounds,
					const RealVectorFunction& dynamic,
					const RealScalarFunction& activation)
{
    RealScalarFunction derivative=lie_derivative(activation,dynamic);
    Interval derivative_range = derivative.evaluate(set_bounds);

    if (derivative_range.lower() > 0)
    	return true;
    else if (derivative_range.upper() < 0)
    	return false;
    else return indeterminate;
}

bool
ImageSetHybridEvolver::
is_enclosure_to_be_discarded(const ContinuousEnclosureType& enclosure,
        					 const std::map<DiscreteEvent,RealScalarFunction>& urgent_guards,
        					 const RealVectorFunction& dynamic,
							 Semantics semantics) const
{
	bool result = false;

	std::map<DiscreteEvent,tribool> initially_active_events;
	this->compute_initially_active_events(initially_active_events, urgent_guards, enclosure, semantics);
	ARIADNE_LOG(2,"initially_active_events = "<<initially_active_events<<"\n\n");

	tribool has_any_initially_active_blocking_event = initially_active_events[blocking_event];

	// Test for initially active events, and process these as required
	if(definitely(has_any_initially_active_blocking_event)) {
		ARIADNE_LOG(2,"An invariant is definitely initially active: discarding the set.\n");
		result = true;
	} else if(possibly(has_any_initially_active_blocking_event) && semantics==LOWER_SEMANTICS) {
		ARIADNE_LOG(2,"A blocking event is possibly active: checking whether there is a possibly positive crossing.\n");
		bool has_nonneg_crossing = has_nonnegative_crossing(urgent_guards,dynamic,enclosure.bounding_box(),semantics);
		if (has_nonneg_crossing) {
			ARIADNE_LOG(2,"Terminating lower evolution due to possibly initially active invariant with nonnegative crossing.\n");
			result = true;
		}
	}

	return result;
}

// Compute the flow, parameterising space with the set parameters
void ImageSetHybridEvolver::
compute_flow_model(
		FlowSetModelType& flow_set_model,
		BoxType& flow_bounds, Float& time_step,
        VectorFunction dynamic,
        const SetModelType& starting_set_model,
        const TimeModelType& starting_time_model,
        Float finishing_time,
        Semantics semantics) const
{
    ARIADNE_LOG(3,"compute_flow_model(....)\n");
    const int MAXIMUM_BOUNDS_DIAMETER_FACTOR = 8;
    float remaining_time = finishing_time - starting_time_model.range().lower();
    const Float maximum_step_size=min(time_step, remaining_time);
    const Float maximum_bounds_diameter=max(this->_settings->maximum_enclosure_cell)*MAXIMUM_BOUNDS_DIAMETER_FACTOR;

    BoxType starting_set_bounding_box=starting_set_model.range();
    ARIADNE_LOG(3,"starting_set_bounding_box="<<starting_set_bounding_box<<"\n");
    make_lpair(time_step,flow_bounds)=this->getCalculusInterface(semantics).flow_bounds(dynamic,starting_set_bounding_box,maximum_step_size,maximum_bounds_diameter);
    // Compute the flow model
    ARIADNE_LOG(3,"time_step="<<time_step<<"\n");
    ARIADNE_LOG(3,"flow_bounds="<<flow_bounds<<"\n");
    FlowModelType flow_model=this->getCalculusInterface(semantics).flow_model(dynamic,starting_set_bounding_box,time_step,flow_bounds);
    ARIADNE_LOG(3,"flow_model="<<flow_model<<"\n");
    ScalarTaylorFunction identity_time_expression=ScalarTaylorFunction::variable(BoxType(1u,Interval(-time_step,+time_step)),0u);
    flow_set_model=unchecked_apply(flow_model,combine(starting_set_model.models(),identity_time_expression.model()));
}



void ImageSetHybridEvolver::
compute_eventBlockingTimes_and_nonTransverseEvents(
		std::map<DiscreteEvent,TimeModelType>& event_blocking_times,
        std::set<DiscreteEvent>& non_transverse_events,
        const std::map<DiscreteEvent,RealScalarFunction>& urgent_guards,
        const FlowSetModelType& flow_set_model,
        Semantics semantics) const
{
    ARIADNE_LOG(3,"Computing blocking events.\n");

    uint dimension=flow_set_model.result_size();
    const double SMALL_RELATIVE_TIME = 1./16;
    FlowSetModelType positive_flow_set_model(split(flow_set_model.models(),dimension,1));

    for(std::map<DiscreteEvent,RealScalarFunction>::const_iterator guard_iter=urgent_guards.begin();
        guard_iter!=urgent_guards.end(); ++guard_iter)
    {
        const DiscreteEvent event=guard_iter->first;
        const RealVectorFunction guard(1,guard_iter->second);
        tribool active = this->getCalculusInterface(semantics).active(guard,positive_flow_set_model);
        if(possibly(active)) {
            ARIADNE_LOG(3,"Event "<<event<<" possibly active.\n");
            TimeModelType crossing_time_model;
            Interval normal_derivative;
            try {
                crossing_time_model=this->getCalculusInterface(semantics).scaled_crossing_time(guard,flow_set_model);
                normal_derivative=this->normal_derivative(guard,flow_set_model,crossing_time_model);
                assert(normal_derivative.lower()>0 || normal_derivative.upper()<0);
                if(normal_derivative.lower()>0) {
                    ARIADNE_LOG(3,"Event "<<event<<" inserted into blocking times.\n");
                    event_blocking_times[event]=crossing_time_model;
                }
            }
            catch(DegenerateCrossingException e) {
                ARIADNE_LOG(3,"Degenerate Crossing exception catched.\n");
                BoxType space_domain=project(flow_set_model.domain(),range(0,flow_set_model.argument_size()-1));
                Interval touching_time_interval=this->getCalculusInterface(semantics).scaled_touching_time_interval(guard,flow_set_model);
                TimeModelType touching_time_model=this->getCalculusInterface(semantics).time_model(touching_time_interval,space_domain);
                // Use 1.0 as upper bound above since flow set model has time interval normalised to [-1,+1]
                ARIADNE_LOG(3,"touching_time_interval="<<touching_time_interval<<"\n");
                if(touching_time_interval.upper()>=0 && touching_time_interval.lower()<=1.0) {
                    SetModelType finishing_set_model=partial_evaluate(flow_set_model.models(),dimension,1.0);
                    tribool finishing_set_active=this->getCalculusInterface(semantics).active(guard,finishing_set_model);
                    if(definitely(finishing_set_active)) {
                        ARIADNE_LOG(3,"Event is definitely finally active, inserting it into blocking times.\n");
                        event_blocking_times[event]=touching_time_model;
                    } else if(possibly(finishing_set_active)) {
                        ARIADNE_LOG(3,"Event is possibly finally active.\n");
                        if(touching_time_interval.lower()>SMALL_RELATIVE_TIME) {
                            ARIADNE_LOG(3,"lower touching time is greater than zero, inserting event into blocking times.\n");
                            TaylorModel lower_touching_time_model=
                            		this->getCalculusInterface(semantics).time_model(touching_time_interval.lower(),space_domain);
                            event_blocking_times[finishing_event]=lower_touching_time_model;
                        } else {
                            ARIADNE_LOG(3,"DANGER: we can't determine whether the crossing is completely finished or not..\n");
                            // FIXME: Here we are stuck, we can't determine whether the crossing is completely finished or not.
                            // Just put this in as a blocking event and hope for the best...
                            // event_blocking_times[event]=touching_time_model;
                            // DAVIDE: setting this as a blocking event is not correct,
                            //         because it causes continuous evolution to stop.
                            //         I think it is better to put it into non-transverse-events.
                            non_transverse_events.insert(event);
                        }
                    } else {
                        ARIADNE_LOG(3,"After the flow step, the event is not active again, so the crossing was tangential.\n");
                        // After the flow step, the event is not active again, so the crossing was tangential
                        non_transverse_events.insert(event);
                    }
                }
            } // end non-transverse crossing
        } // end possibly active
    } // end main loop

    return;
}


void ImageSetHybridEvolver::
compute_blockingTime_and_relatedEvents(std::set<DiscreteEvent>& blocking_events,
                      TimeModelType& blocking_time,
                      const std::map<DiscreteEvent,TimeModelType>& event_blocking_times) const
{
    assert(!event_blocking_times.empty());
    blocking_events.insert(event_blocking_times.begin()->first);
    blocking_time=event_blocking_times.begin()->second;

    for(std::map<DiscreteEvent,TimeModelType>::const_iterator iter=++event_blocking_times.begin();
        iter!=event_blocking_times.end(); ++iter)
    {
        const DiscreteEvent event=iter->first;
        const TimeModelType& time=iter->second;
        tribool is_first=(time<blocking_time);
        if(definitely(is_first)) {
            blocking_events.clear();
            blocking_events.insert(event);
            blocking_time=time;
        } else if(possibly(is_first)) {
            std::set<DiscreteEvent> new_blocking_events;
            new_blocking_events.insert(event);
            TimeModelType new_blocking_time=time;
            for(std::set<DiscreteEvent>::const_iterator iter=blocking_events.begin();
                iter!=blocking_events.end(); ++iter)
            {
                const TimeModelType& current_event_time=event_blocking_times.find(*iter)->second;
                if(possibly(current_event_time<time)) {
                    new_blocking_events.insert(*iter);
                    blocking_time=min(current_event_time,blocking_time);
                }
            }
            blocking_time.swap(new_blocking_time);
            blocking_events.swap(new_blocking_events);
        }
    }
    return;
}


void ImageSetHybridEvolver::
compute_activationTimes(std::map<DiscreteEvent,tuple<TimeModelType,TimeModelType> >& activation_times,
                         const std::map<DiscreteEvent,RealScalarFunction>& activations,
                         const FlowSetModelType& flow_set_model,
                         const TimeModelType& blocking_time_model,
                         const Semantics semantics) const
{
    SetModelType initial_set_model=partial_evaluate(flow_set_model.models(),flow_set_model.argument_size()-1,0.0);
    SetModelType final_set_model=partial_evaluate(flow_set_model.models(),flow_set_model.argument_size()-1,1.0);
    TimeModelType zero_time_model=blocking_time_model*0.0;

    for(std::map<DiscreteEvent,RealScalarFunction>::const_iterator iter=activations.begin(); iter!=activations.end(); ++iter) {
        DiscreteEvent event=iter->first;
        RealVectorFunction activation(1,iter->second);

        ARIADNE_LOG(3,"Computing activation time for event "<<event<<"...");

        // Compute whether the event might be enabled on the entire time interval
        tribool active=this->getCalculusInterface(semantics).active(activation,flow_set_model);

        if(definitely(active)) {
            // The event is enabled over the entire time interval
            ARIADNE_LOG(3,"event is enabled over the entire time interval.\n");
            activation_times.insert(make_pair(event,make_tuple(zero_time_model,blocking_time_model)));
        } else if(possibly(active)) {
            ARIADNE_LOG(3,"event is possibly enabled.\n");
            // Compute whether the event is enabled at the beginning and end of the time interval
            tribool initially_active=this->getCalculusInterface(semantics).active(activation,initial_set_model);
            tribool finally_active=this->getCalculusInterface(semantics).active(activation,final_set_model);

            TimeModelType crossing_time_model=this->crossing_time(activation,flow_set_model,semantics);

            TimeModelType lower_crossing_time_model=crossing_time_model-crossing_time_model.error();
            TimeModelType upper_crossing_time_model=crossing_time_model+crossing_time_model.error();
            lower_crossing_time_model.set_error(0);
            upper_crossing_time_model.set_error(0);

            TimeModelType lower_active_time_model, upper_active_time_model;

            // Determine lower activation time
            if(definitely(not(initially_active))) {
                ARIADNE_LOG(3,"event is definitely not initially active.\n");
                switch(semantics) {
                    case UPPER_SEMANTICS: lower_active_time_model=lower_crossing_time_model; break;
                    case LOWER_SEMANTICS: lower_active_time_model=upper_crossing_time_model; break;
                }
            } else if(definitely(initially_active)) {
                ARIADNE_LOG(3,"event is definitely initially active.\n");
                lower_active_time_model=zero_time_model;
            } else {
                ARIADNE_LOG(3,"Event is undeteriminate initially active.\n");
                switch(semantics) {
                    case UPPER_SEMANTICS: lower_active_time_model=zero_time_model; break;
                    case LOWER_SEMANTICS: lower_active_time_model=upper_crossing_time_model; break;
                }
            }

            // Compute upper activation time
            if(definitely(not(finally_active))) {
                ARIADNE_LOG(3,"Event is definitely not finally active.\n");
                switch(semantics) {
                    case UPPER_SEMANTICS: upper_active_time_model=upper_crossing_time_model; break;
                    case LOWER_SEMANTICS: upper_active_time_model=lower_crossing_time_model; break;
                }
            } else if(definitely(finally_active)) {
                ARIADNE_LOG(3,"Event is definitely finally active.\n");
                upper_active_time_model=blocking_time_model;
            } else {
                ARIADNE_LOG(3,"Event is undeteriminate finally active.\n");
                switch(semantics) {
                    case UPPER_SEMANTICS: upper_active_time_model=blocking_time_model; break;
                    case LOWER_SEMANTICS: upper_active_time_model=upper_crossing_time_model; break;
                }
            }

            // In case of lower semantics, event may not be
            switch(semantics) {
                case UPPER_SEMANTICS:
                    activation_times.insert(make_pair(event,make_tuple(lower_active_time_model,upper_active_time_model)));
                    break;
                case LOWER_SEMANTICS:
                    if(lower_active_time_model<upper_active_time_model) {
                        activation_times.insert(make_pair(event,make_tuple(lower_active_time_model,upper_active_time_model)));
                    }
                    break;
            }
        }
    }
}


void
ImageSetHybridEvolver::
_logStepAtVerbosity1(const std::list<HybridTimedSetType>& working_sets,
					 const EnclosureListType& reach_sets,
					 const EventListType& initial_events,
					 const TimeModelType& initial_time_model,
					 const SetModelType& initial_set_model,
					 const DiscreteLocation& initial_location) const
{
    if(verbosity==1) {
        ARIADNE_LOG(1,"#w="<<std::setw(4)<<working_sets.size()
                    <<"#r="<<std::setw(4)<<std::left<<reach_sets.size()
                    <<" s="<<std::setw(3)<<std::left<<initial_events.size()
                    <<" t="<<std::fixed<<initial_time_model.value()
                    <<" r="<<std::setw(7)<<initial_set_model.radius()
                    <<" l="<<std::setw(3)<<std::left<<initial_location
                    <<" c="<<initial_set_model.centre()
                    <<" e="<<initial_events
                    <<"\n");
    }
}

void ImageSetHybridEvolver::
_computeEvolutionForEvents(std::list< HybridTimedSetType >& working_sets,
						   EnclosureListType& intermediate_sets,
						   const DiscreteLocation& location,
						   const std::set<DiscreteEvent>& blocking_events,
						   const EventListType& events,
						   const ActivationTimesType& activation_times,
						   const SetModelType& flow_set_model,
						   const TimeModelType& time_model,
						   const TimeModelType& blocking_time_model,
						   const Float& time_step,
						   bool ignore_activations,
						   Semantics semantics) const
{
    TimeModelType final_time_model=time_model+blocking_time_model*time_step;
    ARIADNE_LOG(2,"final_time_range="<<final_time_model.range()<<"\n");
    SetModelType evolved_set_model=this->getCalculusInterface(semantics).integration_step(flow_set_model,blocking_time_model);
    ARIADNE_LOG(2,"evolved_set_model.argument_size()="<<evolved_set_model.argument_size()<<"\n");
    ARIADNE_LOG(2,"evolved_set_range="<<evolved_set_model.range()<<"\n");
    // Compute evolution for blocking events
    for(std::set<DiscreteEvent>::const_iterator iter=blocking_events.begin(); iter!=blocking_events.end(); ++iter) {
        const DiscreteEvent event=*iter;
        if(event==finishing_event) {
            // TODO: Better estimate to use smaller blocking time
            intermediate_sets.adjoin(make_pair(location,evolved_set_model));
            working_sets.push_back(make_tuple(location,events,evolved_set_model,final_time_model));
        } else {
            EventKind kind = _sys->event_kind(location,event);
            bool is_transition = (kind == URGENT || kind == PERMISSIVE);
        	if(is_transition && !ignore_activations) {
        		intermediate_sets.adjoin(make_pair(location,evolved_set_model));
				SetModelType jump_set_model=apply(_sys->reset_function(location,event),evolved_set_model);
				DiscreteLocation jump_location=_sys->target(location,event);
				std::vector<DiscreteEvent> jump_events=events;
				jump_events.push_back(event);
				working_sets.push_back(make_tuple(jump_location,jump_events,jump_set_model,final_time_model));
        	}
        }
    }
    if (!ignore_activations) {
		// Compute evolution for non-blocking events
		for(std::map< DiscreteEvent, tuple<TimeModelType,TimeModelType> >::const_iterator
			iter=activation_times.begin(); iter!=activation_times.end(); ++iter) {
			const DiscreteEvent event=iter->first;
			const TimeModelType lower_active_time_model=iter->second.first;
			const TimeModelType upper_active_time_model=iter->second.second;
			ARIADNE_LOG(3,"Non blocking event "<<event<<":\n");
			ARIADNE_LOG(3,"  lower_active_time_model="<<lower_active_time_model.range()<<";\n");
			ARIADNE_LOG(3,"  upper_active_time_model="<<upper_active_time_model.range()<<";\n");
			SetModelType active_set_model=this->getCalculusInterface(semantics).reachability_step(flow_set_model,lower_active_time_model,upper_active_time_model);
			ARIADNE_LOG(3,"  active_set="<<active_set_model.range()<<";\n");
			SetModelType jump_set_model=apply(_sys->reset_function(location,event),active_set_model);
			ARIADNE_LOG(3,"  jump_set_model="<<active_set_model.range()<<";\n");
			const TimeModelType active_time_model = this->getCalculusInterface(semantics).reachability_time(time_model+lower_active_time_model*time_step,time_model+upper_active_time_model*time_step);
			ARIADNE_LOG(3,"  active_time_model="<<active_time_model.range()<<".\n");

			DiscreteLocation jump_location=_sys->target(location,event);
			std::vector<DiscreteEvent> jump_events=events;
			jump_events.push_back(event);
			working_sets.push_back(make_tuple(jump_location,jump_events,jump_set_model,active_time_model));
		}
    }
}


void
ImageSetHybridEvolver::
_compute_blocking_info(std::set<DiscreteEvent>& non_transverse_events,
				  std::set<DiscreteEvent>& blocking_events,
				  TimeModelType& blocking_time_model,
				  const TimeModelType& time_step_model,
				  const SetModelType& flow_set_model,
				  const std::map<DiscreteEvent,RealScalarFunction>& urgent_guards,
				  double SMALL_RELATIVE_TIME,
				  Semantics semantics) const
{
    // Compute event blocking times
    std::map<DiscreteEvent, TimeModelType> event_blocking_times;

    event_blocking_times[finishing_event]=time_step_model;
    compute_eventBlockingTimes_and_nonTransverseEvents(event_blocking_times,non_transverse_events,urgent_guards,flow_set_model,semantics);
    ARIADNE_LOG(2,"event_blocking_times="<<event_blocking_times<<"\n");

    std::map<DiscreteEvent,Interval> event_blocking_time_intervals;
    for(std::map<DiscreteEvent, TimeModelType>::const_iterator iter=event_blocking_times.begin();
        iter!=event_blocking_times.end(); ++iter) { event_blocking_time_intervals[iter->first]=iter->second.range(); }

    ARIADNE_LOG(2,"event_blocking_time_intervals="<<event_blocking_time_intervals<<"\n");
    ARIADNE_LOG(2,"non_transverse_events="<<non_transverse_events<<"\n\n");

    // Compute blocking events
    compute_blockingTime_and_relatedEvents(blocking_events,blocking_time_model,event_blocking_times);
    ARIADNE_LOG(2,"blocking_events="<<blocking_events<<"\n");
    ARIADNE_LOG(2,"blocking_time="<<blocking_time_model.range()<<"\n\n");

    // If multiple blocking events and flow time is large, move closer to blocking time
    if(blocking_events.size()!=1 && blocking_time_model.range().upper()>SMALL_RELATIVE_TIME) {
        blocking_events.clear();
        blocking_events.insert(finishing_event);
        if(blocking_time_model.range().lower()>0) {
            blocking_time_model-=blocking_time_model.error();
            blocking_time_model.set_error(0.0);
        }
        blocking_time_model*=(1-SMALL_RELATIVE_TIME);
    }
    ARIADNE_LOG(2,"(adjusted)blocking_events="<<blocking_events<<"\n");
    ARIADNE_LOG(2,"(adjusted)blocking_time="<<blocking_time_model.range()<<"\n\n");
}


void
ImageSetHybridEvolver::
_compute_activation_info(std::map<DiscreteEvent,RealScalarFunction>& activations,
						 ActivationTimesType& activation_times,
						 const std::set<DiscreteEvent>& non_transverse_events,
						 const SetModelType& flow_set_model,
						 const TimeModelType& blocking_time_model,
						 const std::map<DiscreteEvent,RealScalarFunction>& urgent_guards,
						 const std::map<DiscreteEvent,RealScalarFunction>& invariants,
						 const Semantics semantics) const
{

    // Treat non-transverse urgent events as non-urgent in upper semantics
    for(std::set<DiscreteEvent>::const_iterator iter=non_transverse_events.begin(); iter!=non_transverse_events.end(); ++iter) {
    	bool is_transition = (invariants.find(*iter) == invariants.end());
        if(is_transition)
            activations[*iter]=urgent_guards.find(*iter)->second;
    }

    this->compute_activationTimes(activation_times,activations,flow_set_model,blocking_time_model,semantics);

    // Display activation time ranges
    std::map< DiscreteEvent,tuple<Interval> > activation_time_intervals;
    for(ActivationTimesType::const_iterator iter=activation_times.begin(); iter!=activation_times.end(); ++iter)
        activation_time_intervals.insert(make_pair(iter->first,make_tuple(iter->second.second.range())));
    ARIADNE_LOG(2,"activation_time_intervals="<<activation_time_intervals<<"\n\n");
}


void
ImageSetHybridEvolver::
_compute_and_adjoin_reachableSet(EnclosureListType& reach_sets,
								SetModelType& reachable_set,
								const DiscreteLocation& location,
							    const SetModelType& flow_set_model,
							    const TimeModelType& zero_time_model,
		 	 	 	 	 	 	const TimeModelType& blocking_time_model,
		 	 	 	 	 	 	Semantics semantics) const
{
    ARIADNE_LOG(4,"flow_set_model="<<flow_set_model<<"\n");
    ARIADNE_LOG(4,"zero_time_model="<<zero_time_model<<"\n");
    ARIADNE_LOG(4,"blocking_time_model="<<blocking_time_model<<"\n");
    reachable_set=this->getCalculusInterface(semantics).reachability_step(flow_set_model,zero_time_model,blocking_time_model);
    reach_sets.adjoin(make_pair(location,reachable_set));

	ARIADNE_LOG(2,"reachable_set="<<reachable_set<<"\n");
    ARIADNE_LOG(2,"reachable_set.argument_size()="<<reachable_set.argument_size()<<"\n");
    ARIADNE_LOG(2,"reachable_set.range()="<<reachable_set.range()<<"\n");
}

void
ImageSetHybridEvolver::
_logEvolutionStepInitialState(const EventListType& previous_events,
							  const TimeModelType& time_model,
							  const DiscreteLocation& location,
							  const SetModelType& set_model,
							  const RealVectorFunction& dynamic,
							  const std::map<DiscreteEvent,RealScalarFunction>& invariants,
							  const std::map<DiscreteEvent,RealScalarFunction>& blocking_guards,
							  const std::map<DiscreteEvent,RealScalarFunction>& permissive_guards) const
{
    ARIADNE_LOG(2,"previous events = "<<previous_events<<" ");
    ARIADNE_LOG(2,"time_range = "<<time_model.range()<<" ");
    ARIADNE_LOG(2,"time_model generators = "<<time_model.argument_size()<<" ");
    ARIADNE_LOG(2,"location = "<<location<<" ");
    ARIADNE_LOG(2,"box = "<<set_model.range()<<" ");
    ARIADNE_LOG(2,"generators = "<<set_model.argument_size()<<" ");
    ARIADNE_LOG(2,"radius = "<<radius(set_model.range())<<"\n\n");

    ARIADNE_LOG(2,"dynamic = "<<dynamic<<"\n");
    ARIADNE_LOG(2,"invariants = "<<invariants<<"\n");
    ARIADNE_LOG(2,"blocking_guards = "<<blocking_guards<<"\n");
    ARIADNE_LOG(2,"permissive_guards = "<<permissive_guards<<"\n\n");
}



void
ImageSetHybridEvolver::
_evolution_add_initialSet(
		std::list< HybridTimedSetType >& working_sets,
		const EnclosureType& initial_set,
		Semantics semantics) const
{
    ARIADNE_LOG(6,"initial_set = "<<initial_set<<"\n");
    DiscreteLocation initial_location;
    ContinuousEnclosureType initial_continuous_set;
    make_lpair(initial_location,initial_continuous_set)=initial_set;
    ARIADNE_LOG(6,"initial_location = "<<initial_location<<"\n");
    SetModelType initial_set_model=this->getCalculusInterface(semantics).set_model(initial_continuous_set);

	// Check for non-zero maximum step size
	ARIADNE_ASSERT_MSG(this->_settings->hybrid_maximum_step_size[initial_location] > 0, "Error: the maximum step size for location " << initial_location.name() << " is zero.");
	// Check for match between the enclosure cell size and the set size
	ARIADNE_ASSERT_MSG(this->_settings->maximum_enclosure_cell.size() == initial_set_model.size(), "Error: mismatch between the maximum_enclosure_cell size and the set size.");

    ARIADNE_LOG(6,"initial_set_model = "<<initial_set_model<<"\n");
    TimeModelType initial_time_model=this->getCalculusInterface(semantics).time_model(0.0,Box(initial_set_model.argument_size()));
    ARIADNE_LOG(6,"initial_time_model = "<<initial_time_model<<"\n");
    TimedSetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
    ARIADNE_LOG(6,"initial_timed_set_model = "<<initial_timed_set_model<<"\n");
    working_sets.push_back(make_tuple(initial_location,EventListType(),initial_set_model,initial_time_model));
}

bool
ImageSetHybridEvolver::
_is_enclosure_too_large(const SetModelType& initial_set_model) const
{
	const Vector<Interval> initial_set_model_range = initial_set_model.range();

	for (uint i=0;i<initial_set_model_range.size();++i)
		if (initial_set_model_range[i].width() > this->_settings->maximum_enclosure_cell[i])
			return true;

	return false;
}

void
ImageSetHybridEvolver::
_add_subdivisions(std::list< HybridTimedSetType >& working_sets,
				  const array< TimedSetModelType >& subdivisions,
				  const DiscreteLocation& initial_location,
				  const EventListType& initial_events,
				  const uint dimension) const
{
    ARIADNE_LOG(3,"subdivisions.size()="<<subdivisions.size()<<"\n");
    for(uint i=0; i!=subdivisions.size(); ++i) {
        TimedSetModelType const& subdivided_timed_set_model=subdivisions[i];
        ARIADNE_LOG(3,"subdivided_timed_set_model.range()="<<subdivided_timed_set_model.range()<<"\n");
        SetModelType subdivided_set_model=Vector<TaylorModel>(project(subdivided_timed_set_model.models(),range(0,dimension)));
        TimeModelType subdivided_time_model=subdivided_timed_set_model[dimension];
        ARIADNE_LOG(3,"subdivided_set_model.range()="<<subdivided_set_model.range()<<"\n");
        ARIADNE_LOG(3,"subdivided_set_model.radius()*10000="<<radius(subdivided_set_model.range())*10000<<"\n");
        ARIADNE_LOG(3,"subdivided_time_model.range()="<<subdivided_time_model.range()<<"\n");
        working_sets.push_back(make_tuple(initial_location,initial_events,subdivided_set_model,subdivided_time_model));
    }
}

void
ImageSetHybridEvolver::
_add_models_subdivisions_autoselect(std::list< HybridTimedSetType >& working_sets,
		  	  	  	  	  	  	  	const SetModelType& initial_set_model,
		  	  	  	  	  	  	  	const TimeModelType& initial_time_model,
		  	  	  	  	  	  	  	const DiscreteLocation& initial_location,
		  	  	  	  	  	  	  	const EventListType& initial_events,
		  	  	  	  	  	  	  	Semantics semantics) const
{
    uint nd=initial_set_model.dimension();
    SetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
    array< TimedSetModelType > subdivisions=this->getCalculusInterface(semantics).subdivide(initial_timed_set_model);
    _add_subdivisions(working_sets,subdivisions,initial_location,initial_events,nd);
}


void
ImageSetHybridEvolver::
_add_models_subdivisions_time(std::list< HybridTimedSetType >& working_sets,
		  	  	  	  	  	  const SetModelType& initial_set_model,
		  	  	  	  	  	  const TimeModelType& initial_time_model,
		  	  	  	  	  	  const DiscreteLocation& initial_location,
		  	  	  	  	  	  const EventListType& initial_events,
		  	  	  	  	  	  Semantics semantics) const
{
    uint nd=initial_set_model.dimension();
    SetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
    array< TimedSetModelType > subdivisions=this->getCalculusInterface(semantics).subdivide(initial_timed_set_model,nd);
    _add_subdivisions(working_sets,subdivisions,initial_location,initial_events,nd);
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

	// Gets the size of the continuous space
	const uint css = hgrid.locations_begin()->second.lengths().size();

	// Initializes the result
	Vector<Float> result(css);

	// For each location and dimension of the space
	for (HybridGrid::locations_const_iterator hg_it = hgrid.locations_begin(); hg_it != hgrid.locations_end(); hg_it++) {
		ARIADNE_ASSERT_MSG(hg_it->second.lengths().size() == css,
				"The maximum enclosure cell is obtained under the assumption that the continuous state space dimension is the same for all locations.");
		for (uint i=0;i<css;i++)
			if (hg_it->second.lengths()[i] > result[i])
				result[i] = hg_it->second.lengths()[i];
	}

	return RATIO*result/(1<<maximum_grid_depth);
}


std::map<DiscreteLocation,Float>
getHybridMaximumStepSize(
		const HybridFloatVector& hmad,
		const HybridGrid& hgrid,
		int maximum_grid_depth,
		Semantics semantics)
{
	// We choose a coefficient for upper semantics such that an enclosure at maximum size is able to cross
	// urgent transitions in one step. For lower semantics we prefer to have a finer result.
	Float coefficient = (semantics == UPPER_SEMANTICS ? 2.0 : 1.0);

	// Initialize the hybrid maximum step size
	std::map<DiscreteLocation,Float> hmss;

	// For each couple DiscreteLocation,Vector<Float>
	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
	{
		const uint dim = hfv_it->second.size();
		// Initializes the maximum step size
		Float mss = 0.0;
		// For each dimension of the space, if the derivative is not zero,
		// evaluates the ratio between the minimum cell length and the derivative itself
		for (uint i=0;i<dim;i++)
			if (hfv_it->second[i] > 0)
				mss = max(mss,hgrid[hfv_it->first].lengths()[i]/(1 << maximum_grid_depth)/hfv_it->second[i]);

		// Inserts the value (twice the value since the maximum enclosure is set as ~2 the grid cell)
		hmss.insert(std::pair<DiscreteLocation,Float>(hfv_it->first,coefficient*mss));
	}

	return hmss;
}


}  // namespace Ariadne

