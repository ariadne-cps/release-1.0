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

#include "logging.h"

#include "hybrid_time.h"
#include "hybrid_automaton.h"
#include "hybrid_evolver-image.h"
#include "interruptible.h"

#include "exceptions.h"

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
    : EvolverBase<SystemType,LocalisedTaylorSet>(system)
    , _settings(new SettingsType(system))
    , _toolbox(new TaylorCalculus())
{
    this->charcode = "e";
}


void
ImageSetHybridEvolver::
tune_settings(
		const HybridBoxes& domain,
		const HybridGrid& grid,
		AccuracyType accuracy)
{
    HybridFloatVector hmad = getHybridMidpointAbsoluteDerivatives(*_sys,domain);

    ARIADNE_LOG(1, "Tuning settings for evolution...");
	this->_settings->_reference_enclosure_widths = getMinimumGridCellWidths(grid,accuracy);
	ARIADNE_LOG(2, "Reference enclosure widths: " << this->_settings->_reference_enclosure_widths);
	this->_settings->set_maximum_step_size(getHybridMaximumStepSize(hmad,grid,accuracy));
	ARIADNE_LOG(2, "Fixed maximum step size: " << this->_settings->maximum_step_size());
}

ContinuousStepResult
ImageSetHybridEvolver::
compute_integration_step_result(const SetModelType& starting_set,
               const DiscreteLocation& location,
               ContinuousEvolutionDirection direction,
               Float step) const {

    RealVectorFunction dynamic = get_directed_dynamic(_sys->dynamic_function(location),direction);

    const int MAXIMUM_BOUNDS_DIAMETER_FACTOR = 8;
    Float maximum_bounds_diameter=max(this->_settings->_reference_enclosure_widths.find(location)->second)*
            MAXIMUM_BOUNDS_DIAMETER_FACTOR*this->_settings->maximum_enclosure_widths_ratio();

    SetModelType flow_model;
    Float effective_step;
    make_lpair(flow_model,effective_step) = compute_flow_and_effective_step(step, dynamic, maximum_bounds_diameter, starting_set);
    TaylorSet finishing_set = partial_evaluate(flow_model.models(),starting_set.argument_size(),1.0);

    return ContinuousStepResult(effective_step,flow_model,finishing_set);
}

ContinuousStepResult
ImageSetHybridEvolver::
_adaptive_step_and_flow(const SetModelType& starting_set,
                        const Float& previous_step,
                        const DiscreteLocation& location,
                        ContinuousEvolutionDirection direction,
                        const Float& remaining_time,
                        const Float& maximum_step,
                        const SetModelType& maximum_flow_model,
                        const SetModelType& maximum_finishing_model) const {

    uint refinement_radius = 3;
    Float improvement_percentage = 0.1;

    Float dim = starting_set.dimension();

    Vector<Float> global_target_widths_ratio_score_terms(dim);
    Vector<Float> global_target_scaled_error_rates(dim);
    Vector<Float> final_widths = this->_settings->_reference_enclosure_widths.find(location)->second;
    for (int i = 0; i < dim; ++i) {
        global_target_widths_ratio_score_terms[i] = (starting_set.widths()[i]/final_widths[i]);
        global_target_scaled_error_rates[i] = (final_widths[i]-starting_set.widths()[i])/final_widths[i]/remaining_time;
    }

    Float resuming_step = min(previous_step,maximum_step);

    std::list<Float> steps;

    if (starting_set.radius() == 0) {
        steps.push_back(resuming_step / std::pow(2,refinement_radius*4));
    } else {
        Float current = resuming_step;
        for (int i=0; i < refinement_radius; ++i) {
            current *= 2;
            if (current < maximum_step)
                steps.push_front(current);
            else {
                break;
            }
        }
        current = resuming_step*3/2;
        if (current < maximum_step)
            steps.push_back(current);

        steps.push_back(resuming_step);
        steps.push_back(resuming_step*3/4);

        current = resuming_step;
        for (int i=0; i < refinement_radius; ++i) {
            steps.push_back(current /= 2);
        }
    }

    std::vector<tuple<ContinuousStepResult,Float,Float> > candidates;

    for (std::list<Float>::const_iterator it = steps.begin(); it != steps.end(); ++it) {

        ContinuousStepResult integration_step_result = (*it == maximum_step ?
                                                        ContinuousStepResult(*it,maximum_flow_model,maximum_finishing_model) :
                                                        compute_integration_step_result(starting_set,location,direction,*it));

        Float exponent = integration_step_result.used_step()/remaining_time;
        Vector<Float> target_width_ratios(dim);
        for (uint i = 0; i < dim; ++i) {
            target_width_ratios[i] = std::pow((Float)global_target_widths_ratio_score_terms[i],exponent);
        }

        Vector<Float> target_scaled_errors(dim);
        for (uint i = 0; i < dim; ++i) {
            if (starting_set.widths()[i] > 0)
                target_scaled_errors[i] = (starting_set.widths()[i]/target_width_ratios[i]-starting_set.widths()[i])/final_widths[i];
        }
        Float target_scaled_error_score = sum(target_scaled_errors);

        Vector<Float> actual_scaled_errors(dim);
        for (uint i = 0; i < dim; ++i) {
            if (starting_set.widths()[i] > 0)
                actual_scaled_errors[i] = (integration_step_result.finishing_set_model().widths()[i]-starting_set.widths()[i])/final_widths[i];
        }
        Float actual_scaled_error_score = sum(actual_scaled_errors);

        candidates.push_back(make_tuple(integration_step_result,target_scaled_error_score,actual_scaled_error_score));
    }

    std::vector<tuple<ContinuousStepResult,Float,Float> >::const_iterator it = candidates.begin();
    tuple<ContinuousStepResult,Float,Float> winner = *it;
    bool target_hit = false;
    for ( ; it != candidates.end(); ++it) {
        Float target_score = it->second;
        Float actual_score = it->third;
        Float current_step = it->first.used_step();
        Float winner_step = winner.first.used_step();
        Float current_relative_score = (actual_score-target_score)/abs(target_score);
        Float winner_relative_score = (winner.third-winner.second)/abs(winner.second);
        Float improvement = (winner_relative_score - current_relative_score)/abs(winner_relative_score);

        // If we improve on the target score for the first time, we set the winner
        if (!target_hit && current_relative_score<0) {
            target_hit = true;
            winner = *it;
        } else {
            if (current_relative_score < winner_relative_score) {
                if (improvement > winner_step/current_step * improvement_percentage)
                    winner = *it;
            }
        }

        cout << "Step " << current_step <<
                ", tgt $ " << target_score <<
                ", act $ " << actual_score <<
                ", rel $ " << current_relative_score <<
                " (" << improvement*100 << "% improv. for " << winner_step/current_step << "x finer step)" <<
                (winner.first.used_step() == current_step ? " <" : "") <<
        endl;

    }

    return winner.first;
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
    _reset_start_time();

    ARIADNE_LOG(1,"Computing evolution up to "<<maximum_hybrid_time.continuous_time()<<" time units and "<<maximum_hybrid_time.discrete_time()<<" steps.");

    std::list<EvolutionData> working_sets;
    _evolution_add_initialSet(working_sets,initial_set,semantics);
    std::map<uint,Vector<Float> > initial_indexed_set_models_widths = _indexed_set_models_widths(working_sets);

    // While there exists a working set, process it
    while(!working_sets.empty()) {

		if (this->_settings->maximum_number_of_working_sets() > 0 && working_sets.size() > this->_settings->maximum_number_of_working_sets())
			throw WorkingSetTooLargeException("too large");

		// Get the least recent working set, pop it and update the corresponding size
		EvolutionData current_set = working_sets.front();
		working_sets.pop_front();

		// Get the members of the current set
		uint set_index = current_set.set_index();
        DiscreteLocation loc=current_set.location();
        StepType previous_step=current_set.previous_step();
        EventListType previous_events=current_set.previous_events();
		SetModelType set_model=current_set.set_model();
		TimeModelType time_model=current_set.time_model();

		// Compute continuous evolution

		Vector<Float> reference_enclosure_widths = this->_settings->_reference_enclosure_widths.find(loc)->second;

		bool isEnclosureTooLarge = _is_enclosure_too_large(loc,set_model,initial_indexed_set_models_widths[set_index]);

		if(time_model.range().lower()>=maximum_hybrid_time.continuous_time() ||
		    previous_events.size()>=uint(maximum_hybrid_time.discrete_time())) {
            ARIADNE_LOG(2,"Final time reached, adjoining result to final sets.");
            final_sets.adjoin(loc,_toolbox->enclosure(set_model));
        } else if (semantics == UPPER_SEMANTICS && this->_settings->enable_subdivisions() && isEnclosureTooLarge) {
            ARIADNE_LOG(2,"Computed set range " << set_model.range() << " widths larger than allowed, subdividing.");
            _add_models_subdivisions_autoselect(working_sets,set_index,set_model,time_model,loc,previous_step,previous_events,semantics);
        } else if((semantics == LOWER_SEMANTICS || !this->_settings->enable_subdivisions()) &&
                  this->_settings->enable_premature_termination_on_enclosure_size() && isEnclosureTooLarge) {
            ARIADNE_LOG(2,"Terminating evolution at time " << time_model.value()
                        << " and set " << set_model.centre() << " (widths " << set_model.widths() << ") due to maximum enclosure bounds being exceeded.");
            if(semantics == UPPER_SEMANTICS)
                final_sets.adjoin(loc,_toolbox->enclosure(set_model));
        } else {
            this->_evolution_step(working_sets,reach_sets,intermediate_sets,current_set,maximum_hybrid_time,
            		ignore_activations,direction,semantics);
        }

		_check_timeout();
    }
}


void
ImageSetHybridEvolver::
_evolution_step(std::list<EvolutionData>& working_sets,
                EnclosureListType& reach_sets,
                EnclosureListType& intermediate_sets,
                const EvolutionData& current_set,
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
    DiscreteLocation location = current_set.location();
    StepType previous_step = current_set.previous_step();
    EventListType events_history = current_set.previous_events();
    SetModelType set_model = current_set.set_model();
    TimeModelType time_model = current_set.time_model();
    uint set_index = current_set.set_index();

    _log_step_summary(working_sets,reach_sets,events_history,time_model,set_model,location,previous_step);

    Float remaining_time = maximum_hybrid_time.continuous_time() - time_model.range().lower();

    if (_settings->enable_reconditioning()) {

        bool has_reconditioned = false;

    	Vector<Float> error_thresholds(set_model.dimension());
    	Vector<Float> reference_enclosure_widths = _settings->reference_enclosure_widths().at(location);
    	for (uint i = 0; i < reference_enclosure_widths.size(); ++i)
    		error_thresholds[i] = reference_enclosure_widths[i]/4;

    	set_model.uniform_error_recondition(error_thresholds);

    	int argument_difference = set_model.argument_size() - time_model.argument_size();
    	if (argument_difference > 0) {
    	    has_reconditioned = true;
    		time_model = embed(time_model,argument_difference);
    	}

    	Array<uint> discarded_parameters = set_model.kuhn_recondition();
    	if (!discarded_parameters.empty()) {
    	    has_reconditioned = true;
    		time_model = recondition(time_model,discarded_parameters,set_model.dimension(),set_model.dimension());
    	}
    }

    // Gets the maximum acceptable step along with the related flow model
    Float proposed_maximum_step = min(remaining_time, _settings->maximum_step_size().at(location));
    ContinuousStepResult maximum_step_and_flow = compute_integration_step_result(set_model,location,direction,proposed_maximum_step);
    Float effective_maximum_step = maximum_step_and_flow.used_step();
    SetModelType flow_set_model = maximum_step_and_flow.flow_set_model();
    SetModelType finishing_set_model = maximum_step_and_flow.finishing_set_model();

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
    			break;
    	}
    }

    // Check to make sure dimensions are correct

    const RealVectorFunction dynamic=get_directed_dynamic(_sys->dynamic_function(location),direction);

    ARIADNE_ASSERT(set_model.argument_size()==time_model.argument_size());
    ARIADNE_ASSERT_MSG(set_model.result_size()==_sys->dimension(location),"set_model="<<set_model<<", location="<<location.name());

    _logEvolutionStepInitialState(events_history,time_model,location,set_model,dynamic,invariants,urgent_guards,permissive_guards);

    ARIADNE_LOG(2,"effective_time_step = "<<effective_maximum_step);
    ARIADNE_LOG(2,"flow_range = "<<flow_set_model.range());
    ARIADNE_LOG(2,"finishing_set_range = "<<maximum_step_and_flow.finishing_set_model().range());

    // Set special events and times; note that the time step is scaled to [0,1]
    TimeModelType zero_time_model = _toolbox->time_model(0.0,Box(time_model.argument_size()));
    TimeModelType time_step_model = _toolbox->time_model(1.0,Box(time_model.argument_size()));

    std::set<DiscreteEvent> blocking_events;
    TimeModelType blocking_time_model;
    std::set<DiscreteEvent> non_transverse_events;
    _compute_blocking_info(non_transverse_events,blocking_events,blocking_time_model,
    				  time_step_model,flow_set_model,urgent_guards,SMALL_RELATIVE_TIME,semantics);

    Float event_reduced_maximum_step = effective_maximum_step;


    // Only if the proposed step if at least twice larger, and if the blocking time is positive (thus excluding sets outside the invariants)
    if (blocking_time_model.range().upper() < 0.5 && blocking_time_model.range().upper() > 0) {

        event_reduced_maximum_step = 1.5*(blocking_time_model * effective_maximum_step).range().upper();
        ContinuousStepResult new_flow_and_step = compute_integration_step_result(set_model,location,direction,event_reduced_maximum_step);
        flow_set_model = new_flow_and_step.flow_set_model();
        event_reduced_maximum_step = new_flow_and_step.used_step();
    }

    Float step = event_reduced_maximum_step;

    if (_settings->enable_error_rate_enforcement()) {
        ContinuousStepResult continuous_step_result = _adaptive_step_and_flow(set_model, previous_step, location, direction, remaining_time, event_reduced_maximum_step, flow_set_model,finishing_set_model);
        flow_set_model = continuous_step_result.flow_set_model();
        step = continuous_step_result.used_step();
    }

    if (step < effective_maximum_step) {
        blocking_events.clear();
        non_transverse_events.clear();
        _compute_blocking_info(non_transverse_events,blocking_events,blocking_time_model,
                          time_step_model,flow_set_model,urgent_guards,SMALL_RELATIVE_TIME,semantics);
    }

    if (is_enclosure_to_be_discarded(set_model,urgent_guards,dynamic,semantics))
        return;

    ActivationTimesType activation_times;
    if (!ignore_activations) {
    	_compute_activation_info(permissive_guards,activation_times,non_transverse_events,
    			flow_set_model,blocking_time_model,urgent_guards,invariants,semantics);
    }

    SetModelType reachable_set;
    _compute_and_adjoin_reachableSet(reach_sets,reachable_set,location,flow_set_model,zero_time_model,blocking_time_model,semantics);

    if(semantics!=LOWER_SEMANTICS || blocking_events.size()==1)
    	_computeEvolutionForEvents(working_sets,intermediate_sets,set_index,location,blocking_events,events_history,
    								activation_times,flow_set_model,time_model,blocking_time_model,step,ignore_activations,semantics);
}


ImageSetHybridEvolver::TimeModelType ImageSetHybridEvolver::
crossing_time(
		VectorFunction guard,
		const FlowSetModelType& flow_set_model,
		Semantics semantics) const
{
    try {
        TimeModelType crossing_time_model=_toolbox->scaled_crossing_time(guard,flow_set_model);
        return crossing_time_model;
    }
    catch(DegenerateCrossingException e) {
        BoxType space_domain=project(flow_set_model.domain(),range(0,flow_set_model.argument_size()-1));
        Interval touching_time_interval=_toolbox->scaled_touching_time_interval(guard,flow_set_model);
        TimeModelType touching_time_model=_toolbox->time_model(touching_time_interval,space_domain);
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
	ARIADNE_LOG(2,"Computing initially active events...");
    tribool blocking_event_initially_active=false;
    for(std::map<DiscreteEvent,RealScalarFunction>::const_iterator iter=urgent_guards.begin(); iter!=urgent_guards.end(); ++iter) {
        RealVectorFunction activation(1,iter->second);

        Vector<Interval> infinite_box(initial_set.size(),Interval(-std::numeric_limits<double>::max(),std::numeric_limits<double>::max()));

        tribool everywhere_active = _toolbox->active(iter->second,infinite_box);
        if (definitely(everywhere_active)) {
        	ARIADNE_LOG(3,"Ignoring urgent guard '" << iter->first.name() << "' with activation true...");
        } else {
			ARIADNE_LOG(3,"Evaluating urgent guard '" << iter->first.name() << "' with activation " << activation << "...");
			tribool initially_active=_toolbox->active(activation,initial_set);
			if(possibly(initially_active)) {
				ARIADNE_LOG(3,"Possibly active.");
				initially_active_events.insert(std::make_pair(iter->first,initially_active));
				blocking_event_initially_active = blocking_event_initially_active || initially_active;
			} else {
				ARIADNE_LOG(3,"Inactive.");
			}
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
	ARIADNE_LOG(3,"Checking nonnegative crossings...");
    for(std::map<DiscreteEvent,RealScalarFunction>::const_iterator iter=urgent_guards.begin(); iter!=urgent_guards.end(); ++iter) {
        RealVectorFunction activation(1,iter->second);
        ARIADNE_LOG(4,"Guard: " << activation);
        tribool is_active = _toolbox->active(activation,set_bounds);
        tribool is_positively_crossing = positively_crossing(set_bounds,dynamic,activation[0]);
        ARIADNE_LOG(4,"Active: " << is_active << "; positively crossing: " << is_positively_crossing);
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
	ARIADNE_LOG(2,"Checking whether the enclosure is to be discarded...")
	std::map<DiscreteEvent,tribool> initially_active_events;
	this->compute_initially_active_events(initially_active_events, urgent_guards, enclosure, semantics);
	ARIADNE_LOG(3,"Initially_active_events = "<<initially_active_events);

	tribool has_any_initially_active_blocking_event = initially_active_events[blocking_event];

	// Test for initially active events, and process these as required
	if(definitely(has_any_initially_active_blocking_event)) {
		ARIADNE_LOG(3,"An invariant is definitely initially active: discarding the set.");
		result = true;
	} else if(possibly(has_any_initially_active_blocking_event) && semantics==LOWER_SEMANTICS) {
		ARIADNE_LOG(3,"A blocking event is possibly active: checking whether there is a possibly positive crossing.");
		bool has_nonneg_crossing = has_nonnegative_crossing(urgent_guards,dynamic,enclosure.bounding_box(),semantics);
		if (has_nonneg_crossing) {
			ARIADNE_LOG(3,"Terminating lower evolution due to possibly initially active invariant with nonnegative crossing.");
			result = true;
		}
	}

	return result;
}

// Compute the flow, parameterising space with the set parameters
void ImageSetHybridEvolver::
compute_flow_model(
		const DiscreteLocation& loc,
		FlowSetModelType& flow_set_model,
		BoxType& flow_bounds, Float& time_step,
        VectorFunction dynamic,
        const SetModelType& starting_set_model,
        const TimeModelType& starting_time_model,
        Float finishing_time) const
{
    ARIADNE_LOG(2,"Computing flow model...");
    const int MAXIMUM_BOUNDS_DIAMETER_FACTOR = 8;
    float remaining_time = finishing_time - starting_time_model.range().lower();
    const Float maximum_step_size=min(time_step, remaining_time);

    const Float maximum_bounds_diameter=max(this->_settings->_reference_enclosure_widths.find(loc)->second)*
            MAXIMUM_BOUNDS_DIAMETER_FACTOR*this->_settings->maximum_enclosure_widths_ratio();

    BoxType starting_set_bounding_box=starting_set_model.range();
    ARIADNE_LOG(3,"starting_set_bounding_box="<<starting_set_bounding_box);
    make_lpair(time_step,flow_bounds)=_toolbox->flow_bounds(dynamic,starting_set_bounding_box,maximum_step_size,maximum_bounds_diameter);
    // Compute the flow model
    ARIADNE_LOG(3,"time_step="<<time_step);
    ARIADNE_LOG(3,"flow_bounds="<<flow_bounds);
    FlowModelType flow_model=_toolbox->flow_model(dynamic,starting_set_bounding_box,time_step,flow_bounds);
    ARIADNE_LOG(3,"flow_model="<<flow_model);
    ScalarTaylorFunction identity_time_expression=ScalarTaylorFunction::variable(BoxType(1u,Interval(-time_step,+time_step)),0u);
    flow_set_model=unchecked_apply(flow_model,combine(starting_set_model.models(),identity_time_expression.model()));
}

std::pair<ImageSetHybridEvolver::FlowSetModelType,Float> ImageSetHybridEvolver::
compute_flow_and_effective_step(
        Float& time_step,
        RealVectorFunction dynamic,
        Float& maximum_bounds_diameter,
        const SetModelType& starting_set_model) const
{
    Box starting_set_bounding_box=starting_set_model.range();
    Box flow_bounds;
    Float effective_time_step;
    make_lpair(effective_time_step,flow_bounds)=_toolbox->flow_bounds(dynamic,starting_set_bounding_box,time_step,maximum_bounds_diameter);
    TaylorCalculus::FlowModelType flow_model= _toolbox->flow_model(dynamic,starting_set_bounding_box,effective_time_step,flow_bounds);
    ScalarTaylorFunction identity_time_expression=ScalarTaylorFunction::variable(Box(1u,Interval(-effective_time_step,+effective_time_step)),0u);
    return std::make_pair(unchecked_apply(flow_model,combine(starting_set_model.models(),identity_time_expression.model())),effective_time_step);
}



void ImageSetHybridEvolver::
compute_eventBlockingTimes_and_nonTransverseEvents(
		std::map<DiscreteEvent,TimeModelType>& event_blocking_times,
        std::set<DiscreteEvent>& non_transverse_events,
        const std::map<DiscreteEvent,RealScalarFunction>& urgent_guards,
        const FlowSetModelType& flow_set_model,
        Semantics semantics) const
{
    ARIADNE_LOG(2,"Computing blocking events...");

    uint dimension=flow_set_model.result_size();
    const double SMALL_RELATIVE_TIME = 1./16;
    FlowSetModelType positive_flow_set_model(split(flow_set_model.models(),dimension,1));

    for(std::map<DiscreteEvent,RealScalarFunction>::const_iterator guard_iter=urgent_guards.begin();
        guard_iter!=urgent_guards.end(); ++guard_iter)
    {
        const DiscreteEvent event=guard_iter->first;
        const RealVectorFunction guard(1,guard_iter->second);
        tribool active = _toolbox->active(guard,positive_flow_set_model);
        if(possibly(active)) {
            ARIADNE_LOG(3,"Event "<<event<<" possibly active.");
            TimeModelType crossing_time_model;
            Interval normal_derivative;
            try {
                crossing_time_model=_toolbox->scaled_crossing_time(guard,flow_set_model);
                normal_derivative=this->normal_derivative(guard,flow_set_model,crossing_time_model);
                assert(normal_derivative.lower()>0 || normal_derivative.upper()<0);
                if(normal_derivative.lower()>0) {
                    ARIADNE_LOG(3,"Event "<<event<<" inserted into blocking times.");
                    event_blocking_times[event]=crossing_time_model;
                }
            }
            catch(DegenerateCrossingException& e) {
                ARIADNE_LOG(3,"Degenerate Crossing exception catched.");
                BoxType space_domain=project(flow_set_model.domain(),range(0,flow_set_model.argument_size()-1));
                Interval touching_time_interval=_toolbox->scaled_touching_time_interval(guard,flow_set_model);
                TimeModelType touching_time_model=_toolbox->time_model(touching_time_interval,space_domain);
                // Use 1.0 as upper bound above since flow set model has time interval normalised to [-1,+1]
                ARIADNE_LOG(3,"touching_time_interval="<<touching_time_interval);
                if(touching_time_interval.upper()>=0 && touching_time_interval.lower()<=1.0) {
                    SetModelType finishing_set_model=partial_evaluate(flow_set_model.models(),dimension,1.0);
                    tribool finishing_set_active=_toolbox->active(guard,finishing_set_model);
                    if(definitely(finishing_set_active)) {
                        ARIADNE_LOG(3,"Event is definitely finally active, inserting it into blocking times.");
                        event_blocking_times[event]=touching_time_model;
                    } else if(possibly(finishing_set_active)) {
                        ARIADNE_LOG(3,"Event is possibly finally active.");
                        if(touching_time_interval.lower()>SMALL_RELATIVE_TIME) {
                            ARIADNE_LOG(3,"lower touching time is greater than zero, inserting event into blocking times.");
                            TaylorModel lower_touching_time_model=_toolbox->time_model(touching_time_interval.lower(),space_domain);
                            event_blocking_times[finishing_event]=lower_touching_time_model;
                        } else {
                            ARIADNE_LOG(3,"DANGER: we can't determine whether the crossing is completely finished or not..");
                            // FIXME: Here we are stuck, we can't determine whether the crossing is completely finished or not.
                            // Just put this in as a blocking event and hope for the best...
                            // event_blocking_times[event]=touching_time_model;
                            // DAVIDE: setting this as a blocking event is not correct,
                            //         because it causes continuous evolution to stop.
                            //         I think it is better to put it into non-transverse-events.
                            non_transverse_events.insert(event);
                        }
                    } else {
                        ARIADNE_LOG(3,"After the flow step, the event is not active again, so the crossing was tangential.");
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
        tribool active=_toolbox->active(activation,flow_set_model);

        if(definitely(active)) {
            // The event is enabled over the entire time interval
            ARIADNE_LOG(3,"Event is enabled over the entire time interval.");
            activation_times.insert(make_pair(event,make_tuple(zero_time_model,blocking_time_model)));
        } else if(possibly(active)) {
            ARIADNE_LOG(3,"Event is possibly enabled.");
            // Compute whether the event is enabled at the beginning and end of the time interval
            tribool initially_active=_toolbox->active(activation,initial_set_model);
            tribool finally_active=_toolbox->active(activation,final_set_model);

            TimeModelType crossing_time_model=this->crossing_time(activation,flow_set_model,semantics);

            TimeModelType lower_crossing_time_model=crossing_time_model-crossing_time_model.error();
            TimeModelType upper_crossing_time_model=crossing_time_model+crossing_time_model.error();
            lower_crossing_time_model.set_error(0);
            upper_crossing_time_model.set_error(0);

            TimeModelType lower_active_time_model, upper_active_time_model;

            // Determine lower activation time
            if(definitely(not(initially_active))) {
                ARIADNE_LOG(3,"Event is definitely not initially active.");
                switch(semantics) {
                    case UPPER_SEMANTICS: lower_active_time_model=lower_crossing_time_model; break;
                    case LOWER_SEMANTICS: lower_active_time_model=upper_crossing_time_model; break;
                }
            } else if(definitely(initially_active)) {
                ARIADNE_LOG(3,"Event is definitely initially active.");
                lower_active_time_model=zero_time_model;
            } else {
                ARIADNE_LOG(3,"Event is undeterminately initially active.");
                switch(semantics) {
                    case UPPER_SEMANTICS: lower_active_time_model=zero_time_model; break;
                    case LOWER_SEMANTICS: lower_active_time_model=upper_crossing_time_model; break;
                }
            }

            // Compute upper activation time
            if(definitely(not(finally_active))) {
                ARIADNE_LOG(3,"Event is definitely not finally active.");
                switch(semantics) {
                    case UPPER_SEMANTICS: upper_active_time_model=upper_crossing_time_model; break;
                    case LOWER_SEMANTICS: upper_active_time_model=lower_crossing_time_model; break;
                }
            } else if(definitely(finally_active)) {
                ARIADNE_LOG(3,"Event is definitely finally active.");
                upper_active_time_model=blocking_time_model;
            } else {
                ARIADNE_LOG(3,"Event is undeterminately finally active.");
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
_log_step_summary(const std::list<EvolutionData>& working_sets,
					 const EnclosureListType& reach_sets,
					 const EventListType& initial_events,
					 const TimeModelType& initial_time_model,
					 const SetModelType& initial_set_model,
					 const DiscreteLocation& initial_location,
					 const Float& time_step) const
{
        ARIADNE_LOG(1,"#w="<<std::setw(4)<<working_sets.size()
                    <<"#r="<<std::setw(4)<<std::left<<reach_sets.size()
                    <<" s="<<std::setw(3)<<std::left<<initial_events.size()
                    <<" ps="<<std::scientific<<std::setw(5)<<std::left<< time_step << std::fixed
                    <<" t="<<std::fixed<<initial_time_model.value()
                    <<" r="<<std::scientific<<std::setw(7)<<initial_set_model.radius()<<std::fixed
                    <<" l="<<std::setw(3)<<std::left<<initial_location
                    <<" c="<<initial_set_model.centre()
                    <<" as="<<initial_set_model.argument_size()
					<<" w="<<initial_set_model.widths()
                    <<" e="<<initial_events);
}

void ImageSetHybridEvolver::
_computeEvolutionForEvents(std::list<EvolutionData>& working_sets,
						   EnclosureListType& intermediate_sets,
						   const uint& set_index,
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
    ARIADNE_LOG(2,"final_time_range="<<final_time_model.range());
    SetModelType evolved_set_model=_toolbox->integration_step(flow_set_model,blocking_time_model);
    ARIADNE_LOG(2,"evolved_set_model.argument_size()="<<evolved_set_model.argument_size());
    ARIADNE_LOG(2,"evolved_set_range="<<evolved_set_model.range());
    // Compute evolution for blocking events
    for(std::set<DiscreteEvent>::const_iterator iter=blocking_events.begin(); iter!=blocking_events.end(); ++iter) {
        const DiscreteEvent event=*iter;
        if(event==finishing_event) {
            // TODO: Better estimate to use smaller blocking time
            intermediate_sets.adjoin(make_pair(location,evolved_set_model));
            working_sets.push_back(EvolutionData(set_index,location,time_step,events,evolved_set_model,final_time_model));
        } else {
            EventKind kind = _sys->event_kind(location,event);
            bool is_transition = (kind == URGENT || kind == PERMISSIVE);
        	if(is_transition && !ignore_activations) {
        		intermediate_sets.adjoin(make_pair(location,evolved_set_model));
				SetModelType jump_set_model=apply(_sys->reset_function(location,event),evolved_set_model);
				DiscreteLocation jump_location=_sys->target(location,event);
				std::vector<DiscreteEvent> jump_events=events;
				jump_events.push_back(event);
				working_sets.push_back(EvolutionData(set_index,jump_location,std::numeric_limits<Float>::infinity(),jump_events,jump_set_model,final_time_model));
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
			ARIADNE_LOG(3,"Non blocking event "<<event<<":");
			ARIADNE_LOG(3,"lower_active_time_model="<<lower_active_time_model.range());
			ARIADNE_LOG(3,"upper_active_time_model="<<upper_active_time_model.range());
			SetModelType active_set_model=_toolbox->reachability_step(flow_set_model,lower_active_time_model,upper_active_time_model);
			ARIADNE_LOG(3,"active_set="<<active_set_model.range());
			SetModelType jump_set_model=apply(_sys->reset_function(location,event),active_set_model);
			ARIADNE_LOG(3,"jump_set_model="<<active_set_model.range());
			const TimeModelType active_time_model = _toolbox->reachability_time(time_model+lower_active_time_model*time_step,time_model+upper_active_time_model*time_step);
			ARIADNE_LOG(3,"active_time_model="<<active_time_model.range());

			DiscreteLocation jump_location=_sys->target(location,event);
			std::vector<DiscreteEvent> jump_events=events;
			jump_events.push_back(event);
			working_sets.push_back(EvolutionData(set_index,jump_location,std::numeric_limits<Float>::infinity(),jump_events,jump_set_model,active_time_model));
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
    ARIADNE_LOG(2,"event_blocking_times="<<event_blocking_times);

    std::map<DiscreteEvent,Interval> event_blocking_time_intervals;
    for(std::map<DiscreteEvent, TimeModelType>::const_iterator iter=event_blocking_times.begin();
        iter!=event_blocking_times.end(); ++iter) { event_blocking_time_intervals[iter->first]=iter->second.range(); }

    ARIADNE_LOG(2,"event_blocking_time_intervals="<<event_blocking_time_intervals);
    ARIADNE_LOG(2,"non_transverse_events="<<non_transverse_events<<"\n");

    // Compute blocking events
    compute_blockingTime_and_relatedEvents(blocking_events,blocking_time_model,event_blocking_times);
    ARIADNE_LOG(2,"blocking_events="<<blocking_events);
    ARIADNE_LOG(2,"blocking_time="<<blocking_time_model.range()<<"\n");

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
    ARIADNE_LOG(2,"(adjusted)blocking_events="<<blocking_events);
    ARIADNE_LOG(2,"(adjusted)blocking_time="<<blocking_time_model.range()<<"\n");
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
    ARIADNE_LOG(2,"activation_time_intervals="<<activation_time_intervals<<"\n");
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
    ARIADNE_LOG(3,"flow_set_model="<<flow_set_model);
    ARIADNE_LOG(3,"zero_time_model="<<zero_time_model);
    ARIADNE_LOG(3,"blocking_time_model="<<blocking_time_model);
    reachable_set=_toolbox->reachability_step(flow_set_model,zero_time_model,blocking_time_model);
    reach_sets.adjoin(make_pair(location,reachable_set));

	ARIADNE_LOG(2,"reachable_set="<<reachable_set);
    ARIADNE_LOG(2,"reachable_set.argument_size()="<<reachable_set.argument_size());
    ARIADNE_LOG(2,"reachable_set.range()="<<reachable_set.range());
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
    ARIADNE_LOG(2,"previous events = "<<previous_events);
    ARIADNE_LOG(2,"time_range = "<<time_model.range());
    ARIADNE_LOG(2,"time_model generators = "<<time_model.argument_size());
    ARIADNE_LOG(2,"location = "<<location);
    ARIADNE_LOG(2,"box = "<<set_model.range());
    ARIADNE_LOG(2,"generators = "<<set_model.argument_size());
    ARIADNE_LOG(2,"radius = "<<radius(set_model.range())<<"\n");

    ARIADNE_LOG(2,"dynamic = "<<dynamic);
    ARIADNE_LOG(2,"invariants = "<<invariants);
    ARIADNE_LOG(2,"blocking_guards = "<<blocking_guards);
    ARIADNE_LOG(2,"permissive_guards = "<<permissive_guards<<"\n");
}



void
ImageSetHybridEvolver::
_evolution_add_initialSet(
		std::list<EvolutionData>& working_sets,
		const EnclosureType& initial_set,
		Semantics semantics) const
{
    ARIADNE_LOG(3,"initial_set = "<<initial_set);
    DiscreteLocation initial_location;
    ContinuousEnclosureType initial_continuous_set;
    make_lpair(initial_location,initial_continuous_set)=initial_set;
    ARIADNE_LOG(3,"initial_location = "<<initial_location);
    SetModelType initial_set_model=_toolbox->set_model(initial_continuous_set);

	// Check for match between the enclosure cell size and the set size
	ARIADNE_ASSERT_MSG(this->_settings->_reference_enclosure_widths.find(initial_location)->second.size() ==
			initial_set_model.size(), "Error: mismatch between the reference_enclosure_widths size and the set size.");

    ARIADNE_LOG(3,"initial_set_model = "<<initial_set_model);
    TimeModelType initial_time_model=_toolbox->time_model(0.0,Box(initial_set_model.argument_size()));
    ARIADNE_LOG(3,"initial_time_model = "<<initial_time_model);
    TimedSetModelType initial_timed_set_model=join(initial_set_model.models(),initial_time_model);
    ARIADNE_LOG(3,"initial_timed_set_model = "<<initial_timed_set_model);
    working_sets.push_back(EvolutionData(working_sets.size(),initial_location,std::numeric_limits<Float>::infinity(),EventListType(),initial_set_model,initial_time_model));
}

std::map<uint,Vector<Float> >
ImageSetHybridEvolver::
_indexed_set_models_widths(std::list<EvolutionData>& working_sets) const
{

	std::map<uint,Vector<Float> > result;

	for (std::list<EvolutionData>::const_iterator set_it = working_sets.begin(); set_it != working_sets.end(); ++set_it) {
		const Vector<Interval> initial_set_model_range = set_it->set_model().range();

		Vector<Float> widths(initial_set_model_range.size());
		for (uint j=0;j<initial_set_model_range.size();++j)
			widths[j] = initial_set_model_range[j].width();

		result.insert(make_pair(set_it->set_index(),widths));
	}

	return result;
}

bool
ImageSetHybridEvolver::
_is_enclosure_too_large(
		const DiscreteLocation& loc,
		const SetModelType& set_model,
		const Vector<Float>& initial_set_model_widths) const
{
	const Vector<Float>& loc_reference_enclosure_widths =
			this->_settings->_reference_enclosure_widths.find(loc)->second;

	// Identify whether we should use the minimum discretised_enclosure_widths
	bool use_initial_set_as_reference = false;
	for (uint i=0;i<initial_set_model_widths.size();++i)
		if (initial_set_model_widths[i] >= loc_reference_enclosure_widths[i]) {
			use_initial_set_as_reference = true;
			break;
		}

	const Vector<Interval> set_model_range = set_model.range();
	for (uint i=0;i<set_model_range.size();++i) {
		if (use_initial_set_as_reference) {
			if (set_model_range[i].width()/initial_set_model_widths[i] >= this->_settings->maximum_enclosure_widths_ratio())
				return true;
		} else {
			if (set_model_range[i].width() >= loc_reference_enclosure_widths[i]*
					this->_settings->maximum_enclosure_widths_ratio())
				return true;
		}
	}

	return false;
}

void
ImageSetHybridEvolver::
_add_subdivisions(std::list<EvolutionData>& working_sets,
				  const array< TimedSetModelType >& subdivisions,
				  const uint& set_index,
				  const DiscreteLocation& initial_location,
				  const Float& previous_step,
				  const EventListType& previous_events,
				  const uint dimension) const
{
    ARIADNE_LOG(3,"subdivisions.size()="<<subdivisions.size());
    for(uint i=0; i!=subdivisions.size(); ++i) {
        TimedSetModelType const& subdivided_timed_set_model=subdivisions[i];
        ARIADNE_LOG(3,"subdivided_timed_set_model.range()="<<subdivided_timed_set_model.range());
        SetModelType subdivided_set_model=Vector<TaylorModel>(project(subdivided_timed_set_model.models(),range(0,dimension)));
        TimeModelType subdivided_time_model=subdivided_timed_set_model[dimension];
        ARIADNE_LOG(3,"subdivided_set_model.range()="<<subdivided_set_model.range());
        ARIADNE_LOG(3,"subdivided_set_model.radius()*10000="<<radius(subdivided_set_model.range())*10000);
        ARIADNE_LOG(3,"subdivided_time_model.range()="<<subdivided_time_model.range());
        working_sets.push_back(EvolutionData(set_index,initial_location,previous_step,previous_events,subdivided_set_model,subdivided_time_model));
    }
}

void
ImageSetHybridEvolver::
_add_models_subdivisions_autoselect(
		std::list<EvolutionData>& working_sets,
		const uint& set_index,
		const SetModelType& set_model,
		const TimeModelType& time_model,
		const DiscreteLocation& location,
		const Float& previous_step,
		const EventListType& previous_events,
		Semantics semantics) const
{
    uint nd=set_model.dimension();
    SetModelType timed_set_model=join(set_model.models(),time_model);
    array< TimedSetModelType > subdivisions=_toolbox->subdivide(timed_set_model);
    _add_subdivisions(working_sets,subdivisions,set_index,location,previous_step,previous_events,nd);
}


void
ImageSetHybridEvolver::
_add_models_subdivisions_time(
		std::list<EvolutionData>& working_sets,
		const uint& set_index,
		const SetModelType& set_model,
		const TimeModelType& time_model,
		const DiscreteLocation& location,
		const Float& previous_step,
		const EventListType& previous_events,
		Semantics semantics) const
{
    uint nd=set_model.dimension();
    SetModelType timed_set_model=join(set_model.models(),time_model);
    array< TimedSetModelType > subdivisions=_toolbox->subdivide(timed_set_model,nd);
    _add_subdivisions(working_sets,subdivisions,set_index,location,previous_step,previous_events,nd);
}


ImageSetHybridEvolverSettings::ImageSetHybridEvolverSettings(const SystemType& sys)
    : _sys(sys)
{
	set_maximum_step_size(std::numeric_limits<Float>::infinity());
	set_reference_enclosure_widths(getMinimumGridCellWidths(HybridGrid(sys.state_space()),0));
	set_maximum_enclosure_widths_ratio(5.0);
	set_enable_error_rate_enforcement(false);
	set_enable_reconditioning(false);
	set_enable_subdivisions(false);
	set_enable_premature_termination_on_enclosure_size(true);
	set_maximum_number_of_working_sets(0);
}

const std::map<DiscreteLocation,Float>&
ImageSetHybridEvolverSettings::maximum_step_size() const {
	return _maximum_step_size;
}

void
ImageSetHybridEvolverSettings::set_maximum_step_size(const Float& value) {
	ARIADNE_ASSERT_MSG(value > 0, "Error: the maximum step size must be greater than zero.");

    HybridSpace hspace(_sys.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        _maximum_step_size[hs_it->first] = value;
    }
}

void
ImageSetHybridEvolverSettings::set_maximum_step_size(const std::map<DiscreteLocation,Float>& value) {

	for (std::map<DiscreteLocation,Float>::const_iterator it = value.begin(); it != value.end(); ++it) {
	    ARIADNE_ASSERT_MSG(it->second > 0, "Error: the maximum step size for location " << it->first.name() << " is zero.");
	}

	_maximum_step_size = value;
}

const HybridFloatVector&
ImageSetHybridEvolverSettings::reference_enclosure_widths() const {
	return _reference_enclosure_widths;
}

void
ImageSetHybridEvolverSettings::set_reference_enclosure_widths(const Float& value) {
    HybridSpace hspace(_sys.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
    	_reference_enclosure_widths[hs_it->first] = Vector<Float>(hs_it->second,value);
    }
}

void
ImageSetHybridEvolverSettings::set_reference_enclosure_widths(const Vector<Float>& value) {
    HybridSpace hspace(_sys.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
    	_reference_enclosure_widths[hs_it->first] = value;
    }
}

void
ImageSetHybridEvolverSettings::set_reference_enclosure_widths(const HybridFloatVector& value) {
	_reference_enclosure_widths = value;
}

const Float&
ImageSetHybridEvolverSettings::maximum_enclosure_widths_ratio() const {
	return _maximum_enclosure_widths_ratio;
}

void
ImageSetHybridEvolverSettings::set_maximum_enclosure_widths_ratio(const Float& value) {
	_maximum_enclosure_widths_ratio = value;
}

const bool&
ImageSetHybridEvolverSettings::enable_error_rate_enforcement() const {
    return _enable_error_rate_enforcement;
}
void
ImageSetHybridEvolverSettings::set_enable_error_rate_enforcement(const bool& value) {
    _enable_error_rate_enforcement = value;
}

const bool&
ImageSetHybridEvolverSettings::enable_reconditioning() const {
	return _enable_reconditioning;
}
void
ImageSetHybridEvolverSettings::set_enable_reconditioning(const bool& value) {
	_enable_reconditioning = value;
}

const bool&
ImageSetHybridEvolverSettings::enable_subdivisions() const {
	return _enable_subdivisions;
}
void
ImageSetHybridEvolverSettings::set_enable_subdivisions(const bool& value) {
	_enable_subdivisions = value;
}

const bool&
ImageSetHybridEvolverSettings::enable_premature_termination_on_enclosure_size() const {
	return _enable_premature_termination_on_enclosure_size;
}
void
ImageSetHybridEvolverSettings::set_enable_premature_termination_on_enclosure_size(const bool& value) {
	_enable_premature_termination_on_enclosure_size = value;
}

const unsigned int&
ImageSetHybridEvolverSettings::maximum_number_of_working_sets() const {
	return _maximum_number_of_working_sets;
}

void
ImageSetHybridEvolverSettings::set_maximum_number_of_working_sets(const unsigned int& value) {
	_maximum_number_of_working_sets = value;
}

std::ostream&
operator<<(std::ostream& os, const ImageSetHybridEvolverSettings& s)
{
    os << "ImageSetHybridEvolverSettings"
       << ",\n  maximum_step_size=" << s.maximum_step_size()
       << ",\n  reference_enclosure_widths=" << s.reference_enclosure_widths()
       << ",\n  maximum_enclosure_widths_ratio=" << s.maximum_enclosure_widths_ratio()
	   << ",\n  enable_error_rate_enforcement=" << s.enable_error_rate_enforcement()
	   << ",\n  enable_reconditioning=" << s.enable_reconditioning()
       << ",\n  enable_subdivisions=" << s.enable_subdivisions()
       << ",\n  enable_premature_termination_on_enclosure_size=" << s.enable_premature_termination_on_enclosure_size()
	   << ",\n  maximum_number_of_working_sets=" << s.maximum_number_of_working_sets()
       << "\n)\n";
    return os;
}

std::ostream&
operator<<(std::ostream& os, const ContinuousStepResult& r)
{
    os << "("
       << "s=" << r.used_step()
       << ", flow=" << r.flow_set_model()
       << ", finishing=" << r.finishing_set_model()
       << ")";
    return os;
}


HybridFloatVector
getMinimumGridCellWidths(
		const HybridGrid& hgrid,
		int maximum_grid_depth)
{
	const Float divider = (1<<maximum_grid_depth);

	HybridFloatVector result;
	for (HybridGrid::locations_const_iterator hg_it = hgrid.locations_begin(); hg_it != hgrid.locations_end(); hg_it++) {
		result.insert(make_pair(hg_it->first,hg_it->second.lengths()/divider));
	}

	return result;
}


std::map<DiscreteLocation,Float>
getHybridMaximumStepSize(
		const HybridFloatVector& hmad,
		const HybridGrid& hgrid,
		int maximum_grid_depth)
{
    const Float divider = (1<<maximum_grid_depth);

	// We choose a coefficient such that an enclosure at maximum size is able to cross
	// urgent transitions in one step.
	const Float coefficient = 2.0;

	std::map<DiscreteLocation,Float> hmss;

	for (HybridFloatVector::const_iterator hfv_it = hmad.begin(); hfv_it != hmad.end(); hfv_it++)
	{
		const uint dim = hfv_it->second.size();
		// For each dimension of the space, if the derivative is not zero,
		// evaluates the ratio between the minimum cell length and the derivative itself
		Float mss = 0.0;
		for (uint i=0;i<dim;i++)
			if (hfv_it->second[i] > 0)
				mss = max(mss,hgrid[hfv_it->first].lengths()[i]/divider/hfv_it->second[i]);

		// If the step size is still zero, it means that the derivatives are zero on each dimension. Hence we can
		// use any large value we want.
		if (mss == 0)
		   mss = std::numeric_limits<Float>::infinity();

		hmss.insert(std::pair<DiscreteLocation,Float>(hfv_it->first,coefficient*mss));
	}

	return hmss;
}


}  // namespace Ariadne

