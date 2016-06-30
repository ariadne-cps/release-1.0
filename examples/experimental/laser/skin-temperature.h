/*****************************************************************************************************
 *            skin_temperature.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of the skin temperature under laser cutting.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef SKIN_TEMPERATURE_H_
#define SKIN_TEMPERATURE_H_

namespace Ariadne {

HybridIOAutomaton getSkinTemperature()
{
    /// Parameters
	RealParameter lambda("lambda",5000.0);
	RealParameter mu("mu",1000000.0);
	RealParameter T0("T0",35.0);
	RealParameter Tevap("Tevap",100.0);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("skin-T");

    /// Create the discrete states
    DiscreteLocation varying("varying");
    DiscreteLocation pre_evaporating("pre_evaporating");
    DiscreteLocation evaporating("evaporating");

    RealVariable p("p");
    RealVariable T("T");

    automaton.add_input_var(p);
    automaton.add_output_var(T);

    // Events
    DiscreteEvent check_pre_evaporating("check_pre_evaporating");
    DiscreteEvent start_evaporating("start_evaporating");
    DiscreteEvent stop_evaporating("stop_evaporating");

    automaton.add_output_event(start_evaporating);
    automaton.add_input_event(stop_evaporating);
    automaton.add_internal_event(check_pre_evaporating);

    automaton.new_mode(varying);
    automaton.new_mode(pre_evaporating);
	automaton.new_mode(evaporating);

	RealExpression dyn_varying = mu*p - lambda*(T-T0);
	RealExpression dyn_evaporating = 0.0;

	automaton.set_dynamics(varying, T, dyn_varying);
	automaton.set_dynamics(evaporating, T, dyn_evaporating);
	automaton.set_dynamics(pre_evaporating, T, dyn_varying);

	/// Transitions
	// Guards
	RealExpression T_greater_Tevap = T - Tevap; // T >= Tevap
	RealExpression Tder_greater_zero = mu*p - lambda*(T+1.0-T0); // T' >= 0

	// Resets
	std::map<RealVariable,RealExpression> reset_evap;
	reset_evap[T] = Tevap;

	automaton.new_forced_transition(check_pre_evaporating,varying,pre_evaporating,Tder_greater_zero);
	automaton.new_forced_transition(start_evaporating,pre_evaporating,evaporating,reset_evap,T_greater_Tevap);
	automaton.new_unforced_transition(stop_evaporating,evaporating,varying);

	return automaton;
}

}

#endif
