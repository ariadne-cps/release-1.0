/*****************************************************************************************************
 *            timeout.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides an automaton with a clock for handling a timeout.
 *
 *****************************************************************************************************/

#include "timeout.h"

#ifndef TIMER_H_
#define TIMER_H_

namespace Ariadne {

HybridIOAutomaton getTimeout()
{
    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("timeout");

    RealParameter stop_time("stop_time",1.0);

    /// Create the discrete states
    DiscreteLocation running("running");
    DiscreteLocation stopped("stopped");

    /// Create the discrete events
    DiscreteEvent start("start");
    DiscreteEvent stop("stop");

    automaton.add_input_event(start);
    automaton.add_output_event(stop);

    RealVariable clk("clk");
    automaton.add_internal_var(clk);

	automaton.new_mode(running);
	automaton.new_mode(stopped);

	RealExpression dyn_running = 1.0;
	automaton.set_dynamics(running, clk, dyn_running);

	RealExpression dyn_stopped = 0.0;
	automaton.set_dynamics(stopped, clk, dyn_stopped);

	RealExpression stop_time_hit = clk - stop_time; // clk >= stop_time

	// Resets
	std::map<RealVariable,RealExpression> reset_zero;
	reset_zero[clk] = 0.0;

	automaton.new_forced_transition(stop,running,stopped,reset_zero,stop_time_hit);
	automaton.new_unforced_transition(start,stopped,running);

	return automaton;
}

}

#endif
