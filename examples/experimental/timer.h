/*****************************************************************************************************
 *            timer.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides an automaton with a single time variable increasing indefinitely, useful to keep track of time.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef TIMER_H_
#define TIMER_H_

namespace Ariadne {

HybridIOAutomaton getTimer()
{
    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("timer");

    /// Create the discrete states
    DiscreteLocation work("work");

    RealVariable t("t");

    automaton.add_output_var(t);

	automaton.new_mode(work);

	RealExpression dyn = 1.0;
	automaton.set_dynamics(work, t, dyn);

	return automaton;
}

}

#endif
