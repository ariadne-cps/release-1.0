/*****************************************************************************************************
 *            skin-exposure.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of a measure of power over the skin point when the laser is close.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef SKIN_EXPOSURE_H_
#define SKIN_EXPOSURE_H_

namespace Ariadne {

HybridIOAutomaton getSkinExposure()
{
    /// Parameters
	RealParameter L("L",0.002);
	RealParameter x0("x0",0.0);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("skin-exposure");

    /// Create the discrete states
    DiscreteLocation close("close");
    DiscreteLocation far("far");

    RealVariable x("x");
    RealVariable p("p");

    automaton.add_input_var(x);
    automaton.add_output_var(p);

    // Events
    DiscreteEvent laser_comes("laser_comes");
    DiscreteEvent laser_leaves("laser_leaves");

    automaton.add_internal_event(laser_comes);
    automaton.add_internal_event(laser_leaves);

	automaton.new_mode(close);
	automaton.new_mode(far);

	RealExpression distance = Ariadne::sqr(x-x0);

	RealExpression dyn_close = -Ariadne::pi<Real>()/L/L * (x-x0) * Ariadne::sin(Ariadne::pi<Real>()/L/L * distance);
	RealExpression dyn_far = 0.0;

	automaton.set_dynamics(close, p, dyn_close);
	automaton.set_dynamics(far, p, dyn_far);

	/// Transitions
	// Guards
	RealExpression distance_greater_L = distance - L*L; // distance >= L
	RealExpression distance_lesser_L = L*L - distance; // distance <= L

	automaton.new_forced_transition(laser_comes,far,close,distance_lesser_L);
	automaton.new_forced_transition(laser_leaves,close,far,distance_greater_L);

	return automaton;
}

}

#endif
