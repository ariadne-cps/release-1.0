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
	RealParameter velocity("velocity",0.46);
	RealParameter L("L",0.00025);
	RealParameter x0("x0",0.00);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("skin-exposure");

    /// Create the discrete states
    DiscreteLocation close("close");
    DiscreteLocation far("far");

    RealVariable x("x");
    RealVariable vx("vx");
    RealVariable p("p");

    automaton.add_input_var(x);
    automaton.add_input_var(vx);
    automaton.add_output_var(p);

    // Events
    DiscreteEvent comes("comes");
    DiscreteEvent leaves("leaves");

    automaton.add_internal_event(comes);
    automaton.add_internal_event(leaves);

	automaton.new_mode(close);
	automaton.new_mode(far);

	RealExpression distance = Ariadne::sqr(x-x0);

	RealExpression dyn_close = -vx*Ariadne::pi<Real>()/L/L * (x-x0) * Ariadne::sin(Ariadne::pi<Real>()/L/L * distance);
	RealExpression dyn_far = 0.0;

	automaton.set_dynamics(far, p, dyn_far);
	automaton.set_dynamics(close, p, dyn_close);

	/// Transitions
	// Guards
	RealExpression distance_greater_L = distance - L*L; // distance >= L
	RealExpression distance_lesser_L = L*L - distance; // distance <= L

	// Resets
	std::map<RealVariable,RealExpression> reset_zero;
	reset_zero[p] = 0.0;

	automaton.new_forced_transition(comes,far,close,reset_zero,distance_lesser_L);
	automaton.new_forced_transition(leaves,close,far,reset_zero,distance_greater_L);

	return automaton;
}

}

#endif
