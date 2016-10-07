/*****************************************************************************************************
 *            spray.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the diffusion of the spray.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef SPRAY_H_
#define SPRAY_H_

namespace Ariadne {

HybridIOAutomaton getSpray()
{
    /// Parameters
	RealParameter L("L",0.005);
	RealParameter x0("x0",0.00);
	RealParameter y0("y0",0.00);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("spray");

    /// Create the discrete states
    DiscreteLocation close("close");
    DiscreteLocation far("far");

    RealVariable x("x");
    RealVariable y("y");
    RealVariable vx("vx");
    RealVariable vy("vy");
    RealVariable s("s");

    automaton.add_input_var(x);
    automaton.add_input_var(y);
    automaton.add_input_var(vx);
    automaton.add_input_var(vy);
    automaton.add_output_var(s);

    // Events
    DiscreteEvent comes("comes");
    DiscreteEvent leaves("leaves");

    automaton.add_internal_event(comes);
    automaton.add_internal_event(leaves);

	automaton.new_mode(close);
	automaton.new_mode(far);

	RealExpression distance = Ariadne::sqr(x-x0) + Ariadne::sqrt(y-y0);

	RealExpression dyn_close = -1.0*Ariadne::pi<Real>()/L/L * (vx*(x-x0) + vy*(y-y0)) * Ariadne::sin(Ariadne::pi<Real>()/L/L * distance);
	RealExpression dyn_far = 0.0;

	automaton.set_dynamics(far, s, dyn_far);
	automaton.set_dynamics(close, s, dyn_close);

	/// Transitions
	// Guards
	RealExpression distance_greater_L = distance - L*L; // distance >= L
	RealExpression distance_lesser_L = L*L - distance; // distance <= L

	// Resets
	std::map<RealVariable,RealExpression> reset_zero;
	reset_zero[s] = 0.0;

	automaton.new_forced_transition(comes,far,close,reset_zero,distance_lesser_L);
	automaton.new_forced_transition(leaves,close,far,reset_zero,distance_greater_L);

	return automaton;
}

}

#endif
