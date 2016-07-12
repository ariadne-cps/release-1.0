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
	RealParameter x0("x0",0.0);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("skin-exposure");

    /// Create the discrete states
    DiscreteLocation close_from_right("close_from_right");
    DiscreteLocation close_from_left("close_from_left");
    DiscreteLocation far_from_right("far_from_right");
    DiscreteLocation far_from_left("far_from_left");

    RealVariable x("x");
    RealVariable p("p");

    automaton.add_input_var(x);
    automaton.add_output_var(p);

    // Events
    DiscreteEvent laser_comes_from_right("laser_comes_from_right");
    DiscreteEvent laser_comes_from_left("laser_comes_from_left");
    DiscreteEvent laser_leaves_from_right("laser_leaves_from_right");
    DiscreteEvent laser_leaves_from_left("laser_leaves_from_left");

    automaton.add_internal_event(laser_comes_from_left);
    automaton.add_internal_event(laser_leaves_from_left);
    automaton.add_internal_event(laser_comes_from_right);
    automaton.add_internal_event(laser_leaves_from_right);

	automaton.new_mode(close_from_right);
	automaton.new_mode(far_from_right);
	automaton.new_mode(close_from_left);
	automaton.new_mode(far_from_left);

	RealExpression distance = Ariadne::sqr(x-x0);

	RealExpression dyn_close_from_right = velocity*Ariadne::pi<Real>()/L/L * (x-x0) * Ariadne::sin(Ariadne::pi<Real>()/L/L * distance);
	RealExpression dyn_close_from_left = -dyn_close_from_right;
	RealExpression dyn_far = 0.0;

	automaton.set_dynamics(close_from_right, p, dyn_close_from_right);
	automaton.set_dynamics(far_from_right, p, dyn_far);
	automaton.set_dynamics(close_from_left, p, dyn_close_from_left);
	automaton.set_dynamics(far_from_left, p, dyn_far);

	/// Transitions
	// Guards
	RealExpression distance_greater_L = distance - L*L; // distance >= L
	RealExpression distance_lesser_L = L*L - distance; // distance <= L

	// Resets
	std::map<RealVariable,RealExpression> reset_zero;
	reset_zero[p] = 0.0;

	automaton.new_forced_transition(laser_comes_from_right,far_from_right,close_from_right,distance_lesser_L);
	automaton.new_forced_transition(laser_comes_from_left,far_from_left,close_from_left,distance_lesser_L);
	automaton.new_forced_transition(laser_leaves_from_right,close_from_right,far_from_left,reset_zero,distance_greater_L);
	automaton.new_forced_transition(laser_leaves_from_left,close_from_left,far_from_right,reset_zero,distance_greater_L);

	return automaton;
}

}

#endif
