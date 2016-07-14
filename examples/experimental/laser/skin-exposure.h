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
    DiscreteLocation close_from_right_in("close_from_right_in");
    DiscreteLocation close_from_right_out("close_from_right_out");
    DiscreteLocation close_from_left_in("close_from_left_in");
    DiscreteLocation close_from_left_out("close_from_left_out");
    DiscreteLocation far_from_right_in("far_from_right_in");
    DiscreteLocation far_from_right_out("far_from_right_out");
    DiscreteLocation far_from_left_in("far_from_left_in");
    DiscreteLocation far_from_left_out("far_from_left_out");

    RealVariable x("x");
    RealVariable p("p");

    automaton.add_input_var(x);
    automaton.add_output_var(p);

    // Events
    DiscreteEvent laser_comes_from_right("laser_comes_from_right");
    DiscreteEvent laser_comes_from_left("laser_comes_from_left");
    DiscreteEvent laser_leaves_from_right("laser_leaves_from_right");
    DiscreteEvent laser_leaves_from_left("laser_leaves_from_left");
    DiscreteEvent laser_crosses_from_left("laser_crosses_from_left");
    DiscreteEvent laser_crosses_from_right("laser_crosses_from_right");
    DiscreteEvent switch_left("switch_right");
    DiscreteEvent switch_right("switch_left");

    automaton.add_internal_event(laser_comes_from_left);
    automaton.add_internal_event(laser_crosses_from_left);
    automaton.add_internal_event(laser_leaves_from_left);
    automaton.add_internal_event(laser_comes_from_right);
    automaton.add_internal_event(laser_crosses_from_right);
    automaton.add_internal_event(laser_leaves_from_right);
    automaton.add_input_event(switch_left);
    automaton.add_input_event(switch_right);

	automaton.new_mode(close_from_right_in);
	automaton.new_mode(close_from_right_out);
	automaton.new_mode(far_from_right_in);
	automaton.new_mode(far_from_right_out);
	automaton.new_mode(close_from_left_in);
	automaton.new_mode(close_from_left_out);
	automaton.new_mode(far_from_left_in);
	automaton.new_mode(far_from_left_out);

	RealExpression distance = Ariadne::sqr(x-x0);

	RealExpression dyn_close_from_right = velocity*Ariadne::pi<Real>()/L/L * (x-x0) * Ariadne::sin(Ariadne::pi<Real>()/L/L * distance);
	RealExpression dyn_close_from_left = -dyn_close_from_right;
	RealExpression dyn_far = 0.0;

	automaton.set_dynamics(close_from_right_in, p, dyn_close_from_right);
	automaton.set_dynamics(close_from_right_out, p, dyn_close_from_right);
	automaton.set_dynamics(far_from_right_in, p, dyn_far);
	automaton.set_dynamics(far_from_right_out, p, dyn_far);
	automaton.set_dynamics(close_from_left_in, p, dyn_close_from_left);
	automaton.set_dynamics(close_from_left_out, p, dyn_close_from_left);
	automaton.set_dynamics(far_from_left_in, p, dyn_far);
	automaton.set_dynamics(far_from_left_out, p, dyn_far);

	/// Transitions
	// Guards
	RealExpression distance_greater_L = distance - L*L; // distance >= L
	RealExpression distance_lesser_L = L*L - distance; // distance <= L
	RealExpression x_greater_x0 = x-x0; // x >= x0
	RealExpression x_lesser_x0 = x0-x; // x <= x0

	// Resets
	std::map<RealVariable,RealExpression> reset_zero;
	reset_zero[p] = 0.0;
	std::map<RealVariable,RealExpression> reset_one;
	reset_one[p] = 1.0;

	automaton.new_forced_transition(laser_comes_from_right,far_from_right_in,close_from_right_in,distance_lesser_L);
	automaton.new_forced_transition(laser_comes_from_left,far_from_left_in,close_from_left_in,distance_lesser_L);
	automaton.new_forced_transition(laser_crosses_from_right,close_from_right_in,close_from_right_out,reset_one,x_lesser_x0);
	automaton.new_forced_transition(laser_crosses_from_left,close_from_left_in,close_from_left_out,reset_one,x_greater_x0);
	automaton.new_forced_transition(laser_leaves_from_right,close_from_right_out,far_from_right_out,reset_zero,distance_greater_L);
	automaton.new_forced_transition(laser_leaves_from_left,close_from_left_out,far_from_left_out,reset_zero,distance_greater_L);
	automaton.new_unforced_transition(switch_left,far_from_left_out,far_from_right_in);
	automaton.new_unforced_transition(switch_right,far_from_right_out,far_from_left_in);

	return automaton;
}

}

#endif
