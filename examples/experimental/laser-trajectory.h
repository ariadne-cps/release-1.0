/*****************************************************************************************************
 *            laser-trajectory.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the trajectory of the laser.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef LASER_TRAJECTORY_H_
#define LASER_TRAJECTORY_H_

namespace Ariadne {

HybridIOAutomaton getLaserTrajectory()
{
    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("trajectory");

    // Parameters
    RealParameter accel("accel",300.0); // Acceleration (modulus) of the laser
    RealParameter half_width("half_width",0.023); // Half width of the cut (excluding the accel/decel tails)

    /// Modes

    DiscreteLocation accelerating_right("accelerating_right");
    DiscreteLocation decelerating_right("decelerating_right");
    DiscreteLocation accelerating_left("accelerating_left");
    DiscreteLocation decelerating_left("decelerating_left");
    DiscreteLocation passing_right("passing_right");
    DiscreteLocation passing_left("passing_left");

	automaton.new_mode(accelerating_right);
	automaton.new_mode(decelerating_right);
	automaton.new_mode(accelerating_left);
	automaton.new_mode(decelerating_left);
	automaton.new_mode(passing_right);
	automaton.new_mode(passing_left);

    // Variables

    RealVariable x("x"); // Position of the laser
    RealVariable vx("vx"); // Speed of the laser

    automaton.add_internal_var(vx);
    automaton.add_output_var(x);

    // Events

    DiscreteEvent accelerate_right("accelerate_right");
    DiscreteEvent stop_accelerating_right("stop_accelerating_right");
    DiscreteEvent decelerate_right("decelerate_right");
    DiscreteEvent accelerate_left("accelerate_left");
    DiscreteEvent stop_accelerating_left("stop_accelerating_left");
    DiscreteEvent decelerate_left("decelerate_left");

	// Dynamics

	RealExpression dyn_vx_accel = accel;
	RealExpression dyn_vx_decel = -accel;
	RealExpression dyn_vx_passing = 0.0;

	automaton.set_dynamics(accelerating_right, vx, dyn_vx_accel);
	automaton.set_dynamics(decelerating_right, vx, dyn_vx_decel);
	automaton.set_dynamics(accelerating_left, vx, dyn_vx_decel);
	automaton.set_dynamics(decelerating_left, vx, dyn_vx_accel);
	automaton.set_dynamics(passing_right, vx, dyn_vx_passing);
	automaton.set_dynamics(passing_left, vx, dyn_vx_passing);

	RealExpression dyn_x = vx;
	automaton.set_dynamics(accelerating_right, x, dyn_x);
	automaton.set_dynamics(decelerating_right, x, dyn_x);
	automaton.set_dynamics(accelerating_left, x, dyn_x);
	automaton.set_dynamics(decelerating_left, x, dyn_x);
	automaton.set_dynamics(passing_right, x, dyn_x);
	automaton.set_dynamics(passing_left, x, dyn_x);

	// Transitions

	RealExpression x_greater_half_width = x - half_width; // x >= half_width
	RealExpression x_lesser_half_width = -x + half_width; // x <= half_width
	RealExpression x_lesser_minus_half_width = -x - half_width; // x <= -half_width
	RealExpression x_greater_minus_half_width = x + half_width; // x >= -half_width
	RealExpression vx_lesser_zero = -vx; // vx <= 0;
	RealExpression vx_greater_zero = vx; // vx >= 0;

	std::map<RealVariable,RealExpression> reset_left;
	reset_left[x] = x;
	reset_left[vx] = 0.0;
	std::map<RealVariable,RealExpression> reset_right;
	reset_right[x] = x;
	reset_right[vx] = 0.0;

	automaton.new_forced_transition(accelerate_right,decelerating_left,accelerating_right,reset_left,vx_greater_zero);
	automaton.new_forced_transition(accelerate_left,decelerating_right,accelerating_left,reset_right,vx_lesser_zero);
	automaton.new_forced_transition(stop_accelerating_right,accelerating_right,passing_right,x_greater_minus_half_width);
	automaton.new_forced_transition(stop_accelerating_left,accelerating_left,passing_left,x_lesser_half_width);
	automaton.new_forced_transition(decelerate_right,passing_right,decelerating_right,x_greater_half_width);
	automaton.new_forced_transition(decelerate_left,passing_left,decelerating_left,x_lesser_minus_half_width);

	return automaton;
}

}

#endif
