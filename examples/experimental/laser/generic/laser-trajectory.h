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
    RealParameter velocity("velocity",0.46); // Velocity in modulus
    RealParameter width("width",0.0046); // Width of the cut

    /// Modes

    DiscreteLocation scanning("scanning");

	automaton.new_mode(scanning);

    // Variables

    RealVariable x("x"); // Position of the laser
    RealVariable vx("vx"); // Velocity of the laser

    automaton.add_output_var(x);
    automaton.add_output_var(vx);

    // Events

    DiscreteEvent switch_left("switch_right");
    DiscreteEvent switch_right("switch_left");

    automaton.add_output_event(switch_left);
    automaton.add_output_event(switch_right);

	// Dynamics

	RealExpression dyn_x = vx;

	automaton.set_dynamics(scanning, x, dyn_x);

	RealExpression dyn_vx = 0.0;

	automaton.set_dynamics(scanning, vx, dyn_vx);

	// Transitions

	RealExpression x_greater_width = x - width; // x >= width
	RealExpression x_lesser_zero = -x; // x <= 0

	std::map<RealVariable,RealExpression> reset_left;
	reset_left[x] = 0.0;
	reset_left[vx] = velocity;
	std::map<RealVariable,RealExpression> reset_right;
	reset_right[x] = width;
	reset_right[vx] = -velocity;

	automaton.new_forced_transition(switch_left,scanning,scanning,reset_right,x_greater_width);
	automaton.new_forced_transition(switch_right,scanning,scanning,reset_left,x_lesser_zero);

	return automaton;
}

}

#endif
