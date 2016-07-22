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
    RealParameter half_width("half_width",0.0023); // Half width of the cut

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

	RealExpression x_greater_half_width = x - half_width; // x >= half_width
	RealExpression x_lesser_minus_half_width = -x - half_width; // x <= -half_width

	std::map<RealVariable,RealExpression> reset_minus_half;
	reset_minus_half[x] = -half_width;
	reset_minus_half[vx] = velocity;
	std::map<RealVariable,RealExpression> reset_half;
	reset_half[x] = half_width;
	reset_half[vx] = -velocity;

	automaton.new_forced_transition(switch_left,scanning,scanning,reset_half,x_greater_half_width);
	automaton.new_forced_transition(switch_right,scanning,scanning,reset_minus_half,x_lesser_minus_half_width);

	return automaton;
}

}

#endif
