/*****************************************************************************************************
 *            trajectory.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the trajectory of the sprayer.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef TRAJECTORY_H_
#define TRAJECTORY_H_

namespace Ariadne {

HybridIOAutomaton getTrajectory()
{
    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("trajectory");

    // Parameters
    RealParameter d("d",0.005); // The distance between the lines of each pass
    RealParameter width("width",0.05); // The width of the scan

    /// Modes

    DiscreteLocation scanning("scanning");
	automaton.new_mode(scanning);

    // Variables

    RealVariable x("x"); // X position of the laser
    RealVariable vx("vx"); // X velocity of the laser
    RealVariable y("y"); // Y position of the laser

    automaton.add_output_var(x);
    automaton.add_output_var(vx);
    automaton.add_output_var(y);

    // Events

    DiscreteEvent switch_left("switch_left");
    DiscreteEvent switch_right("switch_right");

    automaton.add_internal_event(switch_left);
    automaton.add_internal_event(switch_right);

	// Dynamics

	RealExpression dyn_x = vx;
	automaton.set_dynamics(scanning, x, dyn_x);

	RealExpression dyn_vx = 0.0;
	automaton.set_dynamics(scanning, vx, dyn_vx);

	RealExpression dyn_y = 0.0;
	automaton.set_dynamics(scanning, y, dyn_y);

	// Transitions

	RealExpression x_greater_width = x - width; // x >= width
	RealExpression x_lesser_zero = -x; // x <= 0

	std::map<RealVariable,RealExpression> reset_switch_left;
	reset_switch_left[x] = width;
	reset_switch_left[y] = y + d;
	reset_switch_left[vx] = -vx;

	std::map<RealVariable,RealExpression> reset_switch_right;
	reset_switch_right[x] = 0.0;
	reset_switch_right[y] = y + d;
	reset_switch_right[vx] = -vx;

	automaton.new_forced_transition(switch_left,scanning,scanning,reset_switch_left,x_greater_width);
	automaton.new_forced_transition(switch_right,scanning,scanning,reset_switch_right,x_lesser_zero);

	return automaton;
}

}

#endif
