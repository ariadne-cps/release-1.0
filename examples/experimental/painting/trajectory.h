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
    RealParameter angle("angle",Ariadne::pi<Real>()*0.0); // Angle of orientation of the spraying pass
    RealParameter d("d",0.01); // The distance between the lines of each pass

    /// Modes

    DiscreteLocation scanning("scanning");
	automaton.new_mode(scanning);

	DiscreteLocation jumping_line("jumping_line");
	automaton.new_mode(jumping_line);

    // Variables

    RealVariable x("x"); // X position of the laser
    RealVariable vx("vx"); // X velocity of the laser
    RealVariable y("y"); // Y position of the laser
    RealVariable vy("vy"); // Y velocity of the laser

    automaton.add_output_var(x);
    automaton.add_output_var(vx);
    automaton.add_output_var(y);
    automaton.add_output_var(vy);

    // Events

    DiscreteEvent stop("stop");
    DiscreteEvent start("start");

    automaton.add_input_event(stop);
    automaton.add_output_event(start);

	// Dynamics

	RealExpression dyn_x = vx;
	automaton.set_dynamics(jumping_line, x, dyn_x);
	automaton.set_dynamics(scanning, x, dyn_x);
	RealExpression dyn_y = vy;
	automaton.set_dynamics(jumping_line, y, dyn_y);
	automaton.set_dynamics(scanning, y, dyn_y);

	RealExpression dyn_vx = 0.0;
	automaton.set_dynamics(jumping_line, vx, dyn_vx);
	automaton.set_dynamics(scanning, vx, dyn_vx);
	RealExpression dyn_vy = 0.0;
	automaton.set_dynamics(jumping_line, vy, dyn_vy);
	automaton.set_dynamics(scanning, vy, dyn_vy);

	// Transitions

	std::map<RealVariable,RealExpression> reset_jump;
	reset_jump[x] = x - d*Ariadne::sin(angle);
	reset_jump[y] = y + d*Ariadne::cos(angle);
	reset_jump[vx] = -vx;
	reset_jump[vy] = -vy;

	automaton.new_unforced_transition(stop,scanning,jumping_line);
	automaton.new_forced_transition(start,jumping_line,scanning,reset_jump);

	return automaton;
}

}

#endif
