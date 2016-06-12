/*****************************************************************************************************
 *            cutting_depth.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of the cutting depth under a laser.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef CUTTING_DEPTH_H
#define CUTTING_DEPTH_H

namespace Ariadne {

HybridIOAutomaton getCuttingDepth()
{
    /// Parameters
	RealParameter k("k",0.0001);
	RealParameter Tcut_burn("Tcut_burn",50.0);
	RealParameter Tcut_stable("Tcut_stable",45.0);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("cutting-depth");

    /// Create the discrete states
    DiscreteLocation burning("burning");
    DiscreteLocation stable("stable");

    RealVariable T("T");
    RealVariable z("z");

    automaton.add_input_var(T);
    automaton.add_output_var(z);

    // Events
    DiscreteEvent start_burning("start_burning");
    DiscreteEvent stop_burning("stop_burning");

    automaton.add_internal_event(start_burning);
    automaton.add_internal_event(stop_burning);

	automaton.new_mode(burning);
	automaton.new_mode(stable);

	RealExpression dyn_burning = - k * (T-Tcut_stable);
	RealExpression dyn_stable = 0.0;

	automaton.set_dynamics(burning, z, dyn_burning);
	automaton.set_dynamics(stable, z, dyn_stable);

	/// Transitions
	// Guards
	RealExpression T_greater_Tcut_burn = T - Tcut_burn;
	RealExpression T_lesser_Tcut_stable = Tcut_stable - T;

	automaton.new_forced_transition(start_burning,stable,burning,T_greater_Tcut_burn);
	automaton.new_forced_transition(stop_burning,burning,stable,T_lesser_Tcut_stable);

	return automaton;
}

}

#endif
