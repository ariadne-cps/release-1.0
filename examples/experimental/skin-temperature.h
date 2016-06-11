/*****************************************************************************************************
 *            skin_temperature.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of the skin temperature under laser cutting.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef SKIN_TEMPERATURE_H_
#define SKIN_TEMPERATURE_H_

namespace Ariadne {

HybridIOAutomaton getSkinTemperature()
{
    /// Parameters
	RealParameter lambda("lambda",16.0);
	RealParameter mu("mu",1000.0);
	RealParameter T0("T0",35.0);
	RealParameter L("L",0.02);
	RealParameter x0("x0",0.0);
	RealParameter y0("y0",0.05);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("skin-T");

    /// Create the discrete states
    DiscreteLocation heated("heated");
    DiscreteLocation resting("resting");

    RealVariable x("x");
    RealVariable y("y");
    RealVariable T("T");

    automaton.add_input_var(x);
    automaton.add_input_var(y);
    automaton.add_output_var(T);

    // Events
    DiscreteEvent laser_comes("laser_comes");
    DiscreteEvent laser_leaves("laser_leaves");

    automaton.add_internal_event(laser_comes);
    automaton.add_internal_event(laser_leaves);

	automaton.new_mode(heated);
	automaton.new_mode(resting);

	RealExpression distance = Ariadne::sqr(x-x0) + Ariadne::sqr(y-y0);

	RealExpression dyn_heated = - lambda * (T-T0) + 0.5 * mu * (Ariadne::cos(Ariadne::pi<Real>()/L/L * distance) + 1.0);
	RealExpression dyn_resting = - lambda * (T-T0);

	automaton.set_dynamics(heated, T, dyn_heated);
	automaton.set_dynamics(resting, T, dyn_resting);

	/// Transitions
	// Guards
	RealExpression distance_greater_L = distance - L*L;
	RealExpression distance_lesser_L = L*L - distance;

	automaton.new_forced_transition(laser_comes,resting,heated,distance_lesser_L);
	automaton.new_forced_transition(laser_leaves,heated,resting,distance_greater_L);

	return automaton;
}

}

#endif
