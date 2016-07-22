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
	RealParameter k("k",4.1475e-8);
	RealParameter lambda("lambda",3.3456e3);
	RealParameter mu("mu",2.9595e5);
	RealParameter T0("T0",37.0);
	RealParameter Tevap("Teval",100.0);
	RealParameter z_thr("z_thr",15e-6);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("depth");

    /// Create the discrete states
    DiscreteLocation ablating("ablating");
    DiscreteLocation idle("idle");
    DiscreteLocation carbonization("carbonization");

    RealVariable p("p");
    RealVariable z("z"); // The cutting depth
    RealVariable zi("zi"); // The depth for each pass of the laser

    automaton.add_input_var(p);
    automaton.add_internal_var(zi);
    automaton.add_output_var(z);

    // Events
    DiscreteEvent start_evaporating("start_evaporating");
    DiscreteEvent stop_evaporating("stop_evaporating");
    DiscreteEvent start_carbonization("start_carbonization");

    automaton.add_output_event(stop_evaporating);
    automaton.add_input_event(start_evaporating);
    automaton.add_output_event(start_carbonization);

	automaton.new_mode(ablating);
	automaton.new_mode(carbonization);
	automaton.new_mode(idle);

	/// Invariants
	RealExpression invalid = 1.0;
	automaton.new_invariant(carbonization,invalid);
	/// Dynamics

	RealExpression z_der = k*(mu*p - lambda*(Tevap-T0));

	RealExpression dyn_ablating = z_der;

	RealExpression dyn_z_idle = 0.0;
	RealExpression dyn_zi_idle = - lambda*zi;

	automaton.set_dynamics(ablating, z, dyn_ablating);
	automaton.set_dynamics(ablating, zi, dyn_ablating);
	automaton.set_dynamics(idle, z, dyn_z_idle);
	automaton.set_dynamics(idle, zi, dyn_zi_idle);
	automaton.set_dynamics(carbonization, z, dyn_ablating);
	automaton.set_dynamics(carbonization, zi, dyn_ablating);

	/// Transitions
	// Guards
	RealExpression zi_greater_zi_thr = zi - z_thr; // zi >= z_thr
	RealExpression zi_der_lesser_zero = -z_der; // z' <= 0

	// Resets
	automaton.new_unforced_transition(start_evaporating,idle,ablating);
	automaton.new_forced_transition(stop_evaporating,ablating,idle,zi_der_lesser_zero);
	automaton.new_forced_transition(start_carbonization,ablating,carbonization,zi_greater_zi_thr);

	return automaton;

}

}

#endif
