/***************************************************************************
 *            watertank.h
 *
 *  Copyright  2016  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifndef WATERTANK_H_
#define WATERTANK_H_

#include "ariadne.h"

namespace Ariadne {

HybridIOAutomaton getSystem()
{
    // System variables
    RealVariable a("a"); // Valve aperture
    RealVariable x("x"); // Water level

    /// Tank automaton

    // Containing automaton
    HybridIOAutomaton tank("tank");

    // Parameters to be used in the automaton definition
    RealParameter alpha("alpha",0.02);
    RealParameter bfp("bfp",Interval(0.3,0.32863));

    // Locations for discrete states
    DiscreteLocation flow("flow");

    // Registration of the input/output variables
    tank.add_input_var(a);
    tank.add_output_var(x);

    // Registration of the locations
    tank.new_mode(flow);

    // Input/output variables
    tank.add_input_var(a);
    tank.add_output_var(x);

    // Dynamics
    RealExpression dyn = - alpha * sqrt(x) + bfp * a;

    // Registration of the dynamics
    tank.set_dynamics(flow, x, dyn);

    /// Valve automaton

    // Containing automaton
    HybridIOAutomaton valve("valve");

    // Parameters to be used in the automaton definition
    RealParameter T("T",4.0);

    // Locations for discrete states
    DiscreteLocation idle("idle");
    DiscreteLocation opening("opening");
    DiscreteLocation closing("closing");

    // Registration of the input/output variables
    valve.add_output_var(a);

    // Discrete events for transitions
    DiscreteEvent e_open("open");
    DiscreteEvent e_close("close");
    DiscreteEvent e_idle("idle");

    // Registration of the input/internal events
    valve.add_input_event(e_open);
    valve.add_input_event(e_close);
    valve.add_internal_event(e_idle);

    // Dynamics
    RealExpression dynidle = 0.0;
    RealExpression dynopening = 1.0/T;
    RealExpression dynclosing = -1.0/T;

    // Registration of the locations
    valve.new_mode(idle);
    valve.new_mode(opening);
    valve.new_mode(closing);

    // Registration of the dynamics for each location
    valve.set_dynamics(idle, a, dynidle);
    valve.set_dynamics(opening, a, dynopening);
    valve.set_dynamics(closing, a, dynclosing);

    // Transitions

	valve.new_unforced_transition(e_open, idle, opening);
	valve.new_unforced_transition(e_open, opening, opening);
	valve.new_unforced_transition(e_close, idle, closing);
	valve.new_unforced_transition(e_close, closing, closing);
	// when the valve is fully opened go from opening to idle
	RealExpression a_geq_one = a - 1.0;
	std::map<RealVariable,RealExpression> reset_a_one;
	reset_a_one[a] = 1.0;
	valve.new_forced_transition(e_idle, opening, idle, reset_a_one, a_geq_one);
	// when the valve is fully closed go from closing to idle
	RealExpression a_leq_zero = - a;
	std::map<RealVariable,RealExpression> reset_a_zero;
	reset_a_zero[a] = 0.0;
	valve.new_forced_transition(e_idle, closing, idle, reset_a_zero, a_leq_zero);

    /// Controller automaton

    // Containing automaton
    HybridIOAutomaton controller("controller");

    // Parameters to be used in the automaton definition
    RealParameter hmin("hmin",5.75);
    RealParameter hmax("hmax",7.75);
    RealParameter delta("delta",0.1);

    // Locations for discrete states
    DiscreteLocation rising("rising");
    DiscreteLocation falling("falling");

    // Registration of the input/output variables
    controller.add_input_var(x);

    // Two output events (open and close)
    controller.add_output_event(e_open);
    controller.add_output_event(e_close);

    // Registration of the locations
    controller.new_mode(rising);
    controller.new_mode(falling);

    // Invariants
    RealExpression x_leq_hmax = x - hmax - delta;
    RealExpression x_geq_hmin = hmin - delta - x;

    // Registration of the invariants for each location
    controller.new_invariant(rising, x_leq_hmax);
    controller.new_invariant(falling, x_geq_hmin);

    // Guards
    RealExpression x_geq_hmax = x - hmax + delta;
    RealExpression x_leq_hmin = hmin + delta - x;

    // Transitions
    controller.new_unforced_transition(e_close, rising, falling, x_geq_hmax);
    controller.new_unforced_transition(e_open, falling, rising, x_leq_hmin);
	
    /// Composition

    HybridIOAutomaton tank_valve = compose("tank,valve",tank,valve,flow,idle);
    HybridIOAutomaton system = compose("compositional",tank_valve,controller,DiscreteLocation("flow,idle"),rising);

    return system;
}


}

#endif /* WATERTANK_H_ */
