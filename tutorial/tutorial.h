/***************************************************************************
 *            tutorial.h
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

#ifndef TUTORIAL_H_
#define TUTORIAL_H_

#include "ariadne.h"

namespace Ariadne {

HybridIOAutomaton getSystem()
{
    // System variables
    RealVariable x("x");
    RealVariable y("y");

    /// Tank automaton

    // Containing automaton
    HybridIOAutomaton tank("tank");

    // Parameters to be used in the automaton definition
    RealParameter a("a",0.02);
    RealParameter b("b",Interval(0.3,0.32863));

    // Locations for discrete states
    DiscreteLocation flow("flow");

    // Registration of the input/output variables
    tank.add_input_var(y);
    tank.add_output_var(x);

    // Registration of the locations
    tank.new_mode(flow);

    // Input/output variables
    tank.add_input_var(y);
    tank.add_output_var(x);

    // Dynamics
    RealExpression dyn = - a * sqrt(x) + b * y;

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
    valve.add_output_var(y);

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
    valve.set_dynamics(idle, y, dynidle);
    valve.set_dynamics(opening, y, dynopening);
    valve.set_dynamics(closing, y, dynclosing);

    // Resets
    std::map< RealVariable, RealExpression> reset_y_identity;
    reset_y_identity[y] = y;
    std::map< RealVariable, RealExpression> reset_y_one;
    reset_y_one[y] = 1.0;
    std::map< RealVariable, RealExpression> reset_y_zero;
    reset_y_zero[y] = 0.0;

    // Guards
    RealExpression y_geq_one = y - 1.0;
    RealExpression y_leq_zero = - y;

    // Registration of the transitions
    valve.new_unforced_transition(e_open, idle, opening, reset_y_identity);
    valve.new_unforced_transition(e_close, idle, closing, reset_y_identity);
    valve.new_forced_transition(e_idle, opening, idle, reset_y_identity, y_geq_one);
    valve.new_forced_transition(e_idle, closing, idle, reset_y_identity, y_leq_zero);

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

#endif /* TUTORIAL_H_ */
