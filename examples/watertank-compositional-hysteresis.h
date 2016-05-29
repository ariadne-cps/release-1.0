/***************************************************************************
 *            watertank-compositional-hysteresis.h
 *
 *  Copyright  2011  Luca Geretti
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

#ifndef WATERTANK_COMPOSITIONAL_HYSTERESIS_H_
#define WATERTANK_COMPOSITIONAL_HYSTERESIS_H_

#include "ariadne.h"

namespace Ariadne {

HybridIOAutomaton getWatertankCompositionalHysteresis()
{
    /// Set the system parameters
	RealParameter alpha("alpha",0.02);
	RealParameter bfp("bfp",Interval(0.3,0.32863));
	RealParameter T("T",4.0);
	RealParameter hmin("hmin",5.75);
	RealParameter hmax("hmax",7.75);
	RealParameter delta("delta",0.1);

    // System variables
    RealVariable x("x");    // water level
    RealVariable a("a");    // valve aperture

    // Create the tank automaton

	    HybridIOAutomaton tank("tank");

		// States
		DiscreteLocation flow("flow");

		// Add the input/output variables
    	tank.add_input_var(a);
    	tank.add_output_var(x);

		// Only one state with no transitions and no invariants
		RealExpression dyn = - alpha * x + bfp * a;
		tank.new_mode(flow);
		tank.set_dynamics(flow, x, dyn);

	// Create the valve automaton

		HybridIOAutomaton valve("valve");

		// States
		DiscreteLocation idle("idle");
		DiscreteLocation opening("opening");
		DiscreteLocation closing("closing");

		// The valve has one output var (the valve aperture)
		valve.add_output_var(a);

		// Two input events (open and close) and one internal event
		DiscreteEvent e_open("open");
		valve.add_input_event(e_open);
		DiscreteEvent e_close("close");
		valve.add_input_event(e_close);
		DiscreteEvent e_idle("idle");
		valve.add_internal_event(e_idle);

		// Three states:
		// Idle (valve either fully closed or fully opened)
		RealExpression dynidle = 0.0;
		valve.new_mode(idle);
		valve.set_dynamics(idle, a, dynidle);
		// Opening (valve is opening)
		valve.new_mode(opening);
		RealExpression dynopening = 1.0/T;
		valve.set_dynamics(opening, a, dynopening);
		// Closing (valve is closing)
		valve.new_mode(closing);
		RealExpression dynclosing = -1.0/T;
		valve.set_dynamics(closing, a, dynclosing);

		// Transitions

		// when open is received, go to opening
		valve.new_unforced_transition(e_open, idle, opening);
		valve.new_unforced_transition(e_open, opening, opening);
		 // when closed is received, go to closing
		valve.new_unforced_transition(e_close, idle, closing);
		valve.new_unforced_transition(e_close, closing, closing);
		// when the valve is fully opened go from opening to idle
		RealExpression y_geq_one = a - 1.0;
		valve.new_forced_transition(e_idle, opening, idle, y_geq_one);
		// when the valve is fully closed go from closing to idle
		RealExpression y_leq_zero = - a;
		valve.new_forced_transition(e_idle, closing, idle, y_leq_zero);

	// Create the controller automaton

	    HybridIOAutomaton controller("controller");

		// States
		DiscreteLocation rising("rising");
		DiscreteLocation falling("falling");

		// The valve has one input var (the water level)
		controller.add_input_var(x);
		// Two output events (open and close)
		controller.add_output_event(e_open);
		controller.add_output_event(e_close);

		// Two states:
		// Rising (water level is increasing)
		controller.new_mode(rising);
		 // Falling (water level is decreasing)
		controller.new_mode(falling);

		// Transitions
		// when the water is greater than hmax, send a close command
		RealExpression x_geq_hmax = x - hmax + delta;
		controller.new_unforced_transition(e_close, rising, falling, x_geq_hmax);
		// Add the invariant x < hmax + delta to rising
		RealExpression x_leq_hmax = x - hmax - delta;
		controller.new_invariant(rising, x_leq_hmax);

		// when the water is lower than hmin, send a open command
		RealExpression x_leq_hmin = hmin + delta - x;
		controller.new_unforced_transition(e_open, falling, rising, x_leq_hmin);
		// Add the invariant x > hmin - delta to falling
		RealExpression x_geq_hmin = hmin - delta - x;
		controller.new_invariant(falling, x_geq_hmin);

	/// Compose the automata
	HybridIOAutomaton tank_valve = compose("tank,valve",tank,valve,flow,idle);
	HybridIOAutomaton system = compose("watertank-hy",tank_valve,controller,DiscreteLocation("flow,idle"),rising);

	return system;
}


}

#endif /* WATERTANK_COMPOSITIONAL_HYSTERESIS_H_ */
