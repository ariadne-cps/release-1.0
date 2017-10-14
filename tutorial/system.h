/***************************************************************************
 *            system.h
 *
 *  This file provides the system definition.
 *  Specifically, this is watertank system in which a tank with a hole in the
 *  bottom receives an input water flow. Such input flow can
 *  be modulated between zero and its maximum by controlling a valve. The
 *  described controller aims at keeping the water level between an upper threshold
 *  and a lower threshold.
 *
 *  Copyright  2017  Luca Geretti
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

#ifndef TUTORIAL_SYSTEM_H_
#define TUTORIAL_SYSTEM_H_

#include "ariadne.h"

namespace Ariadne {

/*
 * Construction of an automaton takes major six steps:
 *
 * 1. Creation of an automaton
 * 2. Registration of variable on the automaton
 * 3. Registration of events on the automaton
 * 4. Registration of locations as modes of the automaton
 * 5. Registration of dynamics for each mode
 * 6. Registration of transitions from one mode to another mode
 *
 * As a 0-th step, we also need to create (valued) labels: system variables, parameters and events
 *
 * The creation of labels was numbered 0 since they may be shared between multiple automata, and we clearly
 * need to define them only once in the given context, before any usage.
 *
 * While it may be convenient to define all the components of a system on the same context,
 * for sufficiently complex automata it becomes preferable to separate automata within
 * dedicated (header) files. In that case, shared labels must be redefined within each context.
 *
 * Finally, a system is obtained by automated composition of its components.
 */

HybridIOAutomaton getSystem()
{
    // 0: System variables

		RealVariable a("a"); // Valve aperture
		RealVariable x("x"); // Water level

    /// Tank automaton

		// 0: Parameters

			RealParameter alpha("alpha",0.02); // The coefficient for output flow
			RealParameter beta("bfp",Interval(0.3,0.32863)); // The coefficient for input flow, defined as an interval, meaning that all the values are considered

		// 1. Automaton

			HybridIOAutomaton tank("tank");

		// 2. Registration of the input/output variables

			tank.add_input_var(a);
			tank.add_output_var(x);

		// 4. Registration of the locations

			DiscreteLocation flow("flow");

			tank.new_mode(flow);

		/// 5. Registration of the dynamics

			tank.set_dynamics(flow, x, - alpha * x + beta * a);

    /// Valve automaton

		// 0. Parameters

			RealParameter T("T",4.0); // Time constant for opening/closing the valve

		// 1. Automaton

			HybridIOAutomaton valve("valve");

		// 2. Registration of the input/output variables

			valve.add_output_var(a);

		// 3 Registration of the input/internal events

			DiscreteEvent e_open("open");
			DiscreteEvent e_close("close");
			DiscreteEvent e_idle("idle");

			valve.add_input_event(e_open);
			valve.add_input_event(e_close);
			valve.add_internal_event(e_idle);

		// 4. Registration of the locations

			DiscreteLocation idle("idle");
			DiscreteLocation opening("opening");
			DiscreteLocation closing("closing");

			valve.new_mode(idle);
			valve.new_mode(opening);
			valve.new_mode(closing);

		// 5. Registration of the dynamics for each location

			valve.set_dynamics(idle, a, 0.0);
			valve.set_dynamics(opening, a, 1.0/T);
			valve.set_dynamics(closing, a, -1.0/T);

		/// 6. Transitions

			// Guards
			// The library assumes that given a guard g, the relation g >= 0 must hold in the current mode to have a transition
			RealExpression a_geq_one = a - 1.0; // a >= 1
			RealExpression a_leq_zero = - a; // a >= 0

			// Resets
			// We need to define a reset for each output variable of the automaton
			std::map<RealVariable,RealExpression> reset_a_one;
			reset_a_one[a] = 1.0; // a = 1
			std::map<RealVariable,RealExpression> reset_a_zero;
			reset_a_zero[a] = 0.0; // a = 0

			// Forced transitions: transitions which implicitly have complementary guards and consequently
			// force the transition to be taken immediately

			// When the valve is fully opened, go from opening to idle
			valve.new_forced_transition(e_idle, opening, idle, reset_a_one, a_geq_one);
			// When the valve is fully closed go from closing to idle
			valve.new_forced_transition(e_idle, closing, idle, reset_a_zero, a_leq_zero);

			// Unforced transitions: transitions which do not have complementary guards hence do not
			// necessarily force an immediate transition

			// Transitions that depend on input events must be unforced and with no guard or reset
			valve.new_unforced_transition(e_open, idle, opening);
			valve.new_unforced_transition(e_close, idle, closing);

    /// Controller automaton

		// 0. Parameters

			RealParameter hmin("hmin",5.75); // Lower threshold
			RealParameter hmax("hmax",7.75); // Upper threshold
			RealParameter delta("delta",0.1); // Indetermination constant

		// 1. Automaton

			HybridIOAutomaton controller("controller");

		// 2. Registration of the input/output variables

			controller.add_input_var(x);

		// 3. Registration of the events

			controller.add_output_event(e_open);
			controller.add_output_event(e_close);

		// 4. Registration of the locations

			DiscreteLocation rising("rising");
			DiscreteLocation falling("falling");

			controller.new_mode(rising);
			controller.new_mode(falling);

		// 5. Transitions

			// Invariants
			// The library assumes that given an invariant i, the relation i <= 0 must hold in the current mode to allow evolution
			RealExpression x_leq_hmax = x - hmax - delta; // x <= hmax + delta
			RealExpression x_geq_hmin = hmin - delta - x; // x >= hmin - delta

			// Registration of the invariants for each location
			controller.new_invariant(rising, x_leq_hmax);
			controller.new_invariant(falling, x_geq_hmin);

			// Guards
			RealExpression x_geq_hmax = x - hmax + delta; // x >= hmax - delta
			RealExpression x_leq_hmin = hmin + delta - x; // x <= hmin + delta

			// Unforced transitions: here they are used paired with invariants to limit the region of activation of the guard
			controller.new_unforced_transition(e_close, rising, falling, x_geq_hmax);
			controller.new_unforced_transition(e_open, falling, rising, x_leq_hmin);
	
    /// Composition

	/* Composition is obtained by progressively composing two automata. The first argument is the name of
	 * the resulting composition. Please note that the name is actually relevant only for the complete system
	 *
	 * The second and this arguments are the components.
	 *
	 * Then in the fourth and fifth argument we must define an initial location for each component,
	 * in order to have a compact composition which excludes states that would not be reachable from such initial location.
     * Pay attention, during iterative composition, on ordering the initial state correctly in respect
     * to the provided fourth and fifth arguments; discrete location names for composite locations are
	 * created by simply using the comma character between the location names of their components.
     */

    HybridIOAutomaton tank_valve = compose("tank,valve",tank,valve,flow,idle);
    HybridIOAutomaton system = compose("tutorial",tank_valve,controller,DiscreteLocation("flow,idle"),rising);

    return system;
}


}

#endif /* TUTORIAL_SYSTEM_H_ */
