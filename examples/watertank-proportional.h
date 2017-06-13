/***************************************************************************
 *            watertank-proportional.h
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

#ifndef WATERTANK_PROPORTIONAL_H_
#define WATERTANK_PROPORTIONAL_H_

#include "ariadne.h"

namespace Ariadne {

HybridIOAutomaton getWatertankProportional()
{
    /// Set the system parameters
	RealParameter alpha("alpha",0.065); // The constant defining the decrease rate of the tank level
	RealParameter tau("tau",1.25); // The characteristic time for the opening/closing of the valve
	RealParameter ref("ref",6.75); // A reference tank level
	RealParameter bfp("bfp",0.3); // The product beta*f(p)
	RealParameter Kp("Kp",10.0); // The gain of the proportional controller
	RealParameter delta("delta",Interval(-0.00,0.00)); // An indeterminacy in guards evaluation

	// System variables
	RealVariable x("x"); // water level
	RealVariable a("a"); // valve aperture

    // Create the tank automaton

	    HybridIOAutomaton tank("tank");

		// States
		DiscreteLocation flow("flow");

		// Add the input/output variables
    	tank.add_input_var(a);
    	tank.add_output_var(x);

		RealExpression flow_dyn = - alpha * sqrt(x) + bfp * a;
		tank.new_mode(flow);
		tank.set_dynamics(flow, x, flow_dyn);

		// Invariants
	    RealExpression x_geq_zero = -x;   // x >= 0
	    tank.new_invariant(flow,x_geq_zero);

	/// Create the controller automaton

		HybridIOAutomaton controller("controller");

		// Add the input/output variables
		controller.add_input_var(x);
		controller.add_output_var(a);

		// States
		DiscreteLocation opening("opening");
		DiscreteLocation stabilising("stabilising");
		DiscreteLocation closing("closing");

		// Create the discrete events
		DiscreteEvent stabilise_after_closing("stabilise_after_closing");
		DiscreteEvent stabilise_after_opening("stabilise_after_opening");
		DiscreteEvent close("close");
		DiscreteEvent open("open");

	    // Create the dynamics
	    RealExpression opening_d = (1-a)/tau;
	    RealExpression closing_d = -a/tau;
	    RealExpression stabilising_d = (Kp*(ref-x) - a)/tau;

	    // Dynamics at the different modes
	    controller.new_mode(opening);
	    controller.set_dynamics(opening, a, opening_d);
	    controller.new_mode(closing);
	    controller.set_dynamics(closing, a, closing_d);
	    controller.new_mode(stabilising);
	    controller.set_dynamics(stabilising, a, stabilising_d);

		// Invariants
	    RealExpression a_geq_zero = -a;   // a >= 0
	    RealExpression a_leq_one = a-1.0;   // a <= 1

	    // Create the guards
	    RealExpression x_lesser_ref_minus_delta = -x-delta+ref; // x <= ref - Delta
	    RealExpression x_greater_ref_minus_delta = x-delta-ref; // x >= ref + Delta
	    RealExpression x_lesser_ref_kp_minus_delta = -x+ref-1.0/Kp-delta; // x <= ref - 1/Kp - Delta
	    RealExpression x_greater_ref_kp_minus_delta = x-ref+1.0/Kp-delta; // x >= ref - 1/Kp + Delta

	    // Transitions
	    controller.new_forced_transition(stabilise_after_closing,closing,stabilising,x_lesser_ref_minus_delta);
	    controller.new_forced_transition(stabilise_after_opening,opening,stabilising,x_greater_ref_kp_minus_delta);
	    controller.new_forced_transition(close,stabilising,closing,x_greater_ref_minus_delta);
	    controller.new_forced_transition(open,stabilising,opening,x_lesser_ref_kp_minus_delta);

	/// Composition
	HybridIOAutomaton system = compose("watertank-pr",tank,controller,flow,stabilising);

	return system;
}


}

#endif /* WATERTANK_PROPORTIONAL_H_ */
