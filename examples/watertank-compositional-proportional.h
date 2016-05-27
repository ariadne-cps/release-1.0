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
	RealParameter alpha("a",0.065); // The constant defining the decrease rate of the tank level
	RealParameter tau("tau",1.25); // The characteristic time for the opening/closing of the valve
	RealParameter ref("ref",6.75); // A reference tank level
	RealParameter bfp("bfp",0.3); // The product beta*f(p)
	RealParameter Kp("Kp",0.6); // The gain of the proportional controller
	RealParameter delta("delta",Interval(-0.02,0.02)); // An indeterminacy in guards evaluation

	// System variables
	RealVariable x("x"); // water level
	RealVariable a("a"); // valve aperture
	RealVariable w("w"); // control signal

    // Create the tank automaton

	    HybridIOAutomaton tank("tank");

		// States
		DiscreteLocation flow("flow");

		// Add the input/output variables
    	tank.add_input_var(a);
    	tank.add_output_var(x);

		RealExpression flow_dyn = - alpha * x + bfp * a;
		tank.new_mode(flow);
		tank.set_dynamics(flow, x, flow_dyn);

		// Invariants
	    RealExpression x_geq_zero = -x;   // x >= 0
	    //tank.new_invariant(flow,x_geq_zero);

    /// Create the valve automaton

		HybridIOAutomaton valve("valve");

		// States
		DiscreteLocation modulate("modulate");

		// Add the input/output variables
		valve.add_input_var(w);
		valve.add_output_var(a);

		RealExpression modulate_dyn = (w-a)/tau;
		valve.new_mode(modulate);
		valve.set_dynamics(modulate, a, modulate_dyn);

		// Invariants
	    RealExpression a_geq_zero = -a;   // a >= 0
	    RealExpression a_leq_one = a-1.0;   // a <= 1
	    //valve.new_invariant(modulate,a_geq_zero);
	    //valve.new_invariant(modulate,a_leq_one);

	/// Create the controller automaton

		HybridIOAutomaton controller("controller");

		// Add the input/output variables
		controller.add_input_var(x);
		controller.add_output_var(w);

		// States
		DiscreteLocation open("open");
		DiscreteLocation stabilising("stabilising");
		DiscreteLocation close("close");

		// Create the discrete events
		DiscreteEvent stabilising_after_closing("stabilising_after_closing");
		DiscreteEvent stabilising_after_opening("stabilising_after_opening");
		DiscreteEvent start_closing("start_closing");
		DiscreteEvent start_opening("start_opening");

	    // Create the dynamics
	    RealExpression open_d = 0.0;
	    RealExpression close_d = 0.0;
	    RealExpression stabilising_d = Kp*(ref-x);

	    // Dynamics at the different modes
	    controller.new_mode(open);
	    controller.set_dynamics(open, w, open_d);
	    controller.new_mode(close);
	    controller.set_dynamics(close, w, close_d);
	    controller.new_mode(stabilising);
	    controller.set_dynamics(stabilising, w, stabilising_d);

		// Invariants
	    RealExpression w_geq_zero = -w;   // w >= 0
	    RealExpression w_leq_one = w-1.0;   // w <= 1
	    //controller.new_invariant(close,w_geq_zero);
	    //controller.new_invariant(open,w_leq_one);

	    // Create the guards
	    RealExpression x_lesser_ref_minus_delta = -x-delta+ref; // x <= ref - Delta
	    RealExpression x_greater_ref_minus_delta = x-delta-ref; // x >= ref + Delta
	    RealExpression x_lesser_ref_kp_minus_delta = -x+ref-1.0/Kp-delta; // x <= ref - 1/Kp - Delta
	    RealExpression x_greater_ref_kp_minus_delta = x-ref+1.0/Kp-delta; // x >= ref - 1/Kp + Delta

	    // Transitions
	    controller.new_forced_transition(stabilising_after_closing,close,stabilising,x_lesser_ref_minus_delta);
	    controller.new_forced_transition(stabilising_after_opening,open,stabilising,x_greater_ref_kp_minus_delta);
	    controller.new_forced_transition(start_closing,stabilising,close,x_greater_ref_minus_delta);
	    controller.new_forced_transition(start_opening,stabilising,open,x_lesser_ref_kp_minus_delta);

	/// Composition
	HybridIOAutomaton valve_controller = compose("valve_controller",valve,controller,modulate,open);
	HybridIOAutomaton system = compose("watertank-pr",valve_controller,tank,DiscreteLocation("modulate,open"),flow);

	return system;
}


}

#endif /* WATERTANK_PROPORTIONAL_H_ */
