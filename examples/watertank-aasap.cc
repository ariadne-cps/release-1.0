/***************************************************************************
 *            watertank-aasap.cc
 *
 *  Copyright  2010  Luca Geretti
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

#include "ariadne.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int verifierVerbosity = 1;
	if (argc > 1)
		verifierVerbosity = atoi(argv[1]);

    /// Set the system parameters
	RealConstant a("a",0.02);
	RealConstant b("b",0.3);
	RealConstant Ta("Ta",4.0);
	RealConstant T("T",0.1);
	RealConstant hmin("hmin",5.5); 
	RealConstant hmax("hmax",8.0);

    // System variables
    RealVariable x("x");    // water level
    RealVariable y("y");    // valve aperture
    RealVariable t_out("t_out");    // valve aperture

	// The parameter to modify, its interval and the tolerance
    RealConstant Delta("Delta",Interval(0.0,0.0));
    RealConstant parameter = Delta;

    // Create the tank automaton

	    HybridIOAutomaton tank("tank");

		// States
		DiscreteState flow("flow");  

		// Add the input/output variables
    	tank.add_input_var(y);
    	tank.add_output_var(x);
		
		// Only one state with no transitions and no invariants
		RealExpression dyn = - a * x + b * y;
		tank.new_mode(flow);
		tank.set_dynamics(flow, x, dyn);

	// Create the valve automaton

		HybridIOAutomaton valve("valve");

		// States
		DiscreteState idle("idle");
		DiscreteState opening("opening");
		DiscreteState closing("closing");
		
		// The valve has one output var (the valve aperture)
		valve.add_output_var(y);

		// Two input events (open and close) and one internal event
		DiscreteEvent e_open("OPEN");
		valve.add_input_event(e_open);
		DiscreteEvent e_close("CLOSE");
		valve.add_input_event(e_close);
		DiscreteEvent e_idle("idle");
		valve.add_internal_event(e_idle);

		// Three states:
		// Idle (valve either fully closed or fully opened)
		RealExpression dynidle = 0.0;
		valve.new_mode(idle);
		//valve.new_invariant(idle, -y);
		//valve.new_invariant(idle, y-1.0);
		valve.set_dynamics(idle, y, dynidle);
		// Opening (valve is opening)
		valve.new_mode(opening);
		//valve.new_invariant(opening, -y);
		//valve.new_invariant(opening, y-1.0);
		RealExpression dynopening = 1.0/Ta;
		valve.set_dynamics(opening, y, dynopening);
		// Closing (valve is closing)
		valve.new_mode(closing);
		//valve.new_invariant(closing, -y);
		//valve.new_invariant(closing, y-1.0);
		RealExpression dynclosing = -1.0/Ta;
		valve.set_dynamics(closing, y, dynclosing);
		
		// Transitions

		// the identity y' = y.
		std::map< RealVariable, RealExpression> reset_y_identity;
		reset_y_identity[y] = y;
		std::map< RealVariable, RealExpression> reset_y_one;
		reset_y_one[y] = 1.0;
		std::map< RealVariable, RealExpression> reset_y_zero;
		reset_y_zero[y] = 0.0;

		// when open is received, go to opening
		valve.new_unforced_transition(e_open, idle, opening, reset_y_identity);
		valve.new_unforced_transition(e_open, opening, opening, reset_y_identity);
		valve.new_unforced_transition(e_open, closing, opening, reset_y_identity);
		 // when closed is received, go to closing
		valve.new_unforced_transition(e_close, idle, closing, reset_y_identity);
		valve.new_unforced_transition(e_close, opening, closing, reset_y_identity);
		valve.new_unforced_transition(e_close, closing, closing, reset_y_identity);
		// when the valve is fully opened go from opening to idle
		RealExpression y_geq_one = y - 1.0;
		valve.new_forced_transition(e_idle, opening, idle, reset_y_identity, y_geq_one);
		// when the valve is fully closed go from closing to idle
		RealExpression y_leq_zero = - y;
		valve.new_forced_transition(e_idle, closing, idle, reset_y_identity, y_leq_zero);

	// Create the evaluator automaton

		    HybridIOAutomaton evaluator("evaluator");

			// States
			DiscreteState deep("deep");
			DiscreteState shallow("shallow");

			// The evaluator has two output events
			DiscreteEvent e_high("HIGH");
			evaluator.add_output_event(e_high);
			DiscreteEvent e_low("LOW");
			evaluator.add_output_event(e_low);

			// An high event has been issued
			evaluator.new_mode(deep);
			// A low event has been issued;
			evaluator.new_mode(shallow);

			// Transitions
			// When the water is greater than hmax, send a high event
			RealExpression x_geq_hmax = x - hmax;
			evaluator.new_forced_transition(e_high, shallow, deep, x_geq_hmax);
			// When the water is lower than hmin, send a open command
			RealExpression x_leq_hmin = hmin - x;
			evaluator.new_forced_transition(e_low, deep, shallow, x_leq_hmin);


	// Create the controller automaton

	    HybridIOAutomaton controller("controller");

		// States
		DiscreteState increase("increase");
		DiscreteState decrease("decrease");
		DiscreteState nothing("nothing");

		// Involved variables
		controller.add_internal_var(t_out);
 
		// Two input events (high and low)
		controller.add_input_event(e_high);
		controller.add_input_event(e_low);
		// Two output events (open and close)
		controller.add_output_event(e_open); 
		controller.add_output_event(e_close);
		
		// Two states:
		// The controller is about to increase the level by opening the valve
		controller.new_mode(increase);
		controller.set_dynamics(increase, t_out, 1.0);
		// The controller is about to decrease the level by closing the valve
		controller.new_mode(decrease);
		controller.set_dynamics(decrease, t_out, 1.0);
		// The controller does nothing
		controller.new_mode(nothing);
		controller.set_dynamics(nothing, t_out, 0.0);

		// Resets

			// Reset t_out to zero
			std::map< RealVariable, RealExpression> reset_t_out_zero;
			reset_t_out_zero[t_out] = 0.0;

		// Transitions

		// When the water level is high, wait before issuing the close command
		controller.new_unforced_transition(e_high, nothing, decrease, reset_t_out_zero);
		// When the water level is low, wait before issuing the open command
		controller.new_unforced_transition(e_low, nothing, increase, reset_t_out_zero);
		// When the clock period has expired, issue the open command
		controller.new_forced_transition(e_open, increase, nothing, reset_t_out_zero, t_out-T);
		// When the clock period has expired, issue the close command
		controller.new_forced_transition(e_close, decrease, nothing, reset_t_out_zero, t_out-T);

	/// Compose the automata
	HybridIOAutomaton tank_valve = compose("tank,valve",tank,valve,flow,idle);
	HybridIOAutomaton tank_valve_evaluator = compose("tank,valve,evaluator",tank_valve,evaluator,DiscreteState("flow,idle"),shallow);
	HybridIOAutomaton system_io = compose("watertank-aasap",tank_valve_evaluator,controller,DiscreteState("flow,idle,shallow"),nothing);

	/// Create the monolithic automaton
	HybridAutomaton system;
	RealSpace space;
	make_lpair<HybridAutomaton,RealSpace>(system,space) = make_monolithic_automaton(system_io);

	/// Add access to the Delta constant
	//system.register_accessible_constant(Delta);

	// Verification information

	// The initial values
	HybridImageSet initial_set;
	initial_set[DiscreteState("flow,idle,shallow,nothing")] = Box(3, 0.0,0.0, 6.0,6.0, 1.0,1.0);

	// The safety constraint
	List<RealVariable> varlist;
	varlist.append(x);
	varlist.append(y);
	varlist.append(t_out);
	RealExpression expr = x;
	List<RealExpression> consexpr;
	consexpr.append(expr);
	VectorFunction cons_f(consexpr,varlist);
	Box codomain(1,5.25,8.25);
	HybridConstraintSet safety_constraint(system.state_space(),ConstraintSet(cons_f,codomain));

	HybridBoxes domain = bounding_boxes(system.state_space(),Box(3,-0.1,0.2,4.0,10.0,-0.5,1.5));

	/// Verification

	// Create an evolver and analyser objects, then set their verbosity
	HybridEvolver evolver;
	HybridReachabilityAnalyser analyser(evolver);
	analyser.settings().lowest_maximum_grid_depth = 0;
	analyser.settings().highest_maximum_grid_depth = 8;
	analyser.free_cores = 0;
	Verifier verifier(analyser);
	verifier.verbosity = verifierVerbosity;

	SafetyVerificationInput verInfo(system, initial_set, domain, safety_constraint);
    verifier.safety(verInfo);

}
