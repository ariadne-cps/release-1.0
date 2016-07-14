/*****************************************************************************************************
 *            composition.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 * Example of composition of two automata.
 *
 *****************************************************************************************************/

#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Create a HybridAutomaton object
    HybridIOAutomaton a2("a2");

    /// Create the discrete states
    DiscreteLocation a2_l1("a2_l1");
    DiscreteLocation a2_l2("a2_l2");

    RealVariable y("y");

    DiscreteEvent a2_e12("a2_e12");
    DiscreteEvent a2_e21("a2_e21");

    a2.add_output_event(a2_e21);
    a2.add_output_event(a2_e12);

    a2.add_output_var(y);

	a2.new_mode(a2_l1);
	a2.new_mode(a2_l2);

	RealExpression dyn_y1 = 1.0;
	RealExpression dyn_y2 = -1.0;
	a2.set_dynamics(a2_l1, y, dyn_y1);
	a2.set_dynamics(a2_l2, y, dyn_y2);

	RealExpression a2_g21 = y-1.0;
	RealExpression a2_g12 = -y-1.0;

	a2.new_forced_transition(a2_e12,a2_l1,a2_l2,a2_g12);
	a2.new_forced_transition(a2_e21,a2_l2,a2_l1,a2_g21);


    /// Create a HybridAutomaton object
    HybridIOAutomaton a1("a1");

    /// Create the discrete states
    DiscreteLocation a1_l1("a1_l1");
    DiscreteLocation a1_l2("a1_l2");

    RealVariable x("x");

    DiscreteEvent a1_e12("a1_e12");
    DiscreteEvent a1_e21("a1_e21");

    a1.add_input_event(a2_e21);
    a1.add_output_event(a1_e12);
    a1.add_output_event(a1_e21);

    a1.add_output_var(x);

	a1.new_mode(a1_l1);
	a1.new_mode(a1_l2);

	RealExpression dyn_x1 = 1.0;
	RealExpression dyn_x2 = -1.0;
	a1.set_dynamics(a1_l1, x, dyn_x1);
	a1.set_dynamics(a1_l2, x, dyn_x2);

	RealExpression a1_g12 = x-1.0;
	RealExpression a1_g21 = -1.0-x;

	a1.new_unforced_transition(a2_e21,a1_l2,a1_l1);
	a1.new_forced_transition(a1_e12,a1_l1,a1_l2,a1_g12);
	a1.new_forced_transition(a1_e21,a1_l2,a1_l1,a1_g21);

	HybridIOAutomaton system = compose("system",a1,a2,a1_l1,a2_l1);

	cout << system << endl;
}
