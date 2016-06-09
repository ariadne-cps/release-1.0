/*****************************************************************************************************
 *            laser.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides a circle with proper resets to keep the trajectory stable.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef CIRCLE_H_
#define CIRCLE_H_

namespace Ariadne {

HybridIOAutomaton getCircle()
{
    /// Parameters
	RealParameter R("R",Interval(3.0,3.0));
	RealParameter w("w",4.0);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("circle-stable");

    /// Create the discrete states
    DiscreteLocation first("1");
    DiscreteLocation second("2");
    DiscreteLocation third("3");
    DiscreteLocation fourth("4");

    RealVariable x("x");
    RealVariable y("y");

    automaton.add_output_var(x);
    automaton.add_output_var(y);

    // Events
    DiscreteEvent first2second("12");
    DiscreteEvent second2third("23");
    DiscreteEvent third2fourth("34");
    DiscreteEvent fourth2first("41");

	automaton.new_mode(first);
	automaton.new_mode(second);
	automaton.new_mode(third);
	automaton.new_mode(fourth);

	RealExpression dyn_x = - 2.0*Ariadne::pi<Real>()*w * y;
	automaton.set_dynamics(first, x, dyn_x);
	automaton.set_dynamics(second, x, dyn_x);
	automaton.set_dynamics(third, x, dyn_x);
	automaton.set_dynamics(fourth, x, dyn_x);
	RealExpression dyn_y = 2.0*Ariadne::pi<Real>()*w * x;
	automaton.set_dynamics(first, y, dyn_y);
	automaton.set_dynamics(second, y, dyn_y);
	automaton.set_dynamics(third, y, dyn_y);
	automaton.set_dynamics(fourth, y, dyn_y);

	/// Transitions
	// Guards
	RealExpression x_greater_zero = x;
	RealExpression y_greater_zero = y;
	RealExpression x_lesser_zero = -x;
	RealExpression y_lesser_zero = -y;
	// Resets
	std::map<RealVariable,RealExpression> reset12;
	reset12[x] = 0.0;
	reset12[y] = R;
	std::map<RealVariable,RealExpression> reset23;
	reset23[x] = -R;
	reset23[y] = 0.0;
	std::map<RealVariable,RealExpression> reset34;
	reset34[x] = 0.0;
	reset34[y] = -R;
	std::map<RealVariable,RealExpression> reset41;
	reset41[x] = R;
	reset41[y] = 0.0;
	automaton.new_forced_transition(first2second,first,second,reset12,x_lesser_zero);
	automaton.new_forced_transition(second2third,second,third,reset23,y_lesser_zero);
	automaton.new_forced_transition(third2fourth,third,fourth,reset34,x_greater_zero);
	automaton.new_forced_transition(fourth2first,fourth,first,reset41,y_greater_zero);

	return automaton;
}

}

#endif
