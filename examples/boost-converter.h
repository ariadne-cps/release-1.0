/***************************************************************************
 *            boost-converter.h
 *
 *  Copyright  2014   Luca Geretti
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

#ifndef BOOST_CONVERTER_H_
#define BOOST_CONVERTER_H_

#include "ariadne.h"

namespace Ariadne {

HybridIOAutomaton getBoostConverter()
{
    /// Discrete events
    DiscreteEvent turn_on("turn_on");
    DiscreteEvent turn_off("turn_off");
    DiscreteEvent current_becomes_zero("current_becomes_zero");

    // Variables
    RealVariable clk("clk");    // Clock
    RealVariable iL("iL");    // Inductor current
    RealVariable vO("vO");    // Output voltage

	/// Converter automaton
	HybridIOAutomaton converter("converter");

    /// Parameters
	RealParameter Vi("Vi",3.3);
	RealParameter R("R",5.0);
	RealParameter C("C",0.001);
	RealParameter L("L",0.0011);

    /// States
    DiscreteLocation incr("incr");
    DiscreteLocation decr("decr");
    DiscreteLocation zero("zero");

    /// Register events
    converter.add_input_event(turn_on);
    converter.add_input_event(turn_off);
    converter.add_internal_event(current_becomes_zero);

    /// Variables
    converter.add_internal_var(iL);
    converter.add_output_var(vO);

    /// Set dynamics

    // iL dynamics
    RealExpression iL_incr = Vi/L;
    RealExpression iL_decr = (Vi-vO)/L;
    RealExpression iL_zero = 0.0;
    // vO dynamics
    RealExpression vO_incr = -vO/(R*C);
    RealExpression vO_decr = iL/C - vO/(R*C);
    RealExpression vO_zero = -vO/(R*C);

    converter.new_mode(incr);
    converter.set_dynamics(incr,iL,iL_incr);
    converter.set_dynamics(incr,vO,vO_incr);

    converter.new_mode(decr);
    converter.set_dynamics(decr,iL,iL_decr);
    converter.set_dynamics(decr,vO,vO_decr);

    converter.new_mode(zero);
    converter.set_dynamics(zero,iL,iL_zero);
    converter.set_dynamics(zero,vO,vO_zero);

    /// Transitions

    // Guards
    RealExpression iL_leq_zero = -iL;     // iL <= 0

    converter.new_unforced_transition(turn_off,incr,decr);
    converter.new_unforced_transition(turn_on,zero,incr);
    converter.new_unforced_transition(turn_on,decr,incr);
    converter.new_forced_transition(current_becomes_zero,decr,zero,iL_leq_zero);

    /// Controller automaton
	HybridIOAutomaton controller("controller");

    /// Parameters
	RealParameter d("d",0.45);
	RealParameter T("T",0.001);

	/// States
	DiscreteLocation below_duty("below_duty");
	DiscreteLocation over_duty("over_duty");

	/// Register events
    controller.add_output_event(turn_on);
    controller.add_output_event(turn_off);

    /// Variables
    controller.add_output_var(clk);

    /// Set dynamics
    RealExpression clk_any = 1.0;

    controller.new_mode(below_duty);
    controller.set_dynamics(below_duty,clk,clk_any);

    controller.new_mode(over_duty);
    controller.set_dynamics(over_duty,clk,clk_any);

    /// Set transitions

    // Resets
    std::map< RealVariable, RealExpression> reset_clk_zero;
    reset_clk_zero[clk] = 0.0;

    // Guards
    RealExpression clk_geq_dT = clk - d*T;    // clk >= d*T
    RealExpression clk_geq_T = clk - T;       // clk >= T

    controller.new_forced_transition(turn_off,below_duty,over_duty,clk_geq_dT);
    controller.new_forced_transition(turn_on,over_duty,below_duty,reset_clk_zero,clk_geq_T);

    /// Composition of automata
    HybridIOAutomaton system = compose("boost",converter,controller,incr,below_duty);

	return system;
}


}

#endif /* BOOST_CONVERTER_H_ */
