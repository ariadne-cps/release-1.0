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
	HybridIOAutomaton system("boost");

    /// Set the system parameters
	RealParameter Vi("Vi",3.3);
	RealParameter R("R",5.0);
	RealParameter C("C",0.001);
	RealParameter L("L",0.0011);
	RealParameter d("d",0.45);
	RealParameter T("T",0.001);

    /// Create three discrete states
    DiscreteLocation incr("incr");
    DiscreteLocation decr("decr");
    DiscreteLocation zero("zero");

    /// Create the discrete events
    DiscreteEvent turn_on_from_zero("turn_on_from_zero");
    DiscreteEvent turn_on_from_decr("turn_on_from_decr");
    DiscreteEvent turn_off("turn_off");
    DiscreteEvent current_becomes_zero("current_is_zero");

    system.add_output_event(turn_on_from_zero);
    system.add_output_event(turn_on_from_decr);
    system.add_output_event(turn_off);
    system.add_output_event(current_becomes_zero);

    // System variables
    RealVariable clk("clk");    // Clock
    RealVariable iL("iL");    // Inductor current
    RealVariable vO("vO");    // Output voltage

    system.add_output_var(clk);
    system.add_output_var(iL);
    system.add_output_var(vO);

    // clk dynamics
    RealExpression clk_any = 1.0;
    // iL dynamics
    RealExpression iL_incr = Vi/L;
    RealExpression iL_decr = (Vi-vO)/L;
    RealExpression iL_zero = 0.0;
    // vO dynamics
    RealExpression vO_incr = -vO/(R*C);
    RealExpression vO_decr = iL/C - vO/(R*C);
    RealExpression vO_zero = -vO/(R*C);

    system.new_mode(incr);
    system.set_dynamics(incr,clk,clk_any);
    system.set_dynamics(incr,iL,iL_incr);
    system.set_dynamics(incr,vO,vO_incr);

    system.new_mode(decr);
    system.set_dynamics(decr,clk,clk_any);
    system.set_dynamics(decr,iL,iL_decr);
    system.set_dynamics(decr,vO,vO_decr);

    system.new_mode(zero);
    system.set_dynamics(zero,clk,clk_any);
    system.set_dynamics(zero,iL,iL_zero);
    system.set_dynamics(zero,vO,vO_zero);

    // Reset functions
    std::map< RealVariable, RealExpression> reset_identity;
    reset_identity[clk] = clk;
    reset_identity[iL] = iL;
    reset_identity[vO] = vO;

    std::map< RealVariable, RealExpression> reset_clk_zero;
    reset_clk_zero[clk] = 0.0;
    reset_clk_zero[iL] = iL;
    reset_clk_zero[vO] = vO;

    // Guards
    RealExpression clk_geq_dT = clk - d*T;    // clk >= d*T
    RealExpression iL_leq_zero = -iL;     // iL <= 0
    RealExpression clk_geq_T = clk - T;       // clk >= T

    system.new_forced_transition(turn_off,incr,decr,reset_identity,clk_geq_dT);
    system.new_forced_transition(current_becomes_zero,decr,zero,reset_identity,iL_leq_zero);
    system.new_forced_transition(turn_on_from_zero,zero,incr,reset_clk_zero,clk_geq_T);
    system.new_forced_transition(turn_on_from_decr,decr,incr,reset_clk_zero,clk_geq_T);

	return system;
}


}

#endif /* BOOST_CONVERTER_H_ */
