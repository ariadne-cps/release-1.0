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

HybridAutomaton getBoostConverter()
{
	HybridAutomaton system("boost");

    /// Set the system parameters
	RealParameter Vi("Vi",3.3);
	RealParameter L("L",1.0);
	RealParameter C("C",1.0);
	RealParameter R("R",10.0);
	RealParameter T("T",1.0);
	RealParameter d("d",0.4);

    /// Create three discrete states
    DiscreteLocation incr("incr");
    DiscreteLocation decr("decr");
    DiscreteLocation zero("zero");

    /// Create the discrete events
    DiscreteEvent turn_on_from_zero("turn_on_from_zero");
    DiscreteEvent turn_on_from_decr("turn_on_from_decr");
    DiscreteEvent turn_off("turn_off");
    DiscreteEvent current_becomes_zero("current_is_zero");

    // System variables
    RealVariable t("t");    // Clock
    RealVariable iL("iL");    // Inductor current
    RealVariable vO("vO");    // Output voltage
    List<RealVariable> varlist;
    varlist.append(t);
    varlist.append(iL);
    varlist.append(vO);

    // t dynamics
    RealExpression t_any = 1.0;

    // iL dynamics
    RealExpression iL_incr = Vi/L;
    RealExpression iL_decr = (Vi-vO)/L;
    RealExpression iL_zero = 0.0;

    // vO dynamics
    RealExpression vO_incr = -vO/(R*C);
    RealExpression vO_decr = iL/C - vO/(R*C);
    RealExpression vO_zero = -vO/(R*C);

    // Dynamics at the different modes
    List<RealExpression> exprlist;
    exprlist.append(t_any);
    exprlist.append(iL_incr);
    exprlist.append(vO_incr);
    VectorFunction dyn_incr(exprlist, varlist);
    exprlist[1] = iL_decr;
    exprlist[2] = vO_decr;
    VectorFunction dyn_decr(exprlist, varlist);
    exprlist[1] = iL_zero;
    exprlist[2] = vO_zero;
    VectorFunction dyn_zero(exprlist, varlist);

    // Reset functions
    RealExpression id_t_r = t;
    RealExpression id_iL_r = iL;
    RealExpression id_vO_r = vO;
    RealExpression zero_t_r = 0.0;
    exprlist[0] = id_t_r;
    exprlist[1] = id_iL_r;
    exprlist[2] = id_vO_r;
    VectorFunction reset_identity(exprlist, varlist);
    exprlist[0] = zero_t_r;
    VectorFunction reset_t_zero(exprlist, varlist);

    // Create the guards.
    // Guards are true when f(x) >= 0
    RealExpression t_geq_dT = t - d*T;       // t >= d*T
    ScalarFunction turn_off_g(t_geq_dT, varlist);
    RealExpression iL_leq_zero = -iL;                 // iL <= 0    
    ScalarFunction zero_current_g(iL_leq_zero, varlist);
    RealExpression t_geq_T = t - T;       // t >= T
    ScalarFunction turn_on_g(t_geq_T, varlist);

    /// Build the automaton
    system.new_mode(incr,dyn_incr);
    system.new_mode(decr,dyn_decr);
    system.new_mode(zero,dyn_zero);

    system.new_forced_transition(turn_off,incr,decr,reset_identity,turn_off_g);
    system.new_forced_transition(current_becomes_zero,decr,zero,reset_identity,zero_current_g);
    system.new_forced_transition(turn_on_from_zero,zero,incr,reset_t_zero,turn_on_g);
    system.new_forced_transition(turn_on_from_decr,decr,incr,reset_t_zero,turn_on_g);

	return system;
}


}

#endif /* BOOST_CONVERTER_H_ */
