/***************************************************************************
 *            monolithic-forced.h
 *
 *  Copyright  2014  Luca Geretti
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

#ifndef MONOLITHIC_FORCED_H_
#define MONOLITHIC_FORCED_H_

#include "ariadne.h"

namespace Ariadne {

HybridAutomaton getSystem()
{
    /// Labeled variables

    // Containing system
	HybridAutomaton system("monolithic-forced");

    // Parameters to be used in the system definition
	RealParameter a("a",0.02);
	RealParameter b("b",0.31);
	RealParameter T("T",4.0);
	RealParameter h("h",6.75);

    // Locations for discrete states
    DiscreteLocation opened("opened");
    DiscreteLocation closed("closed");
    DiscreteLocation opening("opening");
    DiscreteLocation closing("closing");

    // Events for transitions
    DiscreteEvent b_opening("b_opening");
    DiscreteEvent e_opening("e_opening");
    DiscreteEvent b_closing("b_closing");
    DiscreteEvent e_closing("e_closing");

    // State variables
    RealVariable x("x");
    RealVariable y("y");

    /// Dynamics

    // Expressions for the dynamics for x, on every location
    RealExpression x_opening_closing = -a*sqrt(x) + b*y;
    RealExpression x_opened = -a*x + b;
    RealExpression x_closed = -a*x;

    // Expressions for the dynamics for y, on every location
    RealExpression y_opening = 1.0/T;
    RealExpression y_closing = -1.0/T;
    RealExpression y_opened_closed = 0.0;

    // Association of the variables to each expression
    List<RealVariable> varlist;
    varlist.append(x);
    varlist.append(y);
    List<RealExpression> exprlist;
    exprlist.append(x_opened);
    exprlist.append(y_opened_closed);
    VectorFunction dyn_opened(exprlist, varlist);
    exprlist[0] = x_closed;
    VectorFunction dyn_closed(exprlist, varlist);
    exprlist[0] = x_opening_closing;
    exprlist[1] = y_opening;
    VectorFunction dyn_opening(exprlist, varlist);
    exprlist[1] = y_closing;
    VectorFunction dyn_closing(exprlist, varlist);

    // Registration of the dynamics for each location
    system.new_mode(opened,dyn_opened);
    system.new_mode(closing,dyn_closing);
    system.new_mode(closed,dyn_closed);
    system.new_mode(opening,dyn_opening);

    /// Transitions

    // Reset functions
    RealExpression idx = x;
    RealExpression zero = 0.0;
    RealExpression one = 1.0;
    exprlist[0] = idx;
    exprlist[1] = zero;
    VectorFunction reset_y_zero(exprlist, varlist);
    exprlist[1] = one;
    VectorFunction reset_y_one(exprlist, varlist);

    // Guards (where f(x) >= 0 must hold for the guard to be true)
    RealExpression x_leq_min = -x + h;  
    ScalarFunction guard_b_opening(x_leq_min, varlist);
    RealExpression y_geq_one = y - 1.0;
    ScalarFunction guard_e_opening(y_geq_one, varlist);
    RealExpression x_geq_max = x - h;
    ScalarFunction guard_b_closing(x_geq_max, varlist);
    RealExpression y_leq_zero = -y;
    ScalarFunction guard_e_closing(y_leq_zero, varlist);

    // Registration of the transitions
    system.new_forced_transition(b_closing,opened,closing,reset_y_one,guard_b_closing);
    system.new_forced_transition(e_closing,closing,closed,reset_y_zero,guard_e_closing);
    system.new_forced_transition(b_opening,closed,opening,reset_y_zero,guard_b_opening);
    system.new_forced_transition(e_opening,opening,opened,reset_y_one,guard_e_opening);

	return system;
}


}

#endif /* MONOLITHIC_FORCED_H_ */
