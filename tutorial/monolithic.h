/***************************************************************************
 *            monolithic-unforced.h
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

#ifndef MONOLITHIC_UNFORCED_H_
#define MONOLITHIC_UNFORCED_H_

#include "ariadne.h"

namespace Ariadne {

HybridAutomaton getSystem()
{
    /// Labeled variables

    // Containing system
    HybridAutomaton system("monolithic-unforced");

    // Parameters to be used in the system definition
    RealParameter a("a",0.02);
    RealParameter b("b",Interval(0.3,0.32863));
    RealParameter T("T",4.0);
    RealParameter hmin("hmin",5.75);
    RealParameter hmax("hmax",7.75); 
    RealParameter Delta("Delta",0.1);

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
    RealExpression x_opened = -a*sqrt(x) + b;
    RealExpression x_closed = -a*sqrt(x);

    // Expressions for the dynamics for y, on every location
    RealExpression y_opening = 1.0/T;
    RealExpression y_closing = -1.0/T;
    RealExpression y_opened_closed = 0.0;

    // Association of the variables to each expression
    List<RealVariable> varlist;
    varlist.append(x);
    varlist.append(y);
    List<RealExpression> opened_exprlist;
    opened_exprlist.append(x_opened);
    opened_exprlist.append(y_opened_closed);
    VectorFunction dyn_opened(opened_exprlist, varlist);
    List<RealExpression> closed_exprlist;
    closed_exprlist.append(x_closed);
    closed_exprlist.append(y_opened_closed);
    VectorFunction dyn_closed(closed_exprlist, varlist);
    List<RealExpression> opening_exprlist;
    opening_exprlist.append(x_opening_closing);
    opening_exprlist.append(y_opening);
    VectorFunction dyn_opening(opening_exprlist, varlist);
    List<RealExpression> closing_exprlist;
    closing_exprlist.append(x_opening_closing);
    closing_exprlist.append(y_closing);
    VectorFunction dyn_closing(closing_exprlist, varlist);

    // Registration of the dynamics for each location
    system.new_mode(opened,dyn_opened);
    system.new_mode(closing,dyn_closing);
    system.new_mode(closed,dyn_closed);
    system.new_mode(opening,dyn_opening);

    /// Invariants

    // Expressions (where f(x) <= 0 must hold for the invariant to be true)
    RealExpression x_leq_hmax = x - hmax - Delta;
    ScalarFunction inv_opened(x_leq_hmax, varlist);
    RealExpression x_geq_hmin = -x + hmin - Delta;
    ScalarFunction inv_closed(x_geq_hmin, varlist);

    // Registration of the invariants for each location
    system.new_invariant(opened,inv_opened);
    system.new_invariant(closed,inv_closed);

    /// Transitions

    // Reset functions
    RealExpression idx = x;
    RealExpression zero = 0.0;
    RealExpression one = 1.0;
    List<RealExpression> zero_exprlist;
    zero_exprlist.append(idx);
    zero_exprlist.append(zero);
    VectorFunction reset_y_zero(zero_exprlist, varlist);
    List<RealExpression> one_exprlist;
    one_exprlist.append(idx);
    one_exprlist.append(one);
    VectorFunction reset_y_one(one_exprlist, varlist);

    // Guards (where f(x) >= 0 must hold for the guard to be true)
    RealExpression x_leq_hmin = -x + hmin + Delta;
    ScalarFunction guard_b_opening(x_leq_hmin, varlist);
    RealExpression y_geq_one = y - 1.0;
    ScalarFunction guard_e_opening(y_geq_one, varlist);
    RealExpression x_geq_hmax = x - hmax + Delta;
    ScalarFunction guard_b_closing(x_geq_hmax, varlist);
    RealExpression y_leq_zero = -y;
    ScalarFunction guard_e_closing(y_leq_zero, varlist);

    // Registration of the transitions
    system.new_unforced_transition(b_closing,opened,closing,reset_y_one,guard_b_closing);
    system.new_forced_transition(e_closing,closing,closed,reset_y_zero,guard_e_closing);
    system.new_unforced_transition(b_opening,closed,opening,reset_y_zero,guard_b_opening);
    system.new_forced_transition(e_opening,opening,opened,reset_y_one,guard_e_opening);

    return system;
}


}

#endif /* MONOLITHIC_UNFORCED_H_ */
