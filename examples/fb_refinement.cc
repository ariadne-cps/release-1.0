/***************************************************************************
 *            fb_refinement.cc
 *
 *  Copyright  2011  Luca Geretti
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
#include "function.h"
#include "examples.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int verifierVerbosity = 1;
	if (argc > 1)
		verifierVerbosity = atoi(argv[1]);

	HybridAutomaton system("fb");

	DiscreteLocation first("first");
	DiscreteLocation second("second");

	DiscreteEvent first2second("f2s");
	DiscreteEvent second2first("s2f");

	RealConstant u("u",Interval(0.8,1.0));

	RealVariable x("x");
	RealVariable y("y");
	List<RealVariable> varlist;
	varlist.append(x);
	varlist.append(y);

	RealExpression x_first = u;
	RealExpression y_first = u;

	RealExpression x_second = u;
	RealExpression y_second = u;

	List<RealExpression> exprlist;
	exprlist.append(x_first);
	exprlist.append(y_first);
	VectorFunction first_d(exprlist, varlist);
	exprlist[0] = x_second;
	exprlist[1] = y_second;
	VectorFunction second_d(exprlist, varlist);

	RealExpression idx = x;
	RealExpression zero = 0.0;
	RealExpression one = 1.0;
	exprlist[0] = zero;
	exprlist[1] = zero;
	VectorFunction reset_zero(exprlist, varlist);
	exprlist[0] = x+2;
	exprlist[1] = y+2;
	VectorFunction reset_plus_one(exprlist, varlist);

	RealExpression guard_f2s_expr = x+y-3.0;
	ScalarFunction guard_f2s(guard_f2s_expr, varlist);
	RealExpression guard_s2f_expr = x+y-15.0;
	ScalarFunction guard_s2f(guard_s2f_expr, varlist);

	system.new_mode(first,first_d);
	system.new_mode(second,second_d);

	system.new_forced_transition(first2second,first,second,reset_plus_one,guard_f2s);
	system.new_forced_transition(second2first,second,first,reset_zero,guard_s2f);


	// The initial values
	HybridImageSet initial_set;
	initial_set[DiscreteLocation("first")] = Box(2, -0.4,0.4, -0.4,0.4);

	// The domain
	HybridBoxes domain = bounding_boxes(system.state_space(),Box(2,-1.0,11.0,-1.0,11.0));

	// The safety constraint
	RealExpression consexpr_x = x;
	RealExpression consexpr_y = y;
	List<RealExpression> consexpr;
	consexpr.append(consexpr_x);
	consexpr.append(consexpr_y);

	VectorFunction cons_f(consexpr,varlist);
	Box codomain(2,-std::numeric_limits<double>::max(),9.5,-std::numeric_limits<double>::max(), 9.5);
	HybridConstraintSet safety_constraint(system.state_space(),ConstraintSet(cons_f,codomain));

	/// Verification

	Verifier verifier;
	verifier.settings().enable_backward_refinement_for_safety_proving = true;
	verifier.settings().maximum_parameter_depth = 2;
	verifier.settings().plot_results = false;

	SafetyVerificationInput verInput(system, initial_set, domain, safety_constraint);
	cout << verifier.safety(verInput);
}
