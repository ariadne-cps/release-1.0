/***************************************************************************
 *            watertank-compositional-hysteresis.cc
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
#include "watertank-compositional-hysteresis.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int verb = 1;
	if (argc > 1)
		verb = atoi(argv[1]);

	HybridIOAutomaton system = Ariadne::getWatertankCompositionalHysteresis();

	// Verification information

	// The initial values
	HybridBoundedConstraintSet initial_set(system.state_space());
	initial_set[DiscreteLocation("flow,idle,rising")] = Box(2, 6.0,7.5, 1.0,1.0);
	initial_set[DiscreteLocation("flow,idle,falling")] = Box(2, 6.0,7.5, 0.0,0.0);

	HybridBoxes domain = bounding_boxes(system.state_space(),Box(2,4.5,9.0,-0.1,1.1));

	// The safety constraint
	List<RealVariable> varlist;
	RealVariable x("x");
	RealVariable y("y");
	varlist.append(x);
	varlist.append(y);
	RealExpression expr = x;
	List<RealExpression> consexpr;
	consexpr.append(expr);
	VectorFunction cons_f(consexpr,varlist);
	Box codomain(1,5.25,8.25);
	HybridConstraintSet safety_constraint(system.state_space(),ConstraintSet(cons_f,codomain));

	/// Verification

	Verifier verifier;
	verifier.verbosity = verb;
	verifier.settings().maximum_parameter_depth = 3;
	verifier.settings().plot_results = true;

	RealParameterSet parameters;
	parameters.insert(RealParameter("hmin",Interval(5.0,6.0)));
	parameters.insert(RealParameter("hmax",Interval(7.5,8.5)));

	SafetyVerificationInput verInfo(system, initial_set, domain, safety_constraint);

	cout << verifier.safety(verInfo);
	//std::list<ParametricOutcome> results = verifier.parametric_safety(verInfo, parameters);
	//draw(system.name(),results);
}
