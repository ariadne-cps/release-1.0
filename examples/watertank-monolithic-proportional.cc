/***************************************************************************
 *            watertank-monolithic-proportional-verify.cc
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
#include "taylor_calculus.h"
#include "examples.h"

using namespace Ariadne;

int main(int argc,char *argv[]) 
{
	int verifierVerbosity = 1;
	if (argc > 1)
		verifierVerbosity = atoi(argv[1]);

	// The system
	HybridAutomaton system = Ariadne::getWatertankMonolithicProportional();

	// The initial values
	HybridImageSet initial_set;
	initial_set[DiscreteState(1)] = Box(2, 6.75,6.75, 0.0,1.0);
	initial_set[DiscreteState(2)] = Box(2, 6.75,6.75, 0.0,1.0);
	initial_set[DiscreteState(3)] = Box(2, 6.75,6.75, 0.0,1.0);

	// The domain
	HybridBoxes domain = bounding_boxes(system.state_space(),Box(2,-0.1,10.0,-0.1,1.1));

	// The safety constraint
	RealVariable x("x");
	RealVariable y("y");
	List<RealVariable> varlist;
	varlist.append(x);
	varlist.append(y);
	RealExpression expr = x;
	List<RealExpression> consexpr;
	consexpr.append(expr);
	VectorFunction cons_f(consexpr,varlist);
	Box codomain(1,5.25,8.25);
	HybridConstraintSet safety_constraint(system.state_space(),ConstraintSet(cons_f,codomain));

	SafetyVerificationInput verInfo(system, initial_set, domain, safety_constraint);

	/// Verification

	TaylorCalculus outer_integrator(2,2,1e-4);
	TaylorCalculus lower_integrator(4,6,1e-10);
	ImageSetHybridEvolver evolver(outer_integrator,lower_integrator);
	HybridReachabilityAnalyser analyser(evolver);
	analyser.settings().highest_maximum_grid_depth = 7;
	Verifier verifier(analyser);
	verifier.verbosity = verifierVerbosity;
	verifier.settings().plot_results = false;

	/// Analysis parameters
	RealConstantSet parameters;
	parameters.insert(RealConstant("ref",Interval(5.25,8.25)));
	parameters.insert(RealConstant("Kp",Interval(0.2,0.8)));

	std::list<ParametricOutcome> results = verifier.parametric_safety(verInfo, parameters);
	draw(system.name(),results);
}
