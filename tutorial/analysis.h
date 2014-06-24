/***************************************************************************
 *            analysis.h
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

#include "ariadne.h"

#ifndef ANALYSIS_H_
#define ANALYSIS_H_

using namespace Ariadne;

void analyse(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results)
{
	// The domain
	HybridBoxes domain(system.state_space(),Box(2,4.5,9.0,0.0,1.0));

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
	Box codomain(1,5.52,8.25);
	HybridConstraintSet safety_constraint(system.state_space(),ConstraintSet(cons_f,codomain));

	// The parameters
	RealParameterSet parameters;
	parameters.insert(RealParameter("hmin",Interval(5.0,6.0)));
	parameters.insert(RealParameter("hmax",Interval(7.5,8.5)));

	/// Verification

	Verifier verifier;
	verifier.verbosity = verbosity;
	verifier.ttl = 5;
	verifier.settings().plot_results = plot_results;

	SafetyVerificationInput verInput(system, initial_set, domain, safety_constraint);

	std::list<ParametricOutcome> results = verifier.parametric_safety(verInput, parameters);
	draw(system.name(),results);
}

#endif /* ANALYSIS_H_ */
