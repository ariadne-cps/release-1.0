/***************************************************************************
 *            watertank-monolithic-hysteresis.cc
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
#include "function.h"
#include "examples.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int verb = 1;
	if (argc > 1)
		verb = atoi(argv[1]);

	// The system
	HybridAutomaton system = Ariadne::getWatertankMonolithicHysteresis();

    HybridEvolver evolver(system);
    evolver.verbosity = verb;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(2,3.0);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.5;
    }

    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("opened"),Box(2, 6.5,6.5, 1.0,1.0));

    HybridTime evol_limits(80.0,5);

    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evol_limits,UPPER_SEMANTICS);

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");

	/*
	// The initial values
	HybridBoundedConstraintSet initial_set(system.state_space());
	initial_set[DiscreteLocation("opened")] = Box(2, 6.0,7.5, 1.0,1.0);
	initial_set[DiscreteLocation("closed")] = Box(2, 6.0,7.5, 0.0,0.0);

	// The domain
	HybridBoxes domain(system.state_space(),Box(2,4.5,9.0,0.0,1.0));

    int accuracy = 4;
    HybridTime evol_limits(160.0,16);

    HybridReachabilityAnalyser analyser(system,domain,accuracy);
    analyser.verbosity = verb;

    bool plot_results = true;

    HybridDenotableSet upper_reach = analyser.upper_reach(initial_set,evol_limits);

    if (plot_results) {
        PlotHelper plotter(system.name());
        plotter.plot(upper_reach,"upper",accuracy);
    }
	*/
	/*

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
	verifier.verbosity = verb;
	verifier.ttl = 60;
	verifier.settings().plot_results = true;

	SafetyVerificationInput verInput(system, initial_set, domain, safety_constraint);

	verifier.safety(verInput);
	//std::list<ParametricOutcome> results = verifier.parametric_safety(verInput, parameters);
	//draw(system.name(),results);

	*/
}
