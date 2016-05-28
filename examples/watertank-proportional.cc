/***************************************************************************
 *            watertank-proportional.cc
 *
 *  Copyright  2016  Luca Geretti
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
#include "watertank-compositional-proportional.h"

using namespace Ariadne;

int main(int argc,char *argv[]) 
{
	int verb = 1;
	if (argc > 1)
		verb = atoi(argv[1]);

	// The system
	HybridIOAutomaton system = Ariadne::getWatertankProportional();

    HybridEvolver evolver(system);
    evolver.verbosity = verb;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(2,0.004);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.0005;
    }

    cout << system << endl;

    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("stabilising,flow"),Box(2, 0.6,0.6, 6.72,6.72));

    HybridTime evol_limits(80.0,5);

    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evol_limits,UPPER_SEMANTICS);

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");

    /*
	// The initial values
	HybridBoundedConstraintSet initial_set(system.state_space());
	initial_set[DiscreteLocation(1)] = Box(2, 6.75,6.75, 0.0,1.0);
	initial_set[DiscreteLocation(2)] = Box(2, 6.75,6.75, 0.0,1.0);
	initial_set[DiscreteLocation(3)] = Box(2, 6.75,6.75, 0.0,1.0);

	// The domain
	HybridBoxes domain(system.state_space(),Box(2,0.0,10.0,0.0,1.0));

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

	Verifier verifier;
	verifier.verbosity = verb;
	verifier.ttl = 60;
	verifier.settings().plot_results = false;

	/// Analysis parameters
	RealParameterSet parameters;
	parameters.insert(RealParameter("ref",Interval(5.25,8.25)));
	parameters.insert(RealParameter("Kp",Interval(0.2,0.8)));

	std::list<ParametricOutcome> results = verifier.parametric_safety(verInfo, parameters);
	draw(system.name(),results);
	*/
}
