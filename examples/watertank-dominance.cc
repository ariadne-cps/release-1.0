/***************************************************************************
 *            watertank-monolithic-dominance.cc
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
#include "watertank-proportional.h"
#include "watertank-hysteresis.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int verb = 1;
	if (argc > 1)
		verb = atoi(argv[1]);

	// The systems
	HybridIOAutomaton system_hy = Ariadne::getWatertankHysteresis();
	HybridIOAutomaton system_pr = Ariadne::getWatertankProportional();

	// The initial values
	HybridBoundedConstraintSet initial_hy(system_hy.state_space());
	initial_hy[DiscreteLocation("flow,idle,rising")] = Box(2, 1.0,1.0, 6.0,7.5);
	initial_hy[DiscreteLocation("flow,idle,falling")] = Box(2, 0.0,0.0, 6.0,7.5);
	HybridBoundedConstraintSet initial_pr(system_pr.state_space());
	initial_pr[DiscreteLocation("flow,stabilising")] = Box(2, 0.0,1.0, 6.75,6.75);
	initial_pr[DiscreteLocation("flow,opening")] = Box(2, 0.0,1.0, 6.75,6.75);
	initial_pr[DiscreteLocation("flow,closing")] = Box(2, 0.0,1.0, 6.75,6.75);

	// The domains
	HybridBoxes domain_hy(system_hy.state_space(),Box(2,0.0,1.0,4.5,9.0));
	HybridBoxes domain_pr(system_pr.state_space(),Box(2,0.0,1.0,2.0,10.0));

	// The projections
	Vector<uint> projection_hy(1,0);
	Vector<uint> projection_pr(1,0);

	// Construct the bundles
	DominanceVerificationInput hysteresis(system_hy,initial_hy,domain_hy,projection_hy);
	DominanceVerificationInput proportional(system_pr,initial_pr,domain_pr,projection_pr);

	Verifier verifier;
	verifier.verbosity = verb;
	verifier.ttl = 300;
	verifier.settings().maximum_parameter_depth = 5;
	verifier.settings().plot_results = false;

	// The parametric dominance parameters
	RealParameterSet parameters;
	parameters.insert(RealParameter("Kp",Interval(0.2,0.8)));
	parameters.insert(RealParameter("ref",Interval(5.25,8.25)));

	std::list<ParametricOutcome> results = verifier.parametric_dominance(proportional, hysteresis, parameters);
	draw("watertank-dominance",results);
}
