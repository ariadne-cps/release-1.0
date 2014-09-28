/***************************************************************************
 *            boost-converter.cc
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
#include "function.h"
#include "boost-converter.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int verbosity = 1;
	if (argc > 1)
		verbosity = atoi(argv[1]);

    bool plot_results = true;

	// The system
	HybridAutomaton system = Ariadne::getBoostConverter();

/*
    HybridEvolver evolver(system);
    evolver.verbosity = verbosity;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(3,2.0);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.001;
    }
    
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("start"),Box(3, 0.0,0.0, 0.0,0.0, 0.0,0.0));
  
    HybridTime evol_limits(18.0,14);
 
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evol_limits,UPPER_SEMANTICS);

    if (plot_results) {
        PlotHelper plotter(system.name());
        plotter.plot(orbit.reach(),"reach");
    }
*/
/*
    int accuracy = 7;

    HybridBoxes domain(system.state_space(),Box(3, -0.1,1.1, 1.0,2.6, 5.0,6.0));

    HybridBoundedConstraintSet initial_set(system.state_space());
    initial_set[DiscreteLocation("incr")] = Box(3, 0.0,0.0, 1.173,1.173, 5.65,5.65);

    HybridReachabilityAnalyser analyser(system,domain,accuracy);
    analyser.verbosity = verbosity;

    HybridDenotableSet outer_reach = analyser.outer_chain_reach(initial_set);

    if (plot_results) {
        PlotHelper plotter(system.name());
        plotter.plot(outer_reach,"outer",accuracy);
    }
*/

    HybridBoxes domain(system.state_space(),Box(3, 0.0,1.0, 0.0,3.5, 4.0,7.5));

    HybridBoundedConstraintSet initial_set(system.state_space());
    initial_set[DiscreteLocation("incr")] = Box(3, 0.0,0.0, 2.0,2.0, 5.5,5.5);

    RealVariable t("t");
    RealVariable iL("iL");
    RealVariable vO("vO");
    List<RealVariable> varlist;
    varlist.append(t);
    varlist.append(iL);
    varlist.append(vO);
    RealExpression expr = vO;
    List<RealExpression> consexpr;
    consexpr.append(expr);
    VectorFunction cons_f(consexpr,varlist);
    Box codomain(1,4.5,6.5);

    HybridConstraintSet safety_constraint(system.state_space(),ConstraintSet(cons_f,codomain));

    Verifier verifier;
    verifier.verbosity = verbosity;
    verifier.ttl = 3600;
    verifier.settings().plot_results = plot_results;
    verifier.settings().enable_backward_refinement_for_safety_proving = false;
    verifier.settings().use_param_midpoints_for_proving = true;

	RealParameterSet parameters;
	parameters.insert(RealParameter("T",Interval(0.5,1.0)));
	parameters.insert(RealParameter("d",Interval(0.1,0.5)));

    SafetyVerificationInput verInput(system, initial_set, domain, safety_constraint);

    std::list<ParametricOutcome> results = verifier.parametric_safety(verInput, parameters);
    draw(system.name(),results);

}
