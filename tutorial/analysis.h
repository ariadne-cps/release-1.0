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
using namespace std;

// Forward declarations
void finite_time_upper_evolution(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results);
void infinite_time_outer_evolution(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results);
void infinite_time_lower_evolution(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results);
void safety_verification(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results);
void parametric_safety_verification(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results);
HybridConstraintSet getSafetyConstraint(HybridAutomatonInterface& system);

// The main method for the analysis of the system
void analyse(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results)
{
    cout << "1/5: Finite time (upper) evolution... " << endl << flush;
    finite_time_upper_evolution(system,initial_set,verbosity,plot_results);
    cout << "2/5: Infinite time outer evolution... " << endl << flush; 
    infinite_time_outer_evolution(system,initial_set,verbosity,plot_results);
    cout << "3/5: Infinite time lower evolution... " << endl << flush; 
    infinite_time_lower_evolution(system,initial_set,verbosity,plot_results);
    cout << "4/5: Safety verification... " << endl << flush;
    safety_verification(system,initial_set,verbosity,plot_results);
    cout << "5/5: Parametric safety verification... " << endl << flush;
    parametric_safety_verification(system,initial_set,verbosity,plot_results);
}

void finite_time_upper_evolution(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results) {

    int accuracy = 5;

    HybridTime time(40.0,9);

    HybridBoxes domain(system.state_space(),Box(2,4.5,9.0,0.0,1.0));

    HybridReachabilityAnalyser analyser(system,domain,accuracy);
    analyser.verbosity = verbosity;

    HybridDenotableSet upper_reach = analyser.upper_reach(initial_set,time);

    if (plot_results) {
        SystemPlotter plotter(system);
        plotter.plot_set(upper_reach,"upper",accuracy);
    }
}

void infinite_time_outer_evolution(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results) {

    int accuracy = 5;

    HybridBoxes domain(system.state_space(),Box(2,4.5,9.0,0.0,1.0));

    HybridReachabilityAnalyser analyser(system,domain,accuracy);
    analyser.verbosity = verbosity;

    HybridDenotableSet outer_reach = analyser.outer_chain_reach(initial_set);

    if (plot_results) {
        SystemPlotter plotter(system);
        plotter.plot_set(outer_reach,"outer",accuracy);
    }
}


void infinite_time_lower_evolution(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results) {

    int accuracy = 5;

    HybridBoxes domain(system.state_space(),Box(2,4.5,9.0,0.0,1.0));

    HybridReachabilityAnalyser analyser(system,domain,accuracy);
    analyser.verbosity = verbosity;

	HybridDenotableSet lower_reach;
	HybridFloatVector epsilon;
    make_lpair<HybridDenotableSet,HybridFloatVector>(lower_reach,epsilon) = analyser.lower_chain_reach_and_epsilon(initial_set);

    if (plot_results) {
        SystemPlotter plotter(system);
        plotter.plot_set(lower_reach,"lower",accuracy);
    }
}

void safety_verification(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results) {

    HybridBoxes domain(system.state_space(),Box(2,4.5,9.0,0.0,1.0));
    HybridConstraintSet safety_constraint = getSafetyConstraint(system);

	Verifier verifier;
	verifier.verbosity = verbosity;
	verifier.ttl = 60;
	verifier.settings().plot_results = plot_results;

	SafetyVerificationInput verInput(system, initial_set, domain, safety_constraint);

	verifier.safety(verInput);
}

void parametric_safety_verification(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results) {

    HybridBoxes domain(system.state_space(),Box(2,4.5,9.0,0.0,1.0));
    HybridConstraintSet safety_constraint = getSafetyConstraint(system);

	// The parameters
	RealParameterSet parameters;
	parameters.insert(RealParameter("hmin",Interval(5.0,6.0)));
	parameters.insert(RealParameter("hmax",Interval(7.5,8.5)));

    // Initialization of the verifier
	Verifier verifier;
	verifier.verbosity = verbosity;
	verifier.ttl = 50;
	verifier.settings().plot_results = plot_results;
    verifier.settings().maximum_parameter_depth = 2;

	SafetyVerificationInput verInput(system, initial_set, domain, safety_constraint);

	list<ParametricOutcome> results = verifier.parametric_safety(verInput, parameters);

	if (plot_results)
        draw(system.name(),results);
}

HybridConstraintSet getSafetyConstraint(HybridAutomatonInterface& system) {

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

    return HybridConstraintSet(system.state_space(),ConstraintSet(cons_f,codomain));
}

#endif /* ANALYSIS_H_ */
