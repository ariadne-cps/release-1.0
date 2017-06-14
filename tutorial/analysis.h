/***************************************************************************
 *            analysis.h
 *
 *  Provides a sequence of execution of analysis functions.
 *  For each one, the input data is prepared (apart from the common initial set).
 *
 *  Copyright  2017  Luca Geretti
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

/// Forward declarations, used only to properly organize the source file
void finite_time_upper_evolution(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results);
void infinite_time_outer_evolution(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results);
void infinite_time_lower_evolution(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results);
void safety_verification(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results);
void parametric_safety_verification(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results);
HybridConstraintSet getSafetyConstraint(HybridAutomatonInterface& system);

// The main method for the analysis of the system
// Since the analyses are independent, you may comment out any one if you want
// to focus on specific ones.
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

// Performs finite time evolution, using upper semantics.
void finite_time_upper_evolution(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results) {

    HybridEvolver evolver(system);
    evolver.verbosity = verbosity;

    HybridEvolver::EnclosureType initial_enclosure;
    HybridBoxes initial_set_domain = initial_set.domain();
    for (std::map<DiscreteLocation,Box>::const_iterator it = initial_set_domain.locations_begin(); it != initial_set_domain.locations_end(); ++it) {
    	if (!it->second.empty()) {
    		initial_enclosure = HybridEvolver::EnclosureType(it->first,Box(it->second.centre()));
    		break;
    	}
    }
  
    HybridTime evol_limits(30.0,8);
 
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evol_limits,UPPER_SEMANTICS);

    if (plot_results) {
        PlotHelper plotter(system.name());
        plotter.plot(orbit.reach(),"reach");
    }
}

// Performs infinite time outer evolution
void infinite_time_outer_evolution(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results) {

    int accuracy = 5;

    HybridBoxes domain(system.state_space(),Box(2,0.0,1.0,4.5,9.0));

    HybridReachabilityAnalyser analyser(system,domain,accuracy);
    analyser.verbosity = verbosity;

    HybridDenotableSet outer_reach = analyser.outer_chain_reach(initial_set);

    if (plot_results) {
        PlotHelper plotter(system.name());
        plotter.plot(outer_reach,"outer",accuracy);
    }
}

// Performs infinite time epsilon-lower evolution
void infinite_time_lower_evolution(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results) {

    int accuracy = 5;

    HybridBoxes domain(system.state_space(),Box(2,0.0,1.0,4.5,9.0));

    HybridReachabilityAnalyser analyser(system,domain,accuracy);
    analyser.verbosity = verbosity;

    HybridDenotableSet lower_reach;
    HybridFloatVector epsilon;
    make_lpair<HybridDenotableSet,HybridFloatVector>(lower_reach,epsilon) = analyser.lower_chain_reach_and_epsilon(initial_set);

    if (plot_results) {
        PlotHelper plotter(system.name());
        plotter.plot(lower_reach,"lower",accuracy);
    }
}

// Performs verification in respect to a safety specification expresses as a set
void safety_verification(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results) {

	// Creates the domain, necessary to guarantee termination for outer evolution
    HybridBoxes domain(system.state_space(),Box(2,0.0,1.0,4.5,9.0));
    // Creates the safety constraint
    HybridConstraintSet safety_constraint = getSafetyConstraint(system);

    // Initializes the verifier
    Verifier verifier;
    verifier.verbosity = verbosity;
    verifier.settings().plot_results = plot_results;
    // The time (in seconds) after which we stop verification
    verifier.ttl = 140;

    // Collects the verification input
    SafetyVerificationInput verInput(system, initial_set, domain, safety_constraint);

    // Performs verification
    verifier.safety(verInput);
}

// Performs verification in respect to a safety specification expresses as a set,
// but it does such verification within a given parameters space, where hmin and hmax
// are expresses as intervals. Such intervals are then split in order to identify
// a collection of boxes where to perform the safety verification individually.
void parametric_safety_verification(HybridAutomatonInterface& system, HybridBoundedConstraintSet& initial_set, int verbosity, bool plot_results) {

	// Creates the domain, necessary to guarantee termination for outer evolution
    HybridBoxes domain(system.state_space(),Box(2,0.0,1.0,4.5,9.0));
    // Creates the safety constraint
    HybridConstraintSet safety_constraint = getSafetyConstraint(system);

    // The parameters which will be split into disjoint sets
    RealParameterSet parameters;
    parameters.insert(RealParameter("hmin",Interval(5.0,6.0)));
    parameters.insert(RealParameter("hmax",Interval(7.5,8.5)));

    // Initialization of the verifier
    Verifier verifier;
    verifier.verbosity = verbosity;
    verifier.settings().plot_results = plot_results;
    // The time (in seconds) after which we stop verification for this split parameters set and move to another set
    verifier.ttl = 140;
    // The number of consecutive splittings for each parameter, i.e., 2^value.
    // In this case we allow 8x8 = 64 disjoint sets. The larger this number, the more accurate the result for each disjoint set.
    verifier.settings().maximum_parameter_depth = 3;

    // Collects the verification input
    SafetyVerificationInput verInput(system, initial_set, domain, safety_constraint);

    // Performs verification, saving the results as a list for each split set
    list<ParametricOutcome> results = verifier.parametric_safety(verInput, parameters);

    // Plots the list in a 2d mesh
    if (plot_results) {
        PlotHelper plotter(system.name());
        plotter.plot(results,verifier.settings().maximum_parameter_depth);
    }
}

// Constructs the safety constraint for (parametric) safety verification
HybridConstraintSet getSafetyConstraint(HybridAutomatonInterface& system) {

	// The desired constraint is 5.52 <= x <= 8.25 for all locations
	// The construction below may seem convoluted, however it allows for
	// large generality when defining the constraint set.

	// Constructs the variable list, required by the vector function
    RealVariable x("x");
    RealVariable a("a");
    List<RealVariable> varlist;
    varlist.append(x);
    varlist.append(a);
    // Constructs the expression
    RealExpression expr = x;
    List<RealExpression> consexpr;
    consexpr.append(expr);
    VectorFunction cons_f(consexpr,varlist);
    // Constructs the codomain for the expression
    Box codomain(1,5.52,8.25);

    // Constructs a costraint set and then applies it to each location of the system
    return HybridConstraintSet(system.state_space(),ConstraintSet(cons_f,codomain));
}

#endif /* ANALYSIS_H_ */
