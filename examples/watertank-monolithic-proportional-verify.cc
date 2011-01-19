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

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

typedef ImageSetHybridEvolver::EnclosureListType EnclosureListType;


int main(int argc,char *argv[]) 
{
	int analyzerVerbosity = 1;
	if (argc > 1)
		analyzerVerbosity = atoi(argv[1]);

    /// Set the system parameters
	RealConstant a("a",0.02); // The constant defining the decrease rate of the tank level
	RealConstant tau("tau",1.25); // The characteristic time for the opening/closing of the valve
	RealConstant ref("ref",6.75); // A reference tank level
	RealConstant bfp("bfp",Interval(0.125,0.125)); // The product beta*f(p) Interval(0.3,0.32863)
	RealConstant Kp("Kp",5); // The gain of the proportional controller
	RealConstant delta("delta",Interval(-0.0,0.0)); // An indeterminacy in guards evaluation

	// The parameter to modify, its interval and the tolerance
	RealConstant parameter = bfp;
	Interval parameter_interval(0.10,0.125);
	Float tolerance = 1e-2;

    // System variables
	RealVariable x("x"); // water level
	RealVariable y("y"); // valve level
    List<RealVariable> varlist;
    varlist.append(x);
    varlist.append(y);

    // Constants
    ScalarFunction one=ScalarFunction::constant(2,1.0);
    ScalarFunction zero=ScalarFunction::constant(2,0.0);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridAutomaton system;

    /// Create four discrete states
    DiscreteState l1(1);      // Zero saturated
    DiscreteState l2(2);      // Stabilized
    DiscreteState l3(3);      // One saturated

    /// Create the discrete events
    DiscreteEvent e12(12);
    DiscreteEvent e21(21);
    DiscreteEvent e23(23);
    DiscreteEvent e32(32);

    /// Create the dynamics

    RealExpression x_d = -a*x+bfp*y;
    RealExpression y_towardszero_d = -y/tau;
    RealExpression y_controlled_d = (Kp*(ref-x-delta)-y)/tau;
    RealExpression y_towardsone_d = (1-y)/tau;

    // Dynamics at the different modes
    List<RealExpression> exprlist;
    exprlist.append(x_d);
    exprlist.append(y_towardszero_d);
    VectorFunction zerosaturated_d(exprlist, varlist);
    exprlist[1] = y_controlled_d;
    VectorFunction controlled_d(exprlist, varlist);
    exprlist[1] = y_towardsone_d;
    VectorFunction onesaturated_d(exprlist, varlist);

    /// Create the reset
    IdentityFunction reset_id(2);

    /// Create the guards.
    /// Guards are true when f(x) = Ax + b > 0
    /// x <= ref - Delta
    RealExpression x_lesser_ref_minus_delta = -x-delta+ref;
    ScalarFunction guard12(x_lesser_ref_minus_delta,varlist);
    //ScalarAffineFunction guard12(Vector<Float>(4, -1.0,0.0,0.0,-1.0),Rif);
    /// x >= ref - Delta
    RealExpression x_greater_ref_minus_delta = x+delta-ref;
    ScalarFunction guard21(x_greater_ref_minus_delta,varlist);
    //ScalarAffineFunction guard21(Vector<Float>(4, 1.0,0.0,0.0,1.0),-Rif);
    /// x <= ref - 1/Kp - Delta
    RealExpression x_lesser_ref_kp_minus_delta = -x+ref-1.0/Kp-delta;
    ScalarFunction guard23(x_lesser_ref_kp_minus_delta,varlist);
    //ScalarAffineFunction guard23(Vector<Float>(4, -1.0,0.0,0.0,-1.0),(Rif-1.0/Kp));
    /// x >= ref - 1/Kp - Delta
    RealExpression x_greater_ref_kp_minus_delta = x-ref+1.0/Kp+delta;
    ScalarFunction guard32(x_greater_ref_kp_minus_delta,varlist);
    //ScalarAffineFunction guard32(Vector<Float>(4, 1.0,0.0,0.0,1.0),(1.0/Kp - Rif));

    /// Create the invariants.
    /// Invariants are true when f(x) = Ax + b < 0
    /// forced transitions do not need an explicit invariant,
    /// hence we do not need invariants

    /// Build the automaton
    system.new_mode(l1,zerosaturated_d);
    system.new_mode(l2,controlled_d);
    system.new_mode(l3,onesaturated_d);

    system.new_forced_transition(e12,l1,l2,reset_id,guard12);
    system.new_forced_transition(e21,l2,l1,reset_id,guard21);
    system.new_forced_transition(e23,l2,l3,reset_id,guard23);
    system.new_forced_transition(e32,l3,l2,reset_id,guard32);

	// Verification information

	// The initial values
	HybridImageSet initial_set;
	initial_set[l3] = Box(2, 5.5,5.5, 1.0,1.0);

	// The safe region
	HybridBoxes safe_box = bounding_boxes(system.state_space(),Box(2, 5.2, 8.3, -std::numeric_limits<double>::max(), std::numeric_limits<double>::max()));

	// The domain
	HybridBoxes domain = bounding_boxes(system.state_space(),Box(2,4.5,9.0,-0.1,1.1));

	/*
	HybridBoxes domain;
	domain[DiscreteState("flow,opening,rising")] = Box(2,4.5,6.5,-0.1,1.1);
	domain[DiscreteState("flow,closing,falling")] = Box(2,7.0,9.0,-0.1,1.1);
	domain[DiscreteState("flow,idle,falling")] = Box(2,5.0,9.0,-0.1,0.1);
	domain[DiscreteState("flow,idle,rising")] = Box(2,5.0,9.0,0.9,1.1);
	*/

	/// Verification

	// Create an evolver and analyser objects, then set their verbosity
	HybridEvolver evolver;
	evolver.verbosity = 0;
	HybridReachabilityAnalyser analyser(evolver);
	analyser.verbosity = analyzerVerbosity;
	evolver.parameters().enable_subdivisions = false;
	evolver.parameters().enable_set_model_reduction = true;
	analyser.parameters().enable_lower_pruning = true;
	analyser.parameters().lowest_maximum_grid_depth = 0;
	analyser.parameters().highest_maximum_grid_depth = 3;
	analyser.parameters().transient_time = 1e10;
	analyser.parameters().transient_steps = 1;
	analyser.parameters().lock_to_grid_time = 10;
	analyser.parameters().lock_to_grid_steps = 1;
	analyser.plot_verify_results = false;
	analyser.free_cores = 0;
	analyser.chain_reach_dumping = false;

	// The resulting safe and unsafe intervals
	Interval safe_int, unsafe_int;
	// Perform the analysis
	make_lpair(safe_int,unsafe_int) = analyser.safety_unsafety_parametric(system, initial_set, safe_box, domain, parameter, parameter_interval, tolerance);

	cout << "\nResults: " << safe_int << "," << unsafe_int << "\n";

	// Show the result

	if (safe_int == parameter_interval)
		cout << "\nAll the values are safe.\n\n";
	else if (safe_int.empty())
		cout << "\nNo safe value was found.\n\n";
	else if (safe_int.lower() == parameter_interval.lower())
	{
		cout << "\nThe parameter must be <= " << safe_int.upper() << " ( inaccuracy ";
		if (!unsafe_int.empty())
			cout << "<= " << unsafe_int.lower()-safe_int.upper() << ").\n\n";
		else
			cout << "not available).\n\n";
	}
	else if (safe_int.upper() == parameter_interval.upper())
	{
		cout << "\nThe parameter must be >= " << safe_int.lower() << " ( inaccuracy ";
		if (!unsafe_int.empty())
			cout << "<= " << safe_int.lower()-unsafe_int.upper() << ").\n\n";			
		else
			cout << "not available).\n\n";
	}	
	else
		cout << "\nError: the interval could not be verified.\n\n";
}
