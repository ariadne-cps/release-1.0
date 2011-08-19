/***************************************************************************
 *            traffic-arbitrator4.cc
 *
 *  Traffic arbitrator example for 4 cars.
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

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int thisVerbosity = 1;
	if (argc > 1)
		thisVerbosity = atoi(argv[1]);

	HybridAutomaton system("traffic-arbitrator4");

    /// Set the system parameters
	Float a = 10.0;
	Float b = 30.0;
	Float c = 40.0;
	Float d = 50.0;
	Float e = 60.0;
	Float f = 100.0;
	Float rl = 20.0;
	Float ru = 70.0;
	Float alpha_min_prime = 0.0005;
	Float alpha_initial_l = 0.1;
	Float alpha_initial_u = 0.1;
	RealParameter af_fa("af_fa",Interval(a - f, f - a));
	RealParameter dc_eb("dc_eb",Interval(d - c, e - b));
	RealParameter brl_crl("brl_crl",Interval(b - rl, c - rl));
	RealParameter rue_rud("rue_rud",Interval(ru - e, ru - d));

	RealParameter alpha_min("alpha_min",0.002);
	RealParameter alpha_max("alpha_max",1.0);
	RealParameter delta("delta",0.1);

    /// Create the discrete states
    DiscreteLocation cruiseStart("cruiseStart");
    DiscreteLocation recovery12("recovery12");
    DiscreteLocation recovery23("recovery23");
    DiscreteLocation recovery34("recovery34");

    /// Create the discrete events
    DiscreteEvent warning12("warning12");
    DiscreteEvent warning23("warning23");
    DiscreteEvent warning34("warning34");
    DiscreteEvent recovered12("recovered12");
    DiscreteEvent recovered23("recovered23");
    DiscreteEvent recovered34("recovered34");

    // System variables
    RealVariable y12("y12");
    RealVariable y23("y23");
    RealVariable y34("y34");
    List<RealVariable> varlist;
    varlist.append(y12);
    varlist.append(y23);
    varlist.append(y34);

    // Expressions

      // invariant kind (f(x) <= 0)

    // y12 <= alpha_max
    RealExpression y12_leq_alpha_max = y12 - alpha_max;
    // y23 <= alpha_max
    RealExpression y23_leq_alpha_max = y23 - alpha_max;
    // y34 <= alpha_max
    RealExpression y34_leq_alpha_max = y34 - alpha_max;

      // activation kind (f(x) >= 0)

    // y12 <= alpha_min
    RealExpression y12_leq_alpha_min = alpha_min - y12;
    // y12 >= alpha_min
    RealExpression y12_geq_alpha_min = y12 - alpha_min;
    // y23 <= alpha_min
    RealExpression y23_leq_alpha_min = alpha_min - y23;
    // y23 >= alpha_min
    RealExpression y23_geq_alpha_min = y23 - alpha_min;
    // y34 <= alpha_min
    RealExpression y34_leq_alpha_min = alpha_min - y34;
    // y34 >= alpha_min
    RealExpression y34_geq_alpha_min = y34 - alpha_min;

    // Identity
    RealExpression id_y12 = y12;
    RealExpression id_y23 = y23;
    RealExpression id_y34 = y34;
    // Zero
    RealExpression zero = 0.0;

    // Dynamics at the different modes

    List<RealExpression> exprlist;

    exprlist.append(af_fa);
    exprlist.append(af_fa);
    exprlist.append(af_fa);
    VectorFunction cruiseStart_d(exprlist, varlist);

    exprlist[0] = dc_eb;
    exprlist[1] = brl_crl;
    exprlist[2] = zero;
    VectorFunction recovery12_d(exprlist, varlist);

    exprlist[0] = rue_rud;
    exprlist[1] = dc_eb;
    exprlist[2] = brl_crl;
    VectorFunction recovery23_d(exprlist, varlist);

    exprlist[0] = zero;
    exprlist[1] = rue_rud;
    exprlist[2] = dc_eb;
    VectorFunction recovery34_d(exprlist, varlist);

    // Reset functions

    exprlist[0] = id_y12;
    exprlist[0] = id_y23;
    exprlist[0] = id_y34;
    VectorFunction identity(exprlist, varlist);

    // Guards

    ScalarFunction warning12_g(y12_leq_alpha_min, varlist);
    ScalarFunction warning23_g(y23_leq_alpha_min, varlist);
    ScalarFunction warning34_g(y34_leq_alpha_min, varlist);
    ScalarFunction recovered12_g(y12_geq_alpha_min, varlist);
    ScalarFunction recovered23_g(y23_geq_alpha_min, varlist);
    ScalarFunction recovered34_g(y34_geq_alpha_min, varlist);

    // Invariants

    ScalarFunction y12_i(y12_leq_alpha_max, varlist);
    ScalarFunction y23_i(y23_leq_alpha_max, varlist);
    ScalarFunction y34_i(y34_leq_alpha_max, varlist);

    /// Build the automaton

    system.new_mode(cruiseStart,cruiseStart_d);
    system.new_mode(recovery12,recovery12_d);
    system.new_mode(recovery23,recovery23_d);
    system.new_mode(recovery34,recovery34_d);

    system.new_invariant(cruiseStart,y12_i);
    system.new_invariant(cruiseStart,y23_i);
    system.new_invariant(cruiseStart,y34_i);
    system.new_invariant(recovery12,y12_i);
    system.new_invariant(recovery12,y23_i);
    system.new_invariant(recovery12,y34_i);
    system.new_invariant(recovery23,y12_i);
    system.new_invariant(recovery23,y23_i);
    system.new_invariant(recovery23,y34_i);
    system.new_invariant(recovery34,y12_i);
    system.new_invariant(recovery34,y23_i);
    system.new_invariant(recovery34,y34_i);

    system.new_forced_transition(warning12,cruiseStart,recovery12,identity,warning12_g);
    system.new_forced_transition(warning23,cruiseStart,recovery23,identity,warning23_g);
    system.new_forced_transition(warning34,cruiseStart,recovery34,identity,warning34_g);
    system.new_forced_transition(recovered12,recovery12,cruiseStart,identity,recovered12_g);
    system.new_forced_transition(recovered23,recovery23,cruiseStart,identity,recovered23_g);
    system.new_forced_transition(recovered34,recovery34,cruiseStart,identity,recovered34_g);

    // The initial values
    HybridImageSet initial_set;
    initial_set[cruiseStart] = Box(3, alpha_initial_l,alpha_initial_u, alpha_initial_l,alpha_initial_u, alpha_initial_l,alpha_initial_u);

    // The domain
    HybridBoxes domain = bounding_boxes(system.state_space(),
            Box(3,0.0,1.1*alpha_max.value().upper(),0.0,1.1*alpha_max.value().upper(),0.0,1.1*alpha_max.value().upper()));

    // The safety constraint
    Box codomain(3,
            alpha_min_prime,std::numeric_limits<double>::infinity(),
            alpha_min_prime,std::numeric_limits<double>::infinity(),
            alpha_min_prime,std::numeric_limits<double>::infinity());
    HybridConstraintSet safety_constraint(system.state_space(),ConstraintSet(identity,codomain));

    /// Verification

    HybridReachabilityAnalyser analyser(system);
    analyser.settings().highest_maximum_grid_depth = 6;
    Verifier verifier(analyser);
    verifier.set_verbosity(thisVerbosity);
    verifier.settings().plot_results = true;

    SafetyVerificationInput verInfo(system, initial_set, domain, safety_constraint);
    cout << verifier.safety(verInfo) << "\n";
}
