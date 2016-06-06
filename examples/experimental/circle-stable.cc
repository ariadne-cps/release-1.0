/*****************************************************************************************************
 *            circle-stable.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of a circle with proper resets to keep the trajectory stable.
 *
 *****************************************************************************************************/

#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Parameters
	RealParameter R("R",4.0);
	RealParameter w("w",1.0);

    /// Constants
    float EVOL_TIME = 8;   /// Evolution time
    float MAX_ENCL_WIDTH = 0.1;   /// Maximum enclosure width
    float MAX_STEP_SIZE = 1e-2;     /// Maximum step size
    int VERBOSITY = 1;              /// Verbosity of the HybridEvolver
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("circle-stable");

    /// Create the discrete states
    DiscreteLocation first("1");
    DiscreteLocation second("2");
    DiscreteLocation third("3");
    DiscreteLocation fourth("4");

    RealVariable x("x");
    RealVariable y("y");

    automaton.add_output_var(x);
    automaton.add_output_var(y);

    // Events
    DiscreteEvent first2second("12");
    DiscreteEvent second2third("23");
    DiscreteEvent third2fourth("34");
    DiscreteEvent fourth2first("41");

	automaton.new_mode(first);
	automaton.new_mode(second);
	automaton.new_mode(third);
	automaton.new_mode(fourth);

	RealExpression dyn_x = - 2.0*Ariadne::pi<Real>()*w * y;
	automaton.set_dynamics(first, x, dyn_x);
	automaton.set_dynamics(second, x, dyn_x);
	automaton.set_dynamics(third, x, dyn_x);
	automaton.set_dynamics(fourth, x, dyn_x);
	RealExpression dyn_y = 2.0*Ariadne::pi<Real>()*w * x;
	automaton.set_dynamics(first, y, dyn_y);
	automaton.set_dynamics(second, y, dyn_y);
	automaton.set_dynamics(third, y, dyn_y);
	automaton.set_dynamics(fourth, y, dyn_y);

	/// Transitions
	// Guards
	RealExpression x_greater_zero = x;
	RealExpression y_greater_zero = y;
	RealExpression x_lesser_zero = -x;
	RealExpression y_lesser_zero = -y;
	// Resets
	std::map<RealVariable,RealExpression> reset12;
	reset12[x] = 0.0;
	reset12[y] = 1.0;
	std::map<RealVariable,RealExpression> reset23;
	reset23[x] = -1.0;
	reset23[y] = 0.0;
	std::map<RealVariable,RealExpression> reset34;
	reset34[x] = 0.0;
	reset34[y] = -1.0;
	std::map<RealVariable,RealExpression> reset41;
	reset41[x] = 1.0;
	reset41[y] = 0.0;
	automaton.new_forced_transition(first2second,first,second,reset12,x_lesser_zero);
	automaton.new_forced_transition(second2third,second,third,reset23,y_lesser_zero);
	automaton.new_forced_transition(third2fourth,third,fourth,reset34,x_greater_zero);
	automaton.new_forced_transition(fourth2first,fourth,first,reset41,y_greater_zero);

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(automaton);
    evolver.verbosity = VERBOSITY;

    HybridSpace hspace(automaton.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(2,1.0);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.01;
    }

    Box initial_box(2, R.value().lower(), R.value().upper(), 0.0,0.0);
    HybridEvolver::EnclosureType initial_enclosure(first,initial_box);

    HybridTime evolution_time(8*2.0*Ariadne::pi<Real>().upper()*w.value().upper(),32);

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(automaton.name());
    plotter.plot(orbit.reach(),"reach");
}
