/*****************************************************************************************************
 *            sin.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of a sinusoid.
 *
 *****************************************************************************************************/

#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    double A = 2.0;

    /// Constants
    float EVOL_TIME = 1.5;   /// Evolution time
    float MAX_ENCL_WIDTH = 0.1;   /// Maximum enclosure width
    float MAX_STEP_SIZE = 1e-4;     /// Maximum step size
    int VERBOSITY = 1;              /// Verbosity of the HybridEvolver
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("sin");

    /// Create the discrete states
    DiscreteLocation work("work");

    RealVariable x("x");
    RealVariable y("y");

    automaton.add_output_var(x);
    automaton.add_output_var(y);

	automaton.new_mode(work);

	RealExpression dyn_x = 1.0;
	automaton.set_dynamics(work, x, dyn_x);
	RealExpression dyn_y = Ariadne::sin(x);
	automaton.set_dynamics(work, y, dyn_y);


    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(automaton);
    evolver.verbosity = VERBOSITY;

    evolver.settings().set_fixed_maximum_step_size(MAX_STEP_SIZE);
    evolver.settings().set_reference_enclosure_widths(MAX_ENCL_WIDTH);

    Box initial_box(2, 0.0,0.0, 0.0,0.0);
    HybridEvolver::EnclosureType initial_enclosure(work,initial_box);

    HybridTime evolution_time(EVOL_TIME,1);

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(automaton.name());
    plotter.plot(orbit.reach(),"reach");
}
