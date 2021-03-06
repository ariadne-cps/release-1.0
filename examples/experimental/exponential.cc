/*****************************************************************************************************
 *            exponential.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of an exponential.
 *
 *****************************************************************************************************/

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Dynamics parameters
    Vector<Float> dp(11);

    double A = 1.0;

    /// Constants
    float EVOL_TIME = 160.0/A;   /// Evolution time
    float MAX_ENCL_WIDTH = 1e-2;   /// Maximum enclosure width
    float MAX_STEP_SIZE = 1e-2 / A;     /// Maximum step size
    int VERBOSITY = 1;              /// Verbosity of the HybridEvolver
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("exponential");

    /// Create the discrete states
    DiscreteLocation work("work");

    RealVariable t("t");
    RealVariable x("x");

    automaton.add_output_var(t);
    automaton.add_output_var(x);

	automaton.new_mode(work);

	RealExpression dyn_t = 1.0;
	automaton.set_dynamics(work, t, dyn_t);
	RealExpression dyn_x = - A*(x + x*x/2);
	automaton.set_dynamics(work, x, dyn_x);

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(automaton);
    evolver.verbosity = VERBOSITY;

    evolver.settings().set_maximum_step_size(MAX_STEP_SIZE);
    Vector<Float> enclosure_widths(2,0.1,MAX_ENCL_WIDTH);
    evolver.settings().set_reference_enclosure_widths(enclosure_widths);
    evolver.settings().set_maximum_enclosure_widths_ratio(2.0);
    evolver.settings().set_enable_reconditioning(false);
    evolver.settings().set_enable_error_rate_enforcement(true);

    Box initial_box(2, 0.0,0.0, 1.0,1.0);
    HybridEvolver::EnclosureType initial_enclosure(work,initial_box);

    HybridTime evolution_time(EVOL_TIME,1);

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(automaton);
    plotter.plot(orbit.reach(),"reach");
}
