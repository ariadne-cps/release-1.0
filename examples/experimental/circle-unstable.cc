/*****************************************************************************************************
 *            circle-unstable.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of a circle.
 *
 *****************************************************************************************************/

#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Dynamics parameters
    Vector<Float> dp(11);

    double A = 2.0;

    /// Constants
    float EVOL_TIME = 16;   /// Evolution time
    float SCALING = 1e-4;   /// Scaling of both variables
    float MAX_ENCL_WIDTH_RATIO = 100.0; // Ratio for the maximum enclosure in respect to the scaling
    float FIXED_MAXIMUM_STEP_SIZE = 4e-3;     /// Fixed maximum step size
    int VERBOSITY = 1;              /// Verbosity of the HybridEvolver
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("circle-unstable");

    /// Create the discrete states
    DiscreteLocation work("work");

    RealVariable x("x");
    RealVariable y("y");

    automaton.add_output_var(x);
    automaton.add_output_var(y);

	automaton.new_mode(work);

	RealExpression dyn_x = - y;
	automaton.set_dynamics(work, x, dyn_x);
	RealExpression dyn_y = x;
	automaton.set_dynamics(work, y, dyn_y);


    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(automaton);
    evolver.verbosity = VERBOSITY;

    evolver.settings().set_fixed_maximum_step_size(FIXED_MAXIMUM_STEP_SIZE);
    evolver.settings().set_reference_enclosure_widths(SCALING);
    evolver.settings().set_maximum_enclosure_widths_ratio(MAX_ENCL_WIDTH_RATIO);
    evolver.settings().set_enable_adaptive_maximum_step_size(false);
    evolver.settings().set_enable_reconditioning(true);
    evolver.settings().set_enable_boxing_on_contraction(false);

    Box initial_box(2, 0.0,0.0, 1.0,1.0);
    HybridEvolver::EnclosureType initial_enclosure(work,initial_box);

    HybridTime evolution_time(EVOL_TIME,1);

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(automaton.name());
    plotter.plot(orbit.reach(),"reach");
}
