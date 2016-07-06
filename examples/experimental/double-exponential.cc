/*****************************************************************************************************
 *            double-exponential.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of two independent exponentials.
 *
 *****************************************************************************************************/

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    RealParameter A("A",Interval(2.0,2.5));

    /// Constants
    float EVOL_TIME = 40.0/A.value();   /// Evolution time
    float MAX_ENCL_WIDTH = 0.1;   /// Maximum enclosure width
    float MAX_STEP_SIZE = 1e-2 / A.value();     /// Maximum step size
    int VERBOSITY = 1;              /// Verbosity of the HybridEvolver
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("double-exponential");

    /// Create the discrete states
    DiscreteLocation work("work");

    RealVariable x("x");
    RealVariable y("y");

    automaton.add_output_var(y);
    automaton.add_output_var(x);

	automaton.new_mode(work);

	RealExpression dyn_x = - A*(x + x*x/2);
	automaton.set_dynamics(work, x, dyn_x);
	RealExpression dyn_y = - A*(y + y*y/2);
	automaton.set_dynamics(work, y, dyn_y);

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(automaton);
    evolver.verbosity = VERBOSITY;

    evolver.settings().hybrid_maximum_step_size[work] = MAX_STEP_SIZE;
    evolver.settings().minimum_discretised_enclosure_widths[work] = Vector<Float>(2,MAX_ENCL_WIDTH);

    Box initial_box(2, 1.0,1.0, 1.0,1.0);
    HybridEvolver::EnclosureType initial_enclosure(work,initial_box);

    HybridTime evolution_time(EVOL_TIME,1);

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(automaton.name());
    plotter.plot(orbit.reach(),"reach");
}
