/*****************************************************************************************************
 *            circle.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of a circle with proper resets to keep the trajectory stable.
 *
 *****************************************************************************************************/

#include "ariadne.h"
#include "circle.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Constants
    int VERBOSITY = 1;
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton circle = getCircle();

    Real R = circle.parameter_value("R");
    Real w = circle.parameter_value("w");

    /// Create a HybridEvolver object
    HybridEvolver evolver(circle);
    evolver.verbosity = VERBOSITY;

    evolver.settings().set_reference_enclosure_widths(2.0);
    evolver.settings().set_maximum_step_size(0.01);

    Box initial_box(2, R.lower(), R.upper(), 0.0,0.0);
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("1"),initial_box);

    HybridTime evolution_time(1*2.0*Ariadne::pi<Real>().upper()*w.upper(),4);

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(circle.name());
    plotter.plot(orbit.reach(),"reach");
}
