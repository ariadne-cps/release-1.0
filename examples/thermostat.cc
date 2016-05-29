/*****************************************************************************************************
 *            thermostat.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 *****************************************************************************************************/

#include <cstdarg>
#include "ariadne.h"
#include "thermostat.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    int verb = 1;              /// Verbosity of the HybridEvolver
	if (argc > 1)
		verb = atoi(argv[1]);

    /// Create a HybridAutomaton object
    HybridIOAutomaton system = getThermostatSystem();

    HybridEvolver evolver(system);
    evolver.verbosity = verb;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(2,1.0);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 1.0;
    }

    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("flow,start"),Box(2, 8.0,8.0, 0.0,0.0));

    HybridTime evol_limits(1450.0,5);

    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evol_limits,UPPER_SEMANTICS);

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");
}
