/*****************************************************************************************************
 *            home.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 *****************************************************************************************************/

#include <cstdarg>
#include "ariadne.h"
#include "home.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    int verb = 1;
	if (argc > 1)
		verb = atoi(argv[1]);

    HybridIOAutomaton system = getHomeSystem();

    HybridEvolver evolver(system);
    evolver.verbosity = verb;

    evolver.settings().set_reference_enclosure_widths(1.0);
    evolver.settings().set_fixed_maximum_step_size(1.0);

    Real day_len = system.parameter_value("day_len");
    Real Te_max = system.parameter_value("Te_max");
    Real Te_min = system.parameter_value("Te_min");
    Real phi = system.parameter_value("phi");
    Real Te_0 = Te_min + (Te_max - Te_min)*(1.0-cos(phi))/2;

    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("flow,oscillate,night,unregulated"),Box(4, Te_0.lower(),Te_0.upper(), Te_0.lower(),Te_0.upper(), 0.0,0.0, 0.0, 0.0));

    HybridTime evol_limits(2.0*day_len.midpoint(),7);

    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evol_limits,UPPER_SEMANTICS);

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");
}
