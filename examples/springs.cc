/***************************************************************************
 *            springs.cc
 ****************************************************************************/

#include <cstdarg>
#include "ariadne.h"
#include "springs.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{    
    int verb = 1;
	if (argc > 1)
		verb = atoi(argv[1]);

    float x1_0 = 0.0; // Initial position for the first spring
    float x2_0 = 3.0; // Initial position for the second spring

    float EVOL_TIME = 40.0; // Evolution time
    int   EVOL_TRANS = 4; // Evolution transitions
    float MAX_ENCLOSURE_WIDTH = 0.02; // Maximum enclosure width
    float MAX_STEP_SIZE = 0.01; // Maximum integration step size
  
    /// Create a HybridAutomaton object
    HybridIOAutomaton system = getSpringsSystem();

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(system);
    evolver.verbosity = verb;

    /// Set the evolution parameters
    evolver.settings().set_reference_enclosure_widths(MAX_ENCLOSURE_WIDTH);
    evolver.settings().set_hybrid_maximum_step_size(MAX_STEP_SIZE);

    Box initial_box(5, 0.0, 0.0, 0.0,0.0, 0.0,0.0, x1_0,x1_0, x2_0,x2_0);
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("free1,free2"),initial_box);

    HybridTime evolution_time(EVOL_TIME,EVOL_TRANS);
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");
}
