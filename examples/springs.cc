/***************************************************************************
 *            springs.cc
 ****************************************************************************/

#include <cstdarg>
#include "ariadne.h"
#include "springs.h"

using namespace Ariadne;

/// Control variables:
/// x1: position of the first mass
/// x2: position of the second mass
/// v1: speed of the first mass
/// v2: speed of the second mass

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
    HybridIOAutomaton automaton = getSpringsAutomaton();

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(automaton);
    evolver.verbosity = verb;

    /// Set the evolution parameters
    HybridSpace hspace(automaton.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(4,MAX_ENCLOSURE_WIDTH);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = MAX_STEP_SIZE;
    }

    Box initial_box(4, 0.0,0.0, 0.0,0.0, x1_0,x1_0, x2_0,x2_0);
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("free"),initial_box);

    HybridTime evolution_time(EVOL_TIME,EVOL_TRANS);
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);

    PlotHelper plotter(automaton.name());
    plotter.plot(orbit.reach(),"reach");
}
