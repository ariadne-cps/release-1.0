/*****************************************************************************************************
 *            laser.cc
 *
 *  Testing the laser cutting system.
 *
 *  Copyright  2016  Luca Geretti
 *
 *****************************************************************************************************/

#include "ariadne.h"
#include "circle.h"
#include "timer.h"
#include "skin-temperature.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Constants
    int VERBOSITY = 1;
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Get the automata
    HybridIOAutomaton circle = getCircle();
    HybridIOAutomaton timer = getTimer();
    HybridIOAutomaton skin_temperature = getSkinTemperature();

    HybridIOAutomaton timer_circle = compose("timer-circle",timer,circle,DiscreteLocation("work"),DiscreteLocation("1"));
    HybridIOAutomaton system = compose("laser",timer_circle,skin_temperature,DiscreteLocation("work,1"),DiscreteLocation("resting"));

    Real R = circle.parameter_value("R");
    Real w = circle.parameter_value("w");

    /// Create a HybridEvolver object
    HybridEvolver evolver(system);
    evolver.verbosity = VERBOSITY;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(4,2.0);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.0001;
    }

    Box initial_box(4, 0.0,0.0, 0.0,0.0, R.lower(),R.upper(), 0.0,0.0);
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,1,resting"),initial_box);

    HybridTime evolution_time(1*2.0*Ariadne::pi<Real>().upper()*w.upper(),8);

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");
}
