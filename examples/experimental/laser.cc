/*****************************************************************************************************
 *            laser.cc
 *
 *  Testing the laser cutting system.
 *
 *  Copyright  2016  Luca Geretti
 *
 *****************************************************************************************************/

#include "ariadne.h"
#include "laser-trajectory.h"
#include "timer.h"
#include "skin-temperature.h"
#include "cutting-depth.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Constants
    int VERBOSITY = 1;
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Get the automata
    HybridIOAutomaton laser_trajectory = getLaserTrajectory();
    HybridIOAutomaton timer = getTimer();
    HybridIOAutomaton skin_temperature = getSkinTemperature();
    HybridIOAutomaton cutting_depth = getCuttingDepth();

    HybridIOAutomaton timer_traj = compose("timer-traj",timer,laser_trajectory,DiscreteLocation("work"),DiscreteLocation("passing_right"));
    HybridIOAutomaton timer_traj_temperature = compose("timer_traj_temperature",timer_traj,skin_temperature,DiscreteLocation("work,passing_right"),DiscreteLocation("resting"));
    HybridIOAutomaton system = compose("laser",timer_traj_temperature,cutting_depth,DiscreteLocation("work,passing_right,resting"),DiscreteLocation("stable"));

    Real x_i = -laser_trajectory.parameter_value("half_width");
    Real T0 = skin_temperature.parameter_value("T0");
    Real pass_period = 0.05;
    Real vx_i = -4.0*x_i/pass_period;

    /// Create a HybridEvolver object
    HybridEvolver evolver(system);
    evolver.verbosity = VERBOSITY;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(5,2.0);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.0001;
    }

    Box initial_box(5, /*T*/ T0.lower(),T0.upper(), /*t*/ 0.0,0.0, /*vx*/ vx_i.lower(),vx_i.upper(), /*x*/ x_i.lower(),x_i.upper(), /*z*/ 0.0,0.0);
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,passing_right,resting,stable"),initial_box);

    HybridTime evolution_time(pass_period.upper()*8.0,10);

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");
}
