/*****************************************************************************************************
 *            laser.cc
 *
 *  Testing the laser cutting system.
 *
 *  Copyright  2016  Luca Geretti
 *
 *****************************************************************************************************/

#include "ariadne.h"
#include "../timer.h"
#include "laser-trajectory.h"
#include "skin-temperature.h"
#include "skin-exposure.h"
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
    HybridIOAutomaton exposure = getSkinExposure();
    HybridIOAutomaton skin_temperature = getSkinTemperature();
    HybridIOAutomaton cutting_depth = getCuttingDepth();

    HybridIOAutomaton timer_traj = compose("timer-traj",timer,laser_trajectory,DiscreteLocation("work"),DiscreteLocation("passing_right"));
    HybridIOAutomaton timer_traj_exp = compose("timer_traj_exposure",timer_traj,exposure,DiscreteLocation("work,passing_right"),DiscreteLocation("far_from_left"));
    HybridIOAutomaton timer_traj_exp_temp = compose("timer_traj_exp_temp",timer_traj_exp,skin_temperature,DiscreteLocation("work,passing_right,far_from_left"),DiscreteLocation("varying"));
    HybridIOAutomaton system = compose("laser",timer_traj_exp_temp,cutting_depth,DiscreteLocation("work,passing_right,far_from_left,varying"),DiscreteLocation("idle"));

    Real x_i = -laser_trajectory.parameter_value("half_width");
    Real T0 = skin_temperature.parameter_value("T0");
    Real pass_period = 0.05;
    Real vx_i = -4.0*x_i/pass_period;

    /// Create a HybridEvolver object
    HybridEvolver evolver(system);
    evolver.verbosity = VERBOSITY;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(7,2.0);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.00001;
    }

    Box initial_box(7, /*T*/ T0.lower(),T0.upper(), /*p*/ 0.0,0.0, /*t*/ 0.0,0.0, /*vx*/ vx_i.lower(),vx_i.upper(), /*x*/ x_i.lower(),x_i.upper(), /*z*/ 0.0,0.0, /*zi*/ 0.0,0.0);
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,passing_right,far_from_left,varying,idle"),initial_box);

    int num_cycles = 1.0;
    HybridTime evolution_time(pass_period.upper()*20.0,5);

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");
}
