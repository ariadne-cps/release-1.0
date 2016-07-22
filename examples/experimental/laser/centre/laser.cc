/*****************************************************************************************************
 *            laser.cc
 *
 *  Testing the laser cutting system.
 *
 *  Copyright  2016  Luca Geretti
 *
 *****************************************************************************************************/

#include "ariadne.h"

#include "skin-exposure.h"
#include "../common/cutting-depth.h"
#include "../common/laser-trajectory.h"
#include "../common/skin-temperature.h"
#include "../../timer.h"

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
    HybridIOAutomaton timer_traj_exp = compose("timer_traj_exposure",timer_traj,exposure,DiscreteLocation("work,passing_right"),DiscreteLocation("far_from_left_in"));
    HybridIOAutomaton timer_traj_exp_temp = compose("timer_traj_exp_temp",timer_traj_exp,skin_temperature,DiscreteLocation("work,passing_right,far_from_left_in"),DiscreteLocation("varying"));
    HybridIOAutomaton system = compose("laser",timer_traj_exp_temp,cutting_depth,DiscreteLocation("work,passing_right,far_from_left_in,varying"),DiscreteLocation("idle"));

    Real x_i = -laser_trajectory.parameter_value("half_width");
    Real T0 = skin_temperature.parameter_value("T0");
    Real pass_period = 0.012;
    Real vx = 4.0*laser_trajectory.parameter_value("half_width")/pass_period;

    system.substitute(RealParameter("velocity",vx));

    /// Create a HybridEvolver object
    HybridEvolver evolver(system);
    evolver.verbosity = VERBOSITY;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(6,2.0);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.000005;
    }

    Box initial_box(6, /*T*/ T0.lower(),T0.upper(), /*p*/ 0.0,0.0, /*t*/ 0.0,0.0, /*x*/ x_i.lower(),x_i.upper(), /*z*/ 0.0,0.0, /*zi*/ 0.0,0.0);
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,passing_right,far_from_left_in,varying,idle"),initial_box);

    int num_half_cycles = 1;
    HybridTime evolution_time(pass_period.upper()/2*num_half_cycles,7*num_half_cycles);

    //cout << system << endl;


    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");

}
