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
#include "laser-trajectory.h"
#include "../common/cutting-depth.h"
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

    HybridIOAutomaton timer_traj = compose("timer-traj",timer,laser_trajectory,DiscreteLocation("work"),DiscreteLocation("scanning"));
    HybridIOAutomaton timer_traj_exp = compose("timer_traj_exposure",timer_traj,exposure,DiscreteLocation("work,scanning"),DiscreteLocation("far"));
    //HybridIOAutomaton timer_traj_exp_temp = compose("timer_traj_exp_temp",timer_traj_exp,skin_temperature,DiscreteLocation("work,scanning,far"),DiscreteLocation("varying"));
    //HybridIOAutomaton system = compose("laser",timer_traj_exp_temp,cutting_depth,DiscreteLocation("work,scanning,far,varying"),DiscreteLocation("idle"));
    HybridIOAutomaton system = compose("laser",timer_traj_exp,skin_temperature,DiscreteLocation("work,scanning,far"),DiscreteLocation("varying"));

    Real pass_period = 0.191;

    //Real x0(0.0023);
    //Real x0(0.000150603);
    Real x0(0.000150242);
    // 0.00015160: 1.64067e-05
    // 0.00015165: 1.6363e-05
    // 0.00015166: 1.63543e-05
    // 0.00015167: 1.63458e-05
    // 0.000151675: 1.63414e-05
    // 0.00015168: 1.63371e-05

    system.substitute(RealParameter("x0",x0));
    Real T0 = skin_temperature.parameter_value("T0");
    Real velocity = 2.0*laser_trajectory.parameter_value("width")/pass_period;
    system.substitute(RealParameter("velocity",velocity));
    Real vx_i = -velocity;
    Real x_i = 4.0*exposure.parameter_value("L") + x0;

    /// Create a HybridEvolver object
    HybridEvolver evolver(system);
    evolver.verbosity = VERBOSITY;

    //Vector<Float> reference_enclosure_widths(7,1e-1, 1e-3, 1e-2, 1e-3, 1e-6, 1e-8, 1e-8);
    Vector<Float> reference_enclosure_widths(5,1e-1, 1e-3, 1e-2, 1e-3, 1e-6);
    evolver.settings().set_reference_enclosure_widths(reference_enclosure_widths);
    evolver.settings().set_maximum_enclosure_widths_ratio(1e3);
    evolver.settings().set_maximum_step_size(0.00002);
    evolver.settings().set_enable_reconditioning(true);
    evolver.settings().set_enable_error_rate_enforcement(false);

    //Box initial_box(7, /*T*/ T0.lower(),T0.upper(), /*p*/ 0.0,0.0, /*t*/ 0.0,0.0, /*vx*/ vx_i.lower(),vx_i.upper(), /*x*/ x_i.lower(),x_i.upper(), /*z*/ 0.0,0.0, /*zi*/ 0.0,0.0);
    //HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,scanning,far,varying,idle"),initial_box);
    Box initial_box(5, /*T*/ T0.lower(),T0.upper(), /*p*/ 0.0,0.0, /*t*/ 0.0,0.0, /*vx*/ vx_i.lower(),vx_i.upper(), /*x*/ x_i.lower(),x_i.upper());
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,scanning,far,varying"),initial_box);

    int num_half_cycles = 1;
    double evol_time = -8.0*exposure.parameter_value("L")/vx_i.upper();
    //HybridTime evolution_time(pass_period.upper()/4*num_half_cycles,5*num_half_cycles);
    HybridTime evolution_time(evol_time,7);

    //cout << system << endl;

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(system);
    plotter.plot(orbit.reach(),"reach");

    //std::cout << "Final z : " << 2.0*orbit.final().bounding_box()[5] << std::endl;

    //double depth = orbit.reach().bounding_box()[5].upper();
    //std::cout << "Depth of cut : " << depth << std::endl;
    //std::cout << "Maximum value of zi : " << orbit.reach().bounding_box()[6].upper() << std::endl;
    //std::cout << "Carbonization occurred? " << (orbit.reach().bounding_box()[6].upper() > cutting_depth.parameter_value("z_thr").upper() ? "yes" : "no") << std::endl;
}
