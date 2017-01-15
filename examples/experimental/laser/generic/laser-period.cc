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
	int num_points = 2;
	if (argc > 1)
		num_points = atoi(argv[1]);
    int VERBOSITY = 0;
	if (argc > 2)
		VERBOSITY = atoi(argv[2]);

    /// Build the Hybrid System

    /// Get the automata
    HybridIOAutomaton laser_trajectory = getLaserTrajectory();
    HybridIOAutomaton timer = getTimer();
    HybridIOAutomaton exposure = getSkinExposure();
    HybridIOAutomaton skin_temperature = getSkinTemperature();
    HybridIOAutomaton cutting_depth = getCuttingDepth();

    HybridIOAutomaton timer_traj = compose("timer-traj",timer,laser_trajectory,DiscreteLocation("work"),DiscreteLocation("scanning"));
    HybridIOAutomaton timer_traj_exp = compose("timer_traj_exposure",timer_traj,exposure,DiscreteLocation("work,scanning"),DiscreteLocation("far"));
    HybridIOAutomaton timer_traj_exp_temp = compose("timer_traj_exp_temp",timer_traj_exp,skin_temperature,DiscreteLocation("work,scanning,far"),DiscreteLocation("varying"));
    HybridIOAutomaton system = compose("laser",timer_traj_exp_temp,cutting_depth,DiscreteLocation("work,scanning,far,varying"),DiscreteLocation("idle"));

    Real x0 = 0.5*laser_trajectory.parameter_value("width");

    system.substitute(RealParameter("x0",x0));
    Real T0 = skin_temperature.parameter_value("T0");

    Array<Float> periods(num_points);

    for (int i=0; i<num_points;++i) {
    	double pass_period = 0.104 + round(1000.0*(0.060)*i/(num_points-1))/1000.0;
    	std::cout << "#" << i << ": pass period " << pass_period << std::endl;
    }
    for (int i=0; i<num_points;++i) {
    	double pass_period = 0.104 + round(1000.0*(0.060)*i/(num_points-1))/1000.0;
    	std::cout << "#" << i << ": pass period " << pass_period << std::endl;

		Real velocity = 2.0*laser_trajectory.parameter_value("width")/pass_period;
		system.substitute(RealParameter("velocity",velocity));
		Real vx_i = -velocity;
		Real x_i = 4.0*exposure.parameter_value("L") + x0;

		/// Create a HybridEvolver object
		HybridEvolver evolver(system);
		evolver.verbosity = VERBOSITY;

		evolver.settings().set_reference_enclosure_widths(0.5);
	    evolver.settings().set_fixed_maximum_step_size(0.000001);

		Box initial_box(7, /*T*/ T0.lower(),T0.upper(), /*p*/ 0.0,0.0, /*t*/ 0.0,0.0, /*vx*/ vx_i.lower(),vx_i.upper(), /*x*/ x_i.lower(),x_i.upper(), /*z*/ 0.0,0.0, /*zi*/ 0.0,0.0);
		HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,scanning,far,varying,idle"),initial_box);

		int num_half_cycles = 1;
		double evol_time = -8.0*exposure.parameter_value("L")/vx_i.upper();
		//HybridTime evolution_time(pass_period.upper()/4*num_half_cycles,5*num_half_cycles);
		HybridTime evolution_time(evol_time,7);

		//cout << system << endl;

		std::cout << "Computing orbit... " << std::flush;
		HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
		std::cout << "done." << std::endl;

		double depth = orbit.final().bounding_box()[5].upper()*2;
		std::cout << "Depth of cut : " << depth << std::endl;
    }
}
