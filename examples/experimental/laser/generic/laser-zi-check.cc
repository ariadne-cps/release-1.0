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
    bool plot_results = false;
	if (argc > 3)
		plot_results = true;

    Real pass_period = 0.04;

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

    // Prepare the points
    Array<double> points(num_points);
    points[0] = 0.0;
    points[num_points-1] = laser_trajectory.parameter_value("half_width").upper();
    for (int i = 1; i<num_points-1; ++i) {
    	points[i] = (points[num_points-1] + points[i-1])/2.0;
    }

    for (int i = 0; i<num_points; ++i) {
    	double x0 = points[i];
    	std::cout << std::endl << "#" << i << " : Analyzing for x0 = " << x0 << std::endl;
        system.substitute(RealParameter("x0",x0));
        Real T0 = skin_temperature.parameter_value("T0");
        Real vx = 4.0*laser_trajectory.parameter_value("half_width")/pass_period;
        system.substitute(RealParameter("velocity",vx));
        Real x_i = -4.0*exposure.parameter_value("L") + x0;

        /// Create a HybridEvolver object
        HybridEvolver evolver(system);
        evolver.verbosity = VERBOSITY;

        HybridSpace hspace(system.state_space());
        for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
            evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(7,2.0);
            evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.000005;
        }

        Box initial_box(7, /*T*/ T0.lower(),T0.upper(), /*p*/ 0.0,0.0, /*t*/ 0.0,0.0, /*vx*/ vx.lower(),vx.upper(), /*x*/ x_i.lower(),x_i.upper(), /*z*/ 0.0,0.0, /*zi*/ 0.0,0.0);
        HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,scanning,far,varying,idle"),initial_box);

        int num_half_cycles = 1;
        double evol_time = 8.0*exposure.parameter_value("L")/vx.upper();
        //HybridTime evolution_time(pass_period.upper()/4*num_half_cycles,5*num_half_cycles);
        HybridTime evolution_time(evol_time,5);

        //cout << system << endl;

        std::cout << "Computing orbit... " << std::flush;
        HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
        std::cout << "done." << std::endl;

        double depth = orbit.reach().bounding_box()[5].upper();
        std::cout << "Depth of cut : " << depth << std::endl;
        std::cout << "Maximum value of zi : " << orbit.reach().bounding_box()[6].upper() << std::endl;
        std::cout << "Carbonization occurred? " << (depth > cutting_depth.parameter_value("z_thr").upper() ? "yes" : "no") << std::endl;

        if (plot_results) {
        	PlotHelper plotter(system.name());
        	plotter.plot(orbit.reach(),"reach");
        }
    }
}
