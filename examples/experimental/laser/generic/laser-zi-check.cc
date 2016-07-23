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

struct analysis_result {
	double x0;
	double z;
	double zi;
};

analysis_result compute_z(HybridIOAutomaton system, double x0, int verbosity, bool plot_results) {
	double pass_period = 0.03;
    system.substitute(RealParameter("x0",x0));
    Real T0 = system.parameter_value("T0");
    Real velocity = 2.0*system.parameter_value("width")/pass_period;
    system.substitute(RealParameter("velocity",velocity));
    Real vx_i = -velocity;
    Real x_i = 4.0*system.parameter_value("L") + x0;

    /// Create a HybridEvolver object
    HybridEvolver evolver(system);
    evolver.verbosity = verbosity;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(7,0.2);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.000005;
    }

    Box initial_box(7, /*T*/ T0.lower(),T0.upper(), /*p*/ 0.0,0.0, /*t*/ 0.0,0.0, /*vx*/ vx_i.lower(),vx_i.upper(), /*x*/ x_i.lower(),x_i.upper(), /*z*/ 0.0,0.0, /*zi*/ 0.0,0.0);
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,scanning,far,varying,idle"),initial_box);

    int num_half_cycles = 1;
    double evol_time = -8.0*system.parameter_value("L")/vx_i.upper();
    //HybridTime evolution_time(pass_period.upper()/4*num_half_cycles,5*num_half_cycles);
    HybridTime evolution_time(evol_time,7);

    //cout << system << endl;

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    analysis_result result;
    result.x0 = x0;
    result.z = orbit.reach().bounding_box()[5].upper();
    result.zi = orbit.reach().bounding_box()[6].upper();

    std::cout << "Depth of cut : " << result.z << std::endl;
    std::cout << "Maximum value of zi : " << result.zi << std::endl;

    if (plot_results) {
    	PlotHelper plotter(system.name());
    	plotter.plot(orbit.reach(),"reach");
    }

    return result;
}


int main(int argc, char* argv[])
{
    /// Constants
    int accuracy = 11;
	if (argc > 1)
		accuracy = atoi(argv[1]);
    int VERBOSITY = 0;
	if (argc > 2)
		VERBOSITY = atoi(argv[2]);
    bool plot_results = false;
	if (argc > 3)
		plot_results = true;

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


    double left_bound = 0.0;
    double right_bound = 0.0005;//laser_trajectory.parameter_value("width").upper();

    List<analysis_result> results;

    double min_width = pow(10,-accuracy);
    std::cout << "Analyzing down to accuracy " << min_width << std::endl;

    std::cout << "Analyzing initial left at x0 = " << left_bound << std::endl;
    analysis_result left_result = compute_z(system,left_bound,VERBOSITY,plot_results);
    std::cout << "Analyzing initial right at x0 = " << right_bound << std::endl;
    analysis_result right_result = compute_z(system,right_bound,VERBOSITY,plot_results);
    std::cout << "Analyzing initial centre at x0 = " << (left_bound + right_bound)/2.0 << std::endl;
    analysis_result centre_result = compute_z(system,(left_bound + right_bound)/2.0,VERBOSITY,plot_results);

    double eps = 1e-8;

    results.push_back(left_result);
    results.push_back(right_result);
    results.push_back(centre_result);

    bool left = true;
    while (right_result.x0 - left_result.x0 > min_width) {
    	std::cout << std::endl << "#" << results.size() << " : " << Interval(left_result.x0,right_result.x0)
    			<< " (width: " << right_result.x0 - left_result.x0 << ")" << std::endl;

    	if (left) {
			double x0_left = (left_result.x0 + centre_result.x0)/2.0;
			std::cout << "Analyzing left candidate at x0 = " << x0_left << std::endl;
			analysis_result left_candidate = compute_z(system,x0_left,VERBOSITY,plot_results);
			if (left_candidate.zi > centre_result.zi + eps) {
				std::cout << "Left candidate has greater zi value than centre, moving to the left" << std::endl;
				right_result = centre_result;
				centre_result = left_candidate;
				results.push_back(left_candidate);
			} else if (left_candidate.zi >= left_result.zi-eps) {
				std::cout << "Left candidate has no lesser zi value than left, shrinking the left" << std::endl;
				left_result = left_candidate;
				results.push_back(left_candidate);
				double x0_centre = (left_result.x0 + right_result.x0)/2.0;
				std::cout << std::endl << "Analyzing new centre at x0 = " << x0_centre << std::endl;
				centre_result = compute_z(system,x0_centre,VERBOSITY,plot_results);
				results.push_back(centre_result);
			} else
				std::cout << "Left candidate has no greater value than centre, no change" << std::endl;
    	} else {
			double x0_right = (right_result.x0 + centre_result.x0)/2.0;
			std::cout << "Analyzing right candidate at x0 = " << x0_right << std::endl;
			analysis_result right_candidate = compute_z(system,x0_right,VERBOSITY,plot_results);
			if (right_candidate.zi > centre_result.zi + eps) {
				std::cout << "Right candidate has greater zi value than centre, moving to the right" << std::endl;
				left_result = centre_result;
				centre_result = right_candidate;
				results.push_back(right_candidate);
			} else if (right_candidate.zi >= left_result.zi -eps) {
				std::cout << "Right candidate has no lesser zi value than right, shrinking the right" << std::endl;
				right_result = right_candidate;
				results.push_back(right_candidate);
				double x0_centre = (left_result.x0 + right_result.x0)/2.0;
				std::cout << std::endl << "Analyzing new centre candidate at x0 = " << x0_centre << std::endl;
				centre_result = compute_z(system,x0_centre,VERBOSITY,plot_results);
				results.push_back(centre_result);
			} else
				std::cout << "Right candidate has no greater value than centre, no change" << std::endl;
    	}

    	left = !left;
    }

}

