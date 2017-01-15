/*****************************************************************************************************
 *            painting.cc
 *
 *  Testing the robotic spray painting system.
 *
 *  Copyright  2016  Luca Geretti
 *
 *****************************************************************************************************/

#include "ariadne.h"

#include "../timer.h"
#include "spray.h"
#include "trajectory.h"
#include "deposition.h"
#include <iostream>
#include <fstream>

using namespace Ariadne;

struct d_velocity_pair {
	double d;
	double velocity;
};

int count_lines() {
    int number_of_lines = 0;
    std::string line;
    std::ifstream myfile("results.txt");

    while (std::getline(myfile, line))
        ++number_of_lines;

    return number_of_lines;
}

int main(int argc, char* argv[])
{
    /// Build the Hybrid System

    /// Get the automata
    HybridIOAutomaton trajectory = getTrajectory();
    HybridIOAutomaton timer = getTimer();
    HybridIOAutomaton spray = getSpray();
    HybridIOAutomaton deposition = getDeposition();

    HybridIOAutomaton timer_traj = compose("timer-traj",timer,trajectory,DiscreteLocation("work"),DiscreteLocation("scanning"));
    HybridIOAutomaton timer_traj_spray = compose("timer_traj_spray",timer_traj,spray,DiscreteLocation("work,scanning"),DiscreteLocation("far"));
    HybridIOAutomaton system = compose("painting",timer_traj_spray,deposition,DiscreteLocation("work,scanning,far"),DiscreteLocation("accumulating"));

    Real x0(0.015);
    system.substitute(RealParameter("x0",x0));
    Real x_i = 0.0;
    Real y_i = 0.0;
    Real delta(0.00005);

    // Interval d_range(0.0025,0.01);
    Interval velocity_range(0.05,0.4);
    std::list<d_velocity_pair> dv_values;

    for (int i = 25; i<= 100; i+=4) {
    	for (double j = velocity_range.lower(); j <= velocity_range.upper(); j+=0.02) {
    		d_velocity_pair dv;
    		dv.d = 0.01 * i / 100.0;
    		dv.velocity = j;
    		dv_values.push_back(dv);
    	}
    }

    int num_lines = count_lines();
    int total_values = dv_values.size();

    for (int j = 0; j < total_values; ++j) {

    	d_velocity_pair dv = dv_values.front();
    	dv_values.pop_front();

    	if (j < num_lines)
    		continue;

    	double d = dv.d;
    	double velocity = dv.velocity;

		system.substitute(RealParameter("d",d));
		Real vx_i = velocity;

		double length = 4.0*d;
		int num_intervals = ceil(d/2.0/delta);
		Vector<Interval> z_values(num_intervals);

		std::cout << "#" << (j+1) << "/" << total_values <<
				  ": d=" << d << ", v=" << velocity << " (processing " << num_intervals << " intervals: ";

		Interval bounds(-1);

		for (int i = 0; i < num_intervals; ++i) {

			Real y0(length,length+delta);

			/// Create a HybridEvolver object
			HybridEvolver evolver(system);
			evolver.verbosity = 0;

			evolver.settings().set_reference_enclosure_widths(1.0);
			evolver.settings().set_fixed_maximum_step_size(0.0005/velocity);
			evolver.settings().set_maximum_number_of_working_sets(20);

			Box initial_box(6,
					/*s*/ 0.0,0.0,
					/*t*/ 0.0,0.0,
					/*vx*/ vx_i.lower(),vx_i.upper(),
					/*x*/ x_i.lower(),x_i.upper(),
					/*y*/ y_i.lower(),y_i.upper(),
					/*z*/ 0.0,0.0);
			HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,scanning,far,accumulating"),initial_box);

			int num_passes = 8.0;
			double evol_time = num_passes * velocity / trajectory.parameter_value("width").lower();
			HybridTime evolution_time(evol_time,num_passes*3);

			system.substitute(RealParameter("y0",y0));

			try {
				HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
				Interval local_bounds = orbit.final().bounding_box()[5];
				double new_lower = bounds.lower();
				if (bounds.lower() == -1 || bounds.lower() > local_bounds.lower())
					new_lower = local_bounds.lower();
				double new_upper = max(bounds.upper(),local_bounds.upper());
				bounds.set(new_lower, new_upper);

				std::cout << "v" << std::flush;
			} catch (WorkingSetTooLargeException& ex) {
				std::cout << "x" << std::flush;
			}

			length += delta;
		}

		ofstream myfile;
		myfile.open("results.txt", ios::out | ios::app);
		myfile << d << " "
			   << velocity << " "
			   << bounds.lower() << " "
			   << bounds.upper() << std::endl;
		myfile.close();

		std::cout << ")" << std::endl;
    }

}
