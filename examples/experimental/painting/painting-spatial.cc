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

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Constants
    int VERBOSITY = 1;
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

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

    Real velocity(0.3);
    Real d(0.0025);
    system.substitute(RealParameter("d",d));
    Real vx_i = velocity;
    Real x_i = 0.0;
    Real y_i = 0.0;

    Real delta(0.00005);

    double length = 4.0*d;
    int num_intervals = ceil(d/2.0/delta);
    Vector<Interval> z_values(num_intervals);

    for (int i = 0; i < num_intervals; ++i) {

    	Real y0(length,length+delta);
    	std::cout << (i+1) << "/" << num_intervals << ": interval " << y0 << std::flush << std::endl;

        /// Create a HybridEvolver object
        HybridEvolver evolver(system);
        evolver.verbosity = VERBOSITY;

        HybridSpace hspace(system.state_space());
        for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
            evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(6,1.0);
            evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.001;
            evolver.settings().maximum_number_of_working_sets = 20;
        }

        Box initial_box(6,
    			/*s*/ 0.0,0.0,
    			/*t*/ 0.0,0.0,
    			/*vx*/ vx_i.lower(),vx_i.upper(),
    			/*x*/ x_i.lower(),x_i.upper(),
    			/*y*/ y_i.lower(),y_i.upper(),
    			/*z*/ 0.0,0.0);
        HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,scanning,far,accumulating"),initial_box);

        int num_passes = 8.0;
        double evol_time = num_passes * velocity.upper() / trajectory.parameter_value("width").lower();
        HybridTime evolution_time(evol_time,num_passes*3);

    	system.substitute(RealParameter("y0",y0));

    	try {
    		HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    		z_values[i] = orbit.final().bounding_box()[5];
    	} catch (WorkingSetTooLargeException& ex) {
    		z_values[i] = Interval(-1);
    	}

        std::cout << "Final range for deposition: " << z_values[i] << std::flush << std::endl;

    	length += delta;
    }

    for (int i = 0; i < num_intervals; ++i) {
    	std::cout << z_values[i].lower() << " " << z_values[i].upper() << std::endl;
    }

}
