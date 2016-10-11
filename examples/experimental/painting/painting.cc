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
#include "../timeout.h"
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
    HybridIOAutomaton timeout = getTimeout();
    HybridIOAutomaton spray = getSpray();
    HybridIOAutomaton deposition = getDeposition();

    HybridIOAutomaton timer_traj = compose("timer-traj",timer,trajectory,DiscreteLocation("work"),DiscreteLocation("scanning"));
    HybridIOAutomaton timer_traj_timeout = compose("timer_traj_timeout",timer_traj,timeout,DiscreteLocation("work,scanning"),DiscreteLocation("running"));
    HybridIOAutomaton timer_traj_timeout_spray = compose("timer_traj_timeout_spray",timer_traj_timeout,spray,DiscreteLocation("work,scanning,running"),DiscreteLocation("far"));
    HybridIOAutomaton system = compose("painting",timer_traj_timeout_spray,deposition,DiscreteLocation("work,scanning,running,far"),DiscreteLocation("accumulating"));

    Real x0(0.02);
    Real velocity(0.05);
    //Real x0(0.000150603);

    system.substitute(RealParameter("x0",x0));
    Real angle = system.parameter_value("angle");
    Real vx_i = velocity*Ariadne::cos(angle);
    Real vy_i = velocity*Ariadne::sin(angle);
    Real x_i = 0.0;
    Real y_i = 0.0;

    /// Create a HybridEvolver object
    HybridEvolver evolver(system);
    evolver.verbosity = VERBOSITY;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(8,1.0);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.0005;
    }

    Box initial_box(8,
    		/*clk*/ 0.0,0.0,
			/*s*/ 0.0,0.0,
			/*t*/ 0.0,0.0,
			/*vx*/ vx_i.lower(),vx_i.upper(),
			/*vy*/ vy_i.lower(),vy_i.upper(),
			/*x*/ x_i.lower(),x_i.upper(),
			/*y*/ y_i.lower(),y_i.upper(),
			/*z*/ 0.0,0.0);
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,scanning,running,far,accumulating"),initial_box);

    double evol_time = 2.0*timeout.parameter_value("stop_time");
    HybridTime evolution_time(evol_time,6);

    //cout << system << endl;

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");

}
