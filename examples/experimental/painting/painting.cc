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
    Real y0(0.0107,0.01075);

    system.substitute(RealParameter("y0",y0));
    system.substitute(RealParameter("d",d));
    Real vx_i = velocity;
    Real x_i = 0.0;
    Real y_i = 0.0;

    /// Create a HybridEvolver object
    HybridEvolver evolver(system);
    evolver.verbosity = VERBOSITY;

    evolver.settings().set_reference_enclosure_widths(0.7);
    evolver.settings().set_hybrid_maximum_step_size(0.01);

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

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Final range for deposition: " << orbit.final().bounding_box()[5] << std::flush << std::endl;

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");

}
