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

    HybridIOAutomaton timer_traj = compose("timer-traj",timer,trajectory,DiscreteLocation("work"),DiscreteLocation("scanning"));
    HybridIOAutomaton timer_traj_exp = compose("timer_traj_exposure",timer_traj,spray,DiscreteLocation("work,scanning"),DiscreteLocation("far"));
    HybridIOAutomaton timer_traj_exp_temp = compose("timer_traj_exp_temp",timer_traj_exp,skin_temperature,DiscreteLocation("work,scanning,far"),DiscreteLocation("varying"));
    HybridIOAutomaton system = compose("laser",timer_traj_exp_temp,cutting_depth,DiscreteLocation("work,scanning,far,varying"),DiscreteLocation("idle"));

    Real x0(0.0023);
    //Real x0(0.000150603);

    system.substitute(RealParameter("x0",x0));
    system.substitute(RealParameter("velocity",velocity));
    Real vx_i = -velocity;
    Real x_i = 4.0*spray.parameter_value("L") + x0;

    /// Create a HybridEvolver object
    HybridEvolver evolver(system);
    evolver.verbosity = VERBOSITY;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(7,0.5);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.000002;
    }

    Box initial_box(7, /*T*/ T0.lower(),T0.upper(), /*p*/ 0.0,0.0, /*t*/ 0.0,0.0, /*vx*/ vx_i.lower(),vx_i.upper(), /*x*/ x_i.lower(),x_i.upper(), /*z*/ 0.0,0.0, /*zi*/ 0.0,0.0);
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("work,scanning,far,varying,idle"),initial_box);

    int num_half_cycles = 1;
    double evol_time = -8.0*spray.parameter_value("L")/vx_i.upper();
    //HybridTime evolution_time(pass_period.upper()/4*num_half_cycles,5*num_half_cycles);
    HybridTime evolution_time(evol_time,7);

    //cout << system << endl;

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");

    std::cout << "Final z : " << 2.0*orbit.final().bounding_box()[5] << std::endl;

    double depth = orbit.reach().bounding_box()[5].upper();
    std::cout << "Depth of cut : " << depth << std::endl;
    std::cout << "Maximum value of zi : " << orbit.reach().bounding_box()[6].upper() << std::endl;
    std::cout << "Carbonization occurred? " << (orbit.reach().bounding_box()[6].upper() > cutting_depth.parameter_value("z_thr").upper() ? "yes" : "no") << std::endl;
}
