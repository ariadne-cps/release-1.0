/*****************************************************************************************************
 *            integration.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides a variant of the expression with integration.
 *
 *****************************************************************************************************/

#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    int VERBOSITY = 1;
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton trajectory("trajectory");

    // Parameters
    RealParameter velocity("velocity",1.0); // Velocity
    RealParameter half_width("half_width",1.0); // Half width

    /// Modes

    DiscreteLocation passing_right("passing_right");
    DiscreteLocation passing_left("passing_left");

	trajectory.new_mode(passing_right);
	trajectory.new_mode(passing_left);

    // Variables

    RealVariable x("x"); // Position of the laser

    trajectory.add_output_var(x);

    // Events

    DiscreteEvent switch_left("switch_right");
    DiscreteEvent switch_right("switch_left");

	// Dynamics

	RealExpression dyn_x_right = velocity;
	RealExpression dyn_x_left = -velocity;

	trajectory.set_dynamics(passing_right, x, dyn_x_right);
	trajectory.set_dynamics(passing_left, x, dyn_x_left);

	// Transitions

	RealExpression x_greater_half_width = x - half_width; // x >= half_width
	RealExpression x_lesser_minus_half_width = -x - half_width; // x <= -half_width

	std::map<RealVariable,RealExpression> reset_minus_half;
	reset_minus_half[x] = -half_width;
	std::map<RealVariable,RealExpression> reset_half;
	reset_half[x] = half_width;

	trajectory.new_forced_transition(switch_left,passing_right,passing_left,reset_half,x_greater_half_width);
	trajectory.new_forced_transition(switch_right,passing_left,passing_right,reset_minus_half,x_lesser_minus_half_width);

    /// Create a HybridAutomaton object
    HybridIOAutomaton measurements("measurements");

    RealParameter L("L",0.5);

    /// Create the discrete states
    DiscreteLocation far("far");
    DiscreteLocation close("close");

    RealVariable y("y");
    RealVariable z("z");

    DiscreteEvent enters("enters");
    DiscreteEvent exits("exits");

    measurements.add_input_var(x);
    measurements.add_output_var(y);
    //measurements.add_output_var(z);

    measurements.new_mode(far);
    measurements.new_mode(close);

	RealExpression dyn_y_far = 0.0;
	RealExpression dyn_y_close = 0.5*(Ariadne::cos(Ariadne::pi<Real>()*x*x/L/L)+1.0);
	measurements.set_dynamics(far, y, dyn_y_far);
	measurements.set_dynamics(close, y, dyn_y_close);

	RealExpression distance_greater_L = x*x - L*L; // distance >= L
	RealExpression distance_lesser_L = L*L - x*x; // distance <= L

	measurements.new_forced_transition(enters,far,close,distance_lesser_L);
	measurements.new_forced_transition(exits,close,far,distance_greater_L);

	HybridIOAutomaton system = compose("integration",trajectory,measurements,DiscreteLocation("passing_right"),DiscreteLocation("far"));

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(system);
    evolver.verbosity = VERBOSITY;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(2,2.0);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.0001;
    }

    Box initial_box(2, -half_width.value(),-half_width.value(), 0.0,0.0);
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("passing_right,far"),initial_box);

    HybridTime evolution_time(4.0,3);

    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(system.name());
    plotter.plot(orbit.reach(),"reach");
}
