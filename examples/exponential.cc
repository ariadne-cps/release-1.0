/*****************************************************************************************************
 *            exponential.cc
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of a CMOS inverter fed by a sinusoidal (thus analog) input.
 *
 *****************************************************************************************************/

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

int main(int argc, char* argv[])
{
    /// Dynamics parameters
    Vector<Float> dp(11);

    double A = 2.0;

    /// Constants
    float EVOL_TIME = 20.0/A;   /// Evolution time
    float MAX_ENCL_WIDTH = 0.1;   /// Maximum enclosure width
    float MAX_STEP_SIZE = 1e-2 / A;     /// Maximum step size
    int VERBOSITY = 1;              /// Verbosity of the HybridEvolver
	if (argc > 1)
		VERBOSITY = atoi(argv[1]);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton inverter;

    /// Create the discrete states
    DiscreteLocation work("work");

    RealVariable t("t");
    RealVariable x("x");

    inverter.add_output_var(t);
    inverter.add_output_var(x);

	inverter.new_mode(work);

	RealExpression dyn_t = 1.0;
	inverter.set_dynamics(work, t, dyn_t);
	RealExpression dyn_x = - A*(x + x*x/2);
	inverter.set_dynamics(work, x, dyn_x);

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(inverter);
    evolver.verbosity = VERBOSITY;

    evolver.settings().hybrid_maximum_step_size[work] = MAX_STEP_SIZE;
    evolver.settings().minimum_discretised_enclosure_widths[work] = Vector<Float>(2,MAX_ENCL_WIDTH);

    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;

    Box initial_box(2, 0.0,0.0, 1.0,1.0);
    HybridEnclosureType initial_enclosure(work,initial_box);

    HybridTime evolution_time(EVOL_TIME,1);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    Box graphic_box(2, 0.0, EVOL_TIME, -1.0, 1.0);

    plot("exponential", graphic_box, Colour(0.0,0.5,1.0), orbit);
}
