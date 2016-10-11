/*****************************************************************************************************
 *            cutting_depth.h
 *
 *  Copyright  2016  Luca Geretti
 *
 * Provides the behavior of the cutting depth under a laser.
 *
 *****************************************************************************************************/

#include "ariadne.h"

#ifndef CUTTING_DEPTH_H
#define CUTTING_DEPTH_H

namespace Ariadne {

HybridIOAutomaton getDeposition()
{
    /// Parameters
	RealParameter k("k",1.0e-3);
	RealParameter z_thr("z_thr",300e-6);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("deposition");

    /// Create the discrete states
    DiscreteLocation accumulating("accumulating");

    RealVariable s("s"); // The spray exposure
    RealVariable z("z"); // The deposition amount

    automaton.add_input_var(s);
    automaton.add_output_var(z);

	automaton.new_mode(accumulating);

	/// Dynamics

	RealExpression dyn_z = k*s;

	automaton.set_dynamics(accumulating, z, dyn_z);

	return automaton;

}

}

#endif
