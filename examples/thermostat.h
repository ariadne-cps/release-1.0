/***************************************************************************
 *            thermostat.h
 *
 *  Describes a system with a controllable thermostat in a closed environment
 *  subject to an outside temperature source.
 *
 *  Copyright  2016  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifndef THERMOSTAT_H_
#define THERMOSTAT_H_

#include "ariadne.h"

namespace Ariadne {

HybridIOAutomaton getThermostatSystem()
{
    /// Set the system parameters
	RealParameter day("day",1440); // Total time units to reach one day
	RealParameter Tmax("Tmax",16); // Maximum daily temperature in °C
	RealParameter Tmin("Tmin",8); // Minimum daily temperature in °C
	RealParameter phi("phi",-1.0); // Phase of the temperature oscillation

    // System variables
    RealVariable h("h");    // time
    RealVariable Te("Te"); // Exterior temperature

    // Create the day time automaton

    	HybridIOAutomaton timer("timer");

    	// States
    	DiscreteLocation flow("flow");

    	// Events
    	DiscreteEvent new_day("new_day");

		// Add the input/output variables
		timer.add_output_var(h);

		// Add the input/output events
		timer.add_internal_event(new_day);

    	// Dynamics
		RealExpression h_dyn = 1.0;
		timer.new_mode(flow);
		timer.set_dynamics(flow, h, h_dyn);

		// Transitions
		std::map<RealVariable,RealExpression> reset_h_zero;
		reset_h_zero[h] = 0.0;
		RealExpression h_geq_day = h - day; // Guard: h >= day
		timer.new_forced_transition(new_day, flow, flow, reset_h_zero, h_geq_day);

    // Create the exterior temperature automaton

	    HybridIOAutomaton exterior("exterior");

		// States
	    DiscreteLocation start("start");
		DiscreteLocation oscillate("oscillate");

		// Add the input/output variables
		exterior.add_input_var(h);
		exterior.add_output_var(Te);

    	// Events
    	DiscreteEvent init("init");

		// Add the input/output events
    	exterior.add_internal_event(init);

		// Dynamics
		RealExpression Te_dyn_start = 0.0;
		exterior.new_mode(start);
		exterior.set_dynamics(start, Te, Te_dyn_start);
		RealExpression Te_dyn_oscillate = (Tmax - Tmin) * Ariadne::pi<Real>()/day * sin(2.0*Ariadne::pi<Real>()/day*h + phi);
		exterior.new_mode(oscillate);
		exterior.set_dynamics(oscillate, Te, Te_dyn_oscillate);

		// Transitions
		std::map<RealVariable,RealExpression> reset_Te_init;
		reset_Te_init[Te] = Tmin + (Tmax - Tmin) * (1 - cos(phi)) / 2;
		RealExpression h_geq_zero = h; // Guard: h >= 0
		exterior.new_forced_transition(init, start, oscillate, reset_Te_init, h_geq_zero);

	/// Compose the automata
	HybridIOAutomaton system = compose("thermostat",timer,exterior,flow,start);

	return system;
}


}

#endif /* THERMOSTAT_H_ */
