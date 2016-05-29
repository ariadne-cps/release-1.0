/***************************************************************************
 *            home.h
 *
 *  Describes a system with a home environment with a controllable thermostat
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

#ifndef HOME_H_
#define HOME_H_

#include "ariadne.h"

namespace Ariadne {

HybridIOAutomaton getHomeSystem()
{
    /// Set the system parameters
	RealParameter day_len("day_len",1440); // Total time units to reach one day

    // System variables
    RealVariable h("h");    // Time reference
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
		RealExpression h_geq_allday = h - day_len; // Guard: h >= day
		timer.new_forced_transition(new_day, flow, flow, reset_h_zero, h_geq_allday);

    // Create the exterior temperature automaton

	    HybridIOAutomaton exterior("exterior");

	    // Parameters
		RealParameter Te_max("Te_max",16); // Maximum external daily temperature in °C
		RealParameter Te_min("Te_min",8); // Minimum external daily temperature in °C
		RealParameter phi("phi",-1.0); // Phase of the temperature oscillation

		// States
		DiscreteLocation oscillate("oscillate");

		// Add the input/output variables
		exterior.add_input_var(h);
		exterior.add_output_var(Te);

		// Dynamics
		RealExpression Te_dyn = (Te_max - Te_min) * Ariadne::pi<Real>()/day_len * sin(2.0*Ariadne::pi<Real>()/day_len*h + phi);
		exterior.new_mode(oscillate);
		exterior.set_dynamics(oscillate, Te, Te_dyn);

	// Create the thermostat controller automaton

		HybridIOAutomaton thermostat("thermostat");

		// Parameters
		RealParameter t_on_evening("t_on_evening",18.5); // Time to turn on the thermostat when evening arrives, expressed as hours from the start of the time
		RealParameter t_off_night("t_off_night",23.5); // Time to turn off the thermostat when night arrives, expressed as hours from the start of the time
		RealParameter t_on_morning("t_on_morning",6.5); // Time to turn on the thermostat when morning arrives, expressed as hours from the start of the time
		RealParameter t_off_day("t_off_day",8.5); // Time to turn off the thermostat when day arrives, expressed as hours from the start of the time

		RealVariable clk("clk"); // Thermostat clock

		// Add the input/output variables
		thermostat.add_internal_var(clk);

		// States
		DiscreteLocation day("day");
		DiscreteLocation evening("evening");
		DiscreteLocation night("night");
		DiscreteLocation morning("morning");

		// Add the modes
		thermostat.new_mode(day);
		thermostat.new_mode(evening);
		thermostat.new_mode(night);
		thermostat.new_mode(morning);

		// Add the dynamics
		RealExpression clk_d = 1.0;
		thermostat.set_dynamics(day,clk,clk_d);
		thermostat.set_dynamics(evening,clk,clk_d);
		thermostat.set_dynamics(night,clk,clk_d);
		thermostat.set_dynamics(morning,clk,clk_d);
		
		// Events
		DiscreteEvent turn_on_evening("turn_on_evening");
		DiscreteEvent turn_off_night("turn_off_night");
		DiscreteEvent turn_on_morning("turn_on_morning");
		DiscreteEvent turn_off_day("turn_off_day");

		// Add the input/output events
		thermostat.add_output_event(turn_on_evening);
		thermostat.add_output_event(turn_off_night);
		thermostat.add_output_event(turn_on_evening);
		thermostat.add_output_event(turn_off_day);

		// Transitions
		std::map<RealVariable,RealExpression> reset_clk;
		reset_clk[clk] = 0.0;
		RealExpression clk_geq_evening = clk - 60.0*(t_on_evening-t_off_day);
		thermostat.new_forced_transition(turn_on_evening, day, evening, reset_clk, clk_geq_evening);
		RealExpression clk_geq_night = clk - 60.0*(t_off_night-t_on_evening);
		thermostat.new_forced_transition(turn_off_night, evening, night, reset_clk, clk_geq_night);
		RealExpression t_night_morning_diff = t_on_morning - t_off_night;
		if (t_on_morning.value().midpoint() < t_off_night.value().midpoint())
			t_night_morning_diff = t_on_morning + day_len/60 - t_off_night;
		RealExpression clk_geq_morning = clk - 60.0*t_night_morning_diff;
		thermostat.new_forced_transition(turn_on_morning, night, morning, reset_clk, clk_geq_morning);
		RealExpression clk_geq_day = clk - 60.0*(t_off_day-t_on_morning);
		thermostat.new_forced_transition(turn_off_day, morning, day, reset_clk, clk_geq_day);

	/// Compose the automata
	HybridIOAutomaton timer_exterior = compose("timer-exterior",timer,exterior,flow,oscillate);
	HybridIOAutomaton system = compose("home",timer_exterior,thermostat,DiscreteLocation("flow,oscillate"),night);

	return system;
}


}

#endif /* HOME_H_ */
