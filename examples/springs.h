/***************************************************************************
 *            springs.h
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

#ifndef SPRINGS_H_
#define SPRINGS_H_

#include "ariadne.h"

namespace Ariadne {

HybridIOAutomaton getSpringsSystem()
{
    /// Set the system parameters
	RealParameter m1("m1",4.0); // Mass of the first ball
	RealParameter m2("m2",0.75); // Mass of the second ball
    RealParameter k1("k1",2.0); // Elastic constant of the first spring
    RealParameter k2("k2",1.0); // Elastic constant of the second spring
    RealParameter p1("p1",1.0); // Neutral position for the first spring
    RealParameter p2("p2",2.0); // Neutral position for the second spring
    RealParameter st("st",1.9); // Stickyness

	// Variables
	RealVariable x1("x1"); // Position of the first ball
	RealVariable v1("v1"); // Speed of the first ball
	RealVariable x2("x2"); // Position of the second ball
	RealVariable v2("v2"); // Speed of the second ball

    /// Spring 1
    HybridIOAutomaton spring1("spring1");

    	// Register variables
    	spring1.add_input_var(x2);
    	spring1.add_input_var(v2);
    	spring1.add_internal_var(x1);
    	spring1.add_output_var(v1);

		/// Create the discrete states
		DiscreteLocation free1("free1");
		DiscreteLocation stuck1("stuck1");

		// Register modes
		spring1.new_mode(free1);
    	spring1.new_mode(stuck1);

    	/// Create the discrete events
    	DiscreteEvent sticking("sticking");
    	DiscreteEvent unsticking("unsticking");

    	// Register events
    	spring1.add_output_event(sticking);
    	spring1.add_output_event(unsticking);

    	/// Create the dynamics
    	RealExpression free_x1_d = v1;
    	RealExpression free_v1_d = k1*(p1-x1)/m1;
    	spring1.set_dynamics(free1,x1,free_x1_d);
    	spring1.set_dynamics(free1,v1,free_v1_d);
        RealExpression stuck_x1_d = v1;
        RealExpression stuck_v1_d = (k1*p1+k2*p2-(k1+k2)*x1)/(m1+m2);
        spring1.set_dynamics(stuck1,x1,stuck_x1_d);
        spring1.set_dynamics(stuck1,v1,stuck_v1_d);

        /// Create the guards
        RealExpression sticking_g = x1 - x2; // x1 >= x2
        RealExpression unsticking_g = (k1-k2)*x1 + k2*p2 - k1*p1 - st; // (k1-k2)*x1 + k2*p2 - k1*p1 >= st

        /// Create the resets
        std::map<RealVariable,RealExpression> sticking1_r;
        sticking1_r[x1] = x1;
        sticking1_r[v1] = (m1*v1+m2*v2)/(m1+m2);

        /// Transitions
        spring1.new_forced_transition(sticking,free1,stuck1,sticking1_r,sticking_g);
        spring1.new_forced_transition(unsticking,stuck1,free1,unsticking_g);


    /// Spring 2
	HybridIOAutomaton spring2("spring2");

    	// Register variable
		spring2.add_input_var(v1);
		spring2.add_output_var(x2);
		spring2.add_output_var(v2);

		/// Create the discrete states
		DiscreteLocation free2("free2");
		DiscreteLocation stuck2("stuck2");

		// Register modes
		spring2.new_mode(free2);
    	spring2.new_mode(stuck2);

    	// Register events
    	spring2.add_input_event(sticking);
    	spring2.add_input_event(unsticking);

    	// Create the dynamics
    	RealExpression free_x2_d = v2;
    	RealExpression free_v2_d = k2*(p2-x2)/m2;
    	spring2.set_dynamics(free2,x2,free_x2_d);
    	spring2.set_dynamics(free2,v2,free_v2_d);
    	RealExpression stuck_x2_d = v2;
    	RealExpression stuck_v2_d = (k1*p1+k2*p2-(k1+k2)*x2)/(m1+m2);
    	spring2.set_dynamics(stuck2,x2,stuck_x2_d);
    	spring2.set_dynamics(stuck2,v2,stuck_v2_d);

		/// Create the resets
		std::map<RealVariable,RealExpression> sticking2_r;
		sticking2_r[x2] = x2;
		sticking2_r[v2] = (m1*v1+m2*v2)/(m1+m2);

		/// Transitions
		spring2.new_unforced_transition(sticking,free2,stuck2,sticking2_r);
		spring2.new_unforced_transition(unsticking,stuck2,free2);

    HybridIOAutomaton system = compose("springs",spring1,spring2,free1,free2);

	return system;
}


}

#endif /* SPRINGS_H_ */
