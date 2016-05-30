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

HybridIOAutomaton getSpringsAutomaton()
{
    /// Set the system parameters
	RealParameter m1("m1",4.0); // Mass of the first ball
	RealParameter m2("m2",0.75); // Mass of the second ball
    RealParameter k1("k1",2.0); // Elastic constant of the first spring
    RealParameter k2("k2",1.0); // Elastic constant of the second spring
    RealParameter p1("p1",1.0); // Neutral position for the first spring
    RealParameter p2("p2",2.0); // Neutral position for the second spring
    RealParameter st("st",1.9); // Stickyness

    /// Create a HybridAutomaton object
    HybridIOAutomaton automaton("springs");

    // Variables
    RealVariable x1("x1"); // Position of the first ball
    RealVariable x2("x2"); // Position of the second ball
    RealVariable v1("v1"); // Speed of the first ball
    RealVariable v2("v2"); // Speed of the second ball

    // Register variables
    automaton.add_internal_var(x1);
    automaton.add_internal_var(x2);
    automaton.add_internal_var(v1);
    automaton.add_internal_var(v2);

    /// Create the discrete states
    DiscreteLocation free("free");
    DiscreteLocation stuck("stuck");

    // Register modes
    automaton.new_mode(free);
    automaton.new_mode(stuck);

    /// Create the discrete events
    DiscreteEvent sticking("sticking");
    DiscreteEvent unsticking("unsticking");

    automaton.add_internal_event(sticking);
    automaton.add_internal_event(unsticking);

    /// Create the dynamics

    RealExpression free_x1_d = v1;
    RealExpression free_x2_d = v2;
    RealExpression free_v1_d = k1*(p1-x1)/m1;
    RealExpression free_v2_d = k2*(p2-x2)/m2;
    automaton.set_dynamics(free,x1,free_x1_d);
    automaton.set_dynamics(free,x2,free_x2_d);
    automaton.set_dynamics(free,v1,free_v1_d);
    automaton.set_dynamics(free,v2,free_v2_d);

    RealExpression stuck_x1_d = v1;
    RealExpression stuck_x2_d = v2;
    RealExpression stuck_v1_d = (k1*p1+k2*p2-(k1+k2)*x1)/(m1+m2);
    RealExpression stuck_v2_d = stuck_v1_d;
    automaton.set_dynamics(stuck,x1,stuck_x1_d);
    automaton.set_dynamics(stuck,x2,stuck_x2_d);
    automaton.set_dynamics(stuck,v1,stuck_v1_d);
    automaton.set_dynamics(stuck,v2,stuck_v2_d);

    /// Create the resets
    std::map<RealVariable,RealExpression> sticking_r;
    sticking_r[x1] = x1;
    sticking_r[x2] = x2;
    sticking_r[v1] = (m1*v1+m2*v2)/(m1+m2);
    sticking_r[v2] = sticking_r[v1];

    /// Create the guards
    RealExpression sticking_g = x1 - x2; // x1 >= x2
    RealExpression unsticking_g = (k1-k2)*x1 + k2*p2 - k1*p1 - st; // (k1-k2)*x1 + k2*p2 - k1*p1 >= st

    /// Transitions
    automaton.new_forced_transition(sticking,free,stuck,sticking_r,sticking_g);
    automaton.new_forced_transition(unsticking,stuck,free,unsticking_g);

	return automaton;
}


}

#endif /* SPRINGS_H_ */
