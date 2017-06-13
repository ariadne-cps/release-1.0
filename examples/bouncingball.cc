/***************************************************************************
 *            bouncingball.cc
 *
 *  Copyright  2008  Davide Bresolin
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

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;


int main()
{
    /// Set the system parameters
    RealParameter a("a",0.5); // coefficient of absorption
    RealParameter g("g",9.81); // gravity constant

    /// Build the Hybrid System

    HybridIOAutomaton ball("ball");

    DiscreteLocation freefall("freefall");

    /// Create the discrete events
    DiscreteEvent bounce("bounce");

    RealVariable x("x"), vx("vx");

    /// Build the automaton
    ball.add_internal_var(x);
    ball.add_internal_var(vx);

    RealExpression x_d = vx;
    RealExpression vx_d = -g;

    ball.new_mode(freefall);
    ball.set_dynamics(freefall,x,x_d);
    ball.set_dynamics(freefall,vx,vx_d);

    RealExpression guard = -x; // x <= 0

    std::map<RealVariable,RealExpression> reset;
    reset[x] = x;
    reset[vx] = -a * vx;

    ball.new_forced_transition(bounce,freefall,freefall,reset,guard);

    cout << "Automaton = " << ball << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(ball);

    /// Set the evolution parameters
    evolver.settings().set_reference_enclosure_widths(0.05);
    evolver.settings().set_maximum_step_size(1.0/64);
    evolver.verbosity = 1;
    std::cout <<  evolver.settings() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;

    Box initial_box(2, 0.0,0.001, 1.999,2.0);
    HybridEnclosureType initial_enclosure(freefall,initial_box);
    Box bounding_box(2, -10.1,10.1, -0.1,2.1);

    HybridTime evolution_time(4.0,4);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    plot("ball_vx-x",bounding_box, Colour(0.0,0.5,1.0), orbit);
}
