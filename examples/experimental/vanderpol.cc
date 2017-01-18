/***************************************************************************
 *            vanderpol.cc
 *
 *  Copyright  2017  Luca Geretti
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
    /// Create a HybridAutomton object
    HybridIOAutomaton vanderpol("vanderpol");

    /// Create four discrete states
    DiscreteLocation loc("loc");

    RealParameter mu("mu",1.0);

    RealVariable x("x"), y("y");

    /// Build the automaton
    vanderpol.add_internal_var(x);
    vanderpol.add_internal_var(y);

    RealExpression x_d = y;
    RealExpression y_d = mu * (1.0 - Ariadne::sqr(x))*y - x;

    vanderpol.new_mode(loc);
    vanderpol.set_dynamics(loc,x,x_d);
    vanderpol.set_dynamics(loc,y,y_d);

    /// Finished building the automaton

    cout << "Automaton = " << vanderpol << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver(vanderpol);

    /// Set the evolution parameters
    evolver.settings().set_reference_enclosure_widths(1e-2);
    evolver.settings().set_maximum_enclosure_widths_ratio(4000.0);
    evolver.settings().set_fixed_maximum_step_size(1e-2);
    evolver.settings().set_enable_adaptive_maximum_step_size(false);
    evolver.settings().set_enable_boxing_on_contraction(true);
    evolver.verbosity = 1;
    std::cout <<  evolver.settings() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;

    Box initial_box(2, 2.0,2.0, 0.0,0.0);
    HybridEnclosureType initial_enclosure(loc,initial_box);

    HybridTime evolution_time(10.0,4);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    PlotHelper plotter(vanderpol.name());
    plotter.plot(orbit.reach(),"reach");
}