/***************************************************************************
 *            boost-converter.cc
 *
 *  Copyright  2014  Luca Geretti
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

#include "ariadne.h"
#include "function.h"
#include "boost-converter.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int verb = 1;
	if (argc > 1)
		verb = atoi(argv[1]);

    bool plot_results = true;

	// The system
	HybridAutomaton system = Ariadne::getBoostConverter();

    HybridEvolver evolver(system);
    evolver.verbosity = verb;

    HybridSpace hspace(system.state_space());
    for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
        evolver.settings().minimum_discretised_enclosure_widths[hs_it->first] = Vector<Float>(3,1.0);
        evolver.settings().hybrid_maximum_step_size[hs_it->first] = 0.001;
    }
    
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("incr"),Box(3, 0.0,0.0, 1.173,1.173, 5.65,5.65));
  
    HybridTime evol_limits(20.0,14);
 
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evol_limits,UPPER_SEMANTICS);

    if (plot_results) {
        PlotHelper plotter(system.name());
        plotter.plot(orbit.reach(),"reach");
    }
}
