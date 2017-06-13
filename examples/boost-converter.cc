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
#include "boost-converter.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int verbosity = 1;
	if (argc > 1)
		verbosity = atoi(argv[1]);

    bool plot_results = true;

	// The system
	HybridIOAutomaton system = Ariadne::getBoostConverter();

	cout << system << endl;

    HybridEvolver evolver(system);
    evolver.verbosity = verbosity;

    evolver.settings().set_reference_enclosure_widths(1.0);
    evolver.settings().set_maximum_step_size(0.00001);
    
    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("incr,below_duty"),Box(3, 0.0,0.0, 1.0,1.0, 1.0,1.0));
  
    HybridTime evol_limits(18.0,30);
 
    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evol_limits,UPPER_SEMANTICS);

    if (plot_results) {
        PlotHelper plotter(system.name());
        plotter.plot(orbit.reach(),"reach");
    }
}
