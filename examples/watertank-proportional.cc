/***************************************************************************
 *            watertank-proportional.cc
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

#include "ariadne.h"
#include "watertank-proportional.h"

using namespace Ariadne;

int main(int argc,char *argv[]) 
{
	int verb = 1;
	if (argc > 1)
		verb = atoi(argv[1]);

	// The system
	HybridIOAutomaton system = Ariadne::getWatertankProportional();

    HybridEvolver evolver(system);
    evolver.verbosity = verb;

    evolver.settings().set_reference_enclosure_widths(0.004);
    evolver.settings().set_maximum_step_size(0.0005);

    cout << system << endl;

    HybridEvolver::EnclosureType initial_enclosure(DiscreteLocation("stabilising,flow"),Box(2, 0.6,0.6, 6.72,6.72));

    HybridTime evol_limits(80.0,5);

    HybridEvolver::OrbitType orbit = evolver.orbit(initial_enclosure,evol_limits,UPPER_SEMANTICS);

    PlotHelper plotter(system);
    plotter.plot(orbit.reach(),"reach");
}
