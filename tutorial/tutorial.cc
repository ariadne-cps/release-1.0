/***************************************************************************
 *            tutorial.cc
 *
 *  Provides a documented tutorial of system definition and analysis.
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

#include <ariadne.h> // Library header
#include "system.h" // System definition
#include "analysis.h" // Custom analysis routines to be run

// Specifies that classes from the Ariadne library can be used without the Ariadne:: prefix
using namespace Ariadne;

int main(int argc,char *argv[])
{
	// This snippet reads, from the first argument of the executable, the verbosity value to be used
	int verb = 1;
	if (argc > 1)
		verb = atoi(argv[1]);

	// Instructs not to produce any plot results.
	// Set this to true to create plots within a folder named 'tutorial-png' in the current working director
    bool plot_results = true;

    // Loads the system from the system.h file
    HybridIOAutomaton system = Ariadne::getSystem();

    // Outputs the system
    cout << system << endl;

    // Constructs an initial state, in particular from two different locations of the system
    // Please note how the system variables are ordered alphabetically: this is important to
    // understand this when we specify sets, in order to avoid dimension mismatches.
    HybridBoundedConstraintSet initial_set(system.state_space());
    initial_set[DiscreteLocation("flow,idle,rising")] = Box(2, 1.0,1.0, 6.0,7.5);
    initial_set[DiscreteLocation("flow,idle,falling")] = Box(2, 0.0,0.0, 6.0,7.5);

    // Runs the analysis routines set in the analysis.h file
    analyse(system,initial_set,verb,plot_results);
}
