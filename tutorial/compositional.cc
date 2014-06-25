/***************************************************************************
 *            compositional.cc
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
#include "analysis.h"
#include "compositional.h"

using namespace Ariadne;

int main(int argc,char *argv[])
{
	int verbosity = 1;
    bool plot_results = false;

	// The system
	HybridIOAutomaton system = Ariadne::getSystem();

	// The initial values
	HybridBoundedConstraintSet initial_set(system.state_space());
	initial_set[DiscreteLocation("flow,idle,rising")] = Box(2, 6.0,7.5, 1.0,1.0);
	initial_set[DiscreteLocation("flow,idle,falling")] = Box(2, 6.0,7.5, 0.0,0.0);

    analyse(system,initial_set,verbosity,plot_results);
}
