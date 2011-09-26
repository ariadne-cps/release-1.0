/***************************************************************************
 *            discretisation.cc
 *
 *  Copyright  2011  Luca Geretti
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
#include "denotable_set.h"

using namespace Ariadne;


int main(int argc, char* argv[])
{
	if (argc != 3) {
		cout << "Insert dimension and depth!" << "\n";
		return 0;
	}

	int dimension = atoi(argv[1]);
	int depth = atoi(argv[2]);

	Box bx(dimension,Interval(0.0,1.0));

	Grid gr(dimension);

	DenotableSetType ds(gr);

	time_t initial_time = time(NULL);

	ds.adjoin_over_approximation(bx,depth);

	time_t final_time = time(NULL);

	cout << "Took " << final_time - initial_time << " seconds. Look for memory consumption or press enter to end.\n";
	getchar();
}
