/***************************************************************************
 *            denotable_set.cc
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <fstream>
#include <iomanip>

#include "macros.h"
#include "exceptions.h"
#include "stlio.h"
#include "function_set.h"
#include "list_set.h"
#include "denotable_set.h"
#include "taylor_set.h"

#include "set_interface.h"


namespace Ariadne {


Box eps_codomain(const DenotableSetType& denotable_set, const Vector<Float> eps, const VectorFunction& func)
{
	ARIADNE_ASSERT(denotable_set.dimension() == func.argument_size());
	ARIADNE_ASSERT(denotable_set.dimension() == eps.size());

	Box result = Box::empty_box(func.result_size());

	for (DenotableSetType::const_iterator cell_it = denotable_set.begin(); cell_it != denotable_set.end(); ++cell_it) {
		Box cell_box = cell_it->box();
		cell_box.widen(eps);
		result = hull(result,func.evaluate(cell_box));
	}

	return result;
}

tribool covers(const DenotableSetType& covering_set, const DenotableSetType& covered_set, const Vector<Float>& eps)
{
	ARIADNE_ASSERT_MSG(covering_set.dimension() == covered_set.dimension(),"The two sets must have the same dimensions.");
	ARIADNE_ASSERT_MSG(covering_set.dimension() == eps.size(),"The vector eps must have the same dimension of the sets");

	tribool result = true;

	if (covering_set.bounding_box().disjoint(covered_set.bounding_box()))
		return false;

	for (DenotableSetType::const_iterator cell_it = covered_set.begin(); cell_it != covered_set.end(); ++cell_it) {
		Box cell_bx = cell_it->box();
		cell_bx.widen(eps);

		tribool is_superset = covering_set.superset(cell_bx);

		if (indeterminate(is_superset))
			result = indeterminate;
		else if (!possibly(is_superset))
			return false;
	}

	return result;
}

tribool inside(const DenotableSetType& covered_set, const DenotableSetType& covering_set, const Vector<Float>& eps, int accuracy)
{
	ARIADNE_ASSERT_MSG(covering_set.dimension() == covered_set.dimension(),"The two sets must have the same dimensions.");
	ARIADNE_ASSERT_MSG(covering_set.dimension() == eps.size(),"The vector eps must have the same dimension of the sets");

	tribool result = true;

	if (covering_set.bounding_box().disjoint(covered_set.bounding_box()))
		return false;

	DenotableSetType enlarged_covering_set = covering_set;
	for (DenotableSetType::const_iterator cell_it = covering_set.begin(); cell_it != covering_set.end(); ++cell_it) {
		Box cell_bx = cell_it->box();
		cell_bx.widen(eps);
		enlarged_covering_set.adjoin_outer_approximation(cell_bx,accuracy);
	}

	for (DenotableSetType::const_iterator cell_it = covered_set.begin(); cell_it != covered_set.end(); ++cell_it) {
		tribool is_superset = enlarged_covering_set.superset(cell_it->box());

		if (indeterminate(is_superset))
			result = indeterminate;
		else if (!possibly(is_superset))
			return false;
	}

	return result;
}


} // namespace Ariadne

