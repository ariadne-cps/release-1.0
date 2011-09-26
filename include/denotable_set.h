/***************************************************************************
 *            denotable_set.h
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

/*! \file denotable_set.h
 *  \brief Denotable set definition and utility functions.
 */

#ifndef DENOTABLE_SET_H
#define DENOTABLE_SET_H

#include "bdd_set.h"
#include "grid_set.h"

using namespace std;
using namespace Ariadne;

namespace Ariadne {

typedef BDDTreeSet DenotableSetType;
//typedef GridTreeSet DenotableSetType;

//! \brief Evaluates the codomain of \a func applied on the cells of \a denotable_set, each widened by \a eps.
Box eps_codomain(const DenotableSetType& denotable_set, const Vector<Float> eps, const VectorFunction& func);

//! \brief Check whether \a covering_set covers \a covered_set with a tolerance of \a eps.
//! \details Since the cell boxes of \a covered_set, enlarged of \a eps, are checked against \a covering_set,
//! the two sets can feature different grids.
tribool covers(const DenotableSetType& covering_set, const DenotableSetType& covered_set, const Vector<Float>& eps);

//! \brief Check whether \a covering_set covers \a covered_set with a tolerance of \a eps.
//! \details Since the cell boxes of \a covered_set are checked against an overapproximation (using \a accuracy) of the
//! epsilon-enlargement of \a covering_set, the two sets can feature different grids.
tribool inside(const DenotableSetType& covered_set, const DenotableSetType& covering_set, const Vector<Float>& eps, int accuracy);

}

#endif /* DENOTABLE_SET_H */
