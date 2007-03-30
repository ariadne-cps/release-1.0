/***************************************************************************
 *            irregular_grid.cc
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include "geometry/irregular_grid.h"
#include "geometry/irregular_grid.code.h"

#include "numeric/float.h"

namespace Ariadne {
  namespace Geometry {

    template class IrregularGrid<Rational>;

#ifdef ENABLE_FLOAT64
    template class IrregularGrid<Float64>;
#endif
  
#ifdef ENABLE_FLOATMP
    template class IrregularGrid<FloatMP>;
#endif

  }
}
