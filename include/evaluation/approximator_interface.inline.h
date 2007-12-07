/***************************************************************************
 *            approximator_interface.inline.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/grid_cell_list_set.h"
#include "geometry/grid_mask_set.h"

namespace Ariadne {

template<class BS> 
void
Evaluation::ApproximatorInterface<BS>::adjoin_outer_approximation(Geometry::GridMaskSet<R>& gms, 
                                                                  const BS& bs) const
{
  gms.adjoin(this->outer_approximation(bs,gms.grid()));
}


template<class BS> 
void
Evaluation::ApproximatorInterface<BS>::adjoin_outer_approximation(Geometry::GridMaskSet<R>& gms, 
                                                                 const Geometry::ListSet<BS>& ls) const
{
  typedef Geometry::ListSet<BS> list_set_type;
  for(typename list_set_type::const_iterator bs_iter=ls.begin(); bs_iter!=ls.end(); ++bs_iter) {
    this->adjoin_outer_approximation(gms,*bs_iter);
  }
}

} // namespace Ariadne
