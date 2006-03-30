/***************************************************************************
 *            simplex.tpl
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
#include <iostream>
#include <vector>

#include "../geometry/simplex.h"

namespace Ariadne {
  namespace Geometry {

    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Simplex<R>& s) 
    {
//      if(s.empty()) {
//        os << "Empty";
//     }
//      else 
      if(s.dimension() > 0) {
        os << "Simplex(vertices=" << s.vertices() << ") ";
      }

      return os;
    }
    
    template <typename R>
    std::istream& 
    operator>>(std::istream& is, Simplex<R>& s)
    {
      throw std::domain_error("Not implemented");
    }
      
  }
}