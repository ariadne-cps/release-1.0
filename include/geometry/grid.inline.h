/***************************************************************************
 *            grid.inline.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 

namespace Ariadne {
  namespace Geometry {

    template<class R> inline
    const Grid<R>& 
    FiniteGrid<R>::grid() const
    { 
      return *this->_grid_ptr; 
    }
  
    template<class R> inline
    const Combinatoric::LatticeBlock& 
    FiniteGrid<R>::lattice_block() const 
    {
      return this->_lattice_block; 
    }
     

    template<class R> inline
    R 
    FiniteGrid<R>::subdivision_coordinate(dimension_type d, index_type n) const
    {
      return this->_grid_ptr->subdivision_coordinate(d,n);
    }
  
    template<class R> inline
    index_type 
    FiniteGrid<R>::subdivision_interval(dimension_type d, const real_type& x) const 
    {
      return this->_grid_ptr->subdivision_interval(d,x);
    }
      
    template<class R> inline
    index_type 
    FiniteGrid<R>::subdivision_lower_index(dimension_type d, const real_type& x) const 
    {
      return this->_grid_ptr->subdivision_lower_index(d,x);
    }

    template<class R> inline
    index_type 
    FiniteGrid<R>::subdivision_upper_index(dimension_type d, const real_type& x) const 
    {
      return this->_grid_ptr->subdivision_upper_index(d,x);
    }
  




    template<class R> inline
    std::istream& operator>>(std::istream& is, Grid<R>& g)
    {
      return g.read(is);
    }

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Grid<R>& g)
    {
      return g.write(os);
    }

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const FiniteGrid<R>& fg)
    {
      return fg.write(os);
    }

  }
}