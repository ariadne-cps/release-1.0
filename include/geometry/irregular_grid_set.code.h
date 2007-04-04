/***************************************************************************
 *            irregular_grid_set.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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


#include "../base/stlio.h"

#include "rectangle.h"
#include "list_set.h"

namespace Ariadne {
  namespace Geometry {



    template<class R>
    IrregularGridMaskSet<R>::IrregularGridMaskSet(const ListSet< Rectangle<R> >& rls) 
      : _grid(IrregularGrid<R>(rls)), 
        _lattice_set(_grid.block())
    {
      for(typename ListSet< Rectangle<R> >::const_iterator r_iter=rls.begin(); 
          r_iter!=rls.end(); ++r_iter) 
      {
        this->_lattice_set.adjoin(this->_grid.index_block(*r_iter));
      }
    }    
   

 
    template<class R>
    dimension_type 
    IrregularGridMaskSet<R>::dimension() const
    {
      return this->_grid.dimension();
    }


    template<class R>
    const IrregularGrid<R>&
    IrregularGridMaskSet<R>::grid() const
    {
      return this->_grid;
    }


    template<class R>
    const Combinatoric::LatticeMaskSet&
    IrregularGridMaskSet<R>::lattice_set() const
    {
      return this->_lattice_set;
    }


  
    template<class R>
    IrregularGridMaskSet<R>::operator ListSet< Rectangle<R> > () const
    {
      ListSet< Rectangle<R> > result(this->dimension());
      for(Combinatoric::LatticeMaskSet::const_iterator lms_iter=this->_lattice_set.begin(); 
          lms_iter!=this->_lattice_set.end(); ++lms_iter) 
      {
        result.adjoin(this->_grid.rectangle(*lms_iter));
      }
      return result;
    }    
   
    template<class R> 
    tribool
    subset(const Rectangle<R>& r, const IrregularGridMaskSet<R>& igms)
    {
      return Combinatoric::subset(igms.grid().index_block(r),igms.lattice_set());
    }


    template<class R> 
    tribool
    subset(const ListSet< Rectangle<R> >& rls, const IrregularGridMaskSet<R>& igms)
    {
      for(typename ListSet< Rectangle<R> >::const_iterator rls_iter=rls.begin();
          rls_iter!=rls.end(); ++rls_iter)
      {        
        if(!subset(*rls_iter,igms)) {
          return false;
        }
      }
      return true;
    }


    template<class R>
    std::ostream&
    IrregularGridMaskSet<R>::write(std::ostream& os) const
    {
      os << "IrregularGridMaskSet( grid=" << this->_grid 
         << ", lattice_set=" << this->_lattice_set << " )";
      return os;
    }


    template<class R>
    void
    IrregularGridMaskSet<R>::_instantiate() 
    {
      Rectangle<R>* r=0;
      ListSet< Rectangle<R> >* rls=0;
      IrregularGridMaskSet<R>* igms=0;
      subset(*r,*igms);
      subset(*rls,*igms);
    }

  }
} 


