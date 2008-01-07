/***************************************************************************
 *            list_set.template.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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

#include "list_set.h"
#include "base/stlio.h"
#include "geometry/rectangle.h"
#include "geometry/name.h"
#include "geometry/irregular_grid_set.h"

namespace Ariadne {

template<class BS>
Geometry::ListSet<BS>*
Geometry::ListSet<BS>::clone() const
{
  return new ListSet<BS>(*this);
}


template<class BS> 
tribool
Geometry::ListSet<BS>::intersects(const Box<R>& r) const
{
  return !Geometry::disjoint(*this,r);
}


template<class BS> 
tribool
Geometry::ListSet<BS>::disjoint(const Box<R>& r) const
{
  return Geometry::disjoint(*this,r);
}


template<class BS>
tribool
Geometry::ListSet<BS>::superset(const Box<R>& r) const
{
  const ListSet< Box<R> >* rls=
    dynamic_cast< const ListSet< Box<R> > * >(this);
  if(rls) {
    return Geometry::subset(r,IrregularGridMaskSet<R>(*rls));
  } else {
    return indeterminate;
  }
}


template<class BS>
tribool
Geometry::ListSet<BS>::subset(const Box<R>& r) const
{
  return Geometry::subset(*this,r);
}


template<class BS1, class BS2>
tribool
Geometry::disjoint(const ListSet<BS1>& ls1,
                   const ListSet<BS2>& ls2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls1,ls2,"tribool disjoint(ListSet<BS> ls1, ListSet<BS> ls2)");
  tribool result=true;
  for (typename ListSet<BS1>::const_iterator i=ls1.begin(); i!=ls1.end(); ++i) {
    for (typename ListSet<BS2>::const_iterator j=ls2.begin(); j!=ls2.end(); ++j) {
      result = result && disjoint(*i,*j);
      if(!result) { return result; }
    }
  }
  return result;
}



template<class R>
tribool
Geometry::subset(const ListSet< Geometry::Box<R> >& bxls1,
                 const ListSet< Geometry::Box<R> >& bxls2)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class BS1, class BS2>
tribool
Geometry::subset(const ListSet< BS1 >& ls,
                 const BS2& bs)
{
  tribool result=true;
  for(typename ListSet< BS1 >::const_iterator ls_iter=ls.begin();
      ls_iter!=ls.end(); ++ls_iter)
    {
      result = result && Geometry::subset(*ls_iter,bs);
      if(result==false) {
        return result;
      }
    }
  return result;
}



template<class BS>
Geometry::ListSet<BS>
Geometry::join(const ListSet<BS>& ls1,
               const ListSet<BS>& ls2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls1,ls2,"ListSet<BS> join(ListSet<BS> ls1, ListSet<BS> ls2)");
  ListSet<BS> ds_union(ls1);
  ds_union.inplace_union(ls2);
  return ds_union;
}


template<class BS>
Geometry::ListSet<BS>
Geometry::open_intersection(const ListSet<BS>& ls1,
                            const ListSet<BS>& ls2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls1,ls2,"ListSet<BS> open_intersection(ListSet<BS> ls1, ListSet<BS> ls2)");
  ListSet<BS> ds_inter(ls1.dimension());
  for (size_type i=0; i<ls1.size(); i++) {
    for (size_type j=0; j<ls2.size(); j++) {
      if (!disjoint(ls1[i],ls2[j])) {
        ds_inter.push_back(open_intersection(ls1[i],ls2[j]));
      }
    }
  }
  return ds_inter;
}

template<class BS>
Geometry::ListSet<BS>
Geometry::inner_intersection(const ListSet<BS>& ls,
                             const SetInterface<typename BS::real_type>& ops)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls,ops,"ListSet<BS> inner_intersection(ListSet<BS> ls, ListSet<BS> ops)");
  ListSet<BS> ds(ls.dimension());
  for (size_type i=0; i<ls.size(); i++) {
    if(ops.subset(ls[i].bounding_box())) {
      ds.push_back(ls[i]);
    }
  }
  return ds;
}

template<class BS>
Geometry::ListSet<BS>
Geometry::lower_intersection(const ListSet<BS>& ls,
                             const SetInterface<typename BS::real_type>& cls)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls,cls,"ListSet<BS> inner_intersection(ListSet<BS> ls, ListSet<BS> cls)");
  ListSet<BS> ds(ls.dimension());
  for (size_type i=0; i<ls.size(); i++) {
    if(cls.intersects(ls[i].bounding_box())) {
      ds.push_back(ls[i]);
    }
  }
  return ds;
}

template<class BS>
Geometry::ListSet<BS>
Geometry::outer_intersection(const ListSet<BS>& ls,
                             const SetInterface<typename BS::real_type>& cps)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls,cps,"ListSet<BS> inner_intersection(ListSet<BS> ls, ListSet<BS> cps)");
  ListSet<BS> ds(ls.dimension());
  for (size_type i=0; i<ls.size(); i++) {
    if(cps.disjoint(ls[i].bounding_box())) {
    } else {
      ds.push_back(ls[i]);
    }
  }
  return ds;
}




template<class BS>
void
Geometry::ListSet<BS>::_instantiate()
{
}


template<class BS>
std::ostream& 
Geometry::ListSet<BS>::summarize(std::ostream& os) const
{
  const ListSet<BS>& ls=*this;
  os << "ListSet<"<<Geometry::name<BS>()<<">( size="<<ls.size()<<", dimension="<<ls.dimension()<<" )";
  return os;
}

template<class BS>
std::ostream& 
Geometry::ListSet<BS>::write(std::ostream& os) const
{
  const ListSet<BS>& ls=*this;
  os << "ListSet<"<<Geometry::name<BS>()<<">(\n size="<<ls.size()<<",\n";
  os << "  [ ";
  if (ls.size() >0 ) {
    os << ls[0];
  }
  for (size_type i=1; i<ls.size(); i++) {
    os << ",\n    " << ls[i];
  }
  os << "\n  ]\n";
  os << " }" << std::endl;
  return os;
}


template<class BS>
std::istream& 
Geometry::ListSet<BS>::read(std::istream& is)
{
  ListSet<BS>& ls=*this;
  std::vector< BS >& vec(ls._vector);
  is >> vec;
  
  if(vec.size()==0) {
    ls._dimension = 0;
  }
  else {
    ls._dimension=vec[0].dimension();
  }
  return is;
}



}