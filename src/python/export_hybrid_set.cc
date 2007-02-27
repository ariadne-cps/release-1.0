/***************************************************************************
 *            python/export_hybrid_set.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <map>

#include "real_typedef.h"

#include "geometry/hybrid_set.h"
#include "geometry/arbitrary_set.h"

using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;
return_value_policy<reference_existing_object> return_reference_existing_object;

template<class HS, class T> 
inline void hybrid_set_new_location(HS& s, location_type q, const T& t) {
  s.new_location(q,t);
}

template<class HS, class A> 
inline void hybrid_set_set_item(HS& hs, location_type id, const A& x) {
  hs[id]=x;
}

template<class HS> inline 
const typename HS::set_type& hybrid_set_get_item(HS& hs, location_type id) {
  return hs[id];
}

template<class HS, class S> inline 
void hybrid_set_adjoin_set(HS& hs, location_type id, const S& s) {
  hs.adjoin(id,s);
}

template<class HS1, class HS2> inline 
void hybrid_set_adjoin_hybrid_set(HS1& hs1, const HS2& hs2) {
  hs1.adjoin(hs2);
}


template<class R>
void export_hybrid_set() 
{
 
  class_< std::map<location_type,dimension_type> >("DiscreteLocations",no_init)
    .def(self_ns::str(self))
  ;


  class_< HybridSet< ArbitrarySet<R> > >("HybridSet",init<>())
    .def("__len__",&HybridSet< ArbitrarySet<R> >::number_of_locations)
    .def("locations",&HybridSet< ArbitrarySet<R> >::locations)
    .def("__getitem__", &hybrid_set_get_item< HybridSet< ArbitrarySet<R> > >, return_reference_existing_object)
    .def(self_ns::str(self))
  ;
  
  class_< HybridGridMaskSet<R> >("HybridGridMaskSet",init<>())
    .def(init< HybridGridMaskSet<R> >())
    .def("__len__",&HybridGridMaskSet<R>::number_of_locations)
    .def("locations",&HybridGridMaskSet<R>::locations)
    .def("new_location",&hybrid_set_new_location< HybridGridMaskSet<R>, FiniteGrid<R> >)
    .def("new_location",&hybrid_set_new_location< HybridGridMaskSet<R>, GridMaskSet<R> >)
    .def("__getitem__", &hybrid_set_get_item< HybridGridMaskSet<R> >, return_reference_existing_object)
    .def("__setitem__", &hybrid_set_set_item< HybridGridMaskSet<R>, GridMaskSet<R> >)
    .def("adjoin", &hybrid_set_adjoin_set< HybridGridMaskSet<R>, GridCellListSet<R> >)
    .def("adjoin", &hybrid_set_adjoin_set< HybridGridMaskSet<R>, GridMaskSet<R> >)
    .def(self_ns::str(self))
  ;
  
  class_< HybridGridCellListSet<R> >("HybridGridCellListSet",init<>())
    .def(init< HybridGridMaskSet<R> >())
    .def(init< HybridGridCellListSet<R> >())
    .def("__len__",&HybridGridCellListSet<R>::number_of_locations)
    .def("locations",&HybridGridCellListSet<R>::locations)
    .def("new_location",&hybrid_set_new_location< HybridGridCellListSet<R>, Grid<R> >)
    .def("new_location",&hybrid_set_new_location< HybridGridCellListSet<R>, GridCellListSet<R> >)
    .def("new_location",&hybrid_set_new_location< HybridGridCellListSet<R>, GridMaskSet<R> >)
    .def("__getitem__",&hybrid_set_get_item< HybridGridCellListSet<R> >, return_reference_existing_object)
    .def("__setitem__",&hybrid_set_set_item< HybridGridCellListSet<R>, GridCellListSet<R> >)
    .def("adjoin", &hybrid_set_adjoin_set< HybridGridCellListSet<R>, GridCellListSet<R> >)
    .def(self_ns::str(self))
  ;
  
}

template void export_hybrid_set<Real>();
