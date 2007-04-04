/***************************************************************************
 *            python/export_float.cc
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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


#include "python/python_utilities.h"

#include "python/python_float.h"
#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/interval.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void 
export_float() 
{
  typedef R Float;
  typedef Interval<R> IFloat;

  class_<Float>("Float")
    .def(init<int>())
    .def(init<double>())
    .def(init<Float>())
    .def("__neg__", &neg<Float,Float>)
    .def("__add__", &add<IFloat,Float,int,Float,Float>)
    .def("__add__", &add<IFloat,Float,double,Float,Float>)
    .def("__add__", &add<IFloat,Float,Float>)
    .def("__radd__", &add<IFloat,Float,int,Float,Float>)
    .def("__radd__", &add<IFloat,Float,double,Float,Float>)
    .def("__sub__", &sub<IFloat,Float,int,Float,Float>)
    .def("__sub__", &sub<IFloat,Float,double,Float,Float>)
    .def("__sub__", &sub<IFloat,Float,Float>)
    .def("__rsub__", &rsub<IFloat,Float,int,Float,Float>)
    .def("__rsub__", &rsub<IFloat,Float,double,Float,Float>)
    .def("__mul__", &mul<IFloat,Float,int,Float,Float>)
    .def("__mul__", &mul<IFloat,Float,double,Float,Float>)
    .def("__mul__", &mul<IFloat,Float,Float>)
    .def("__rmul__", &mul<IFloat,Float,int,Float,Float>)
    .def("__rmul__", &mul<IFloat,Float,double,Float,Float>)
    .def("__div__", &div<IFloat,Float,int,Float,Float>)
    .def("__div__", &div<IFloat,Float,double,Float,Float>)
    .def("__div__", &div<IFloat,Float,Float>)
    .def("__rdiv__", &rdiv<IFloat,Float,int,Float,Float>)
    .def("__rdiv__", &rdiv<IFloat,Float,double,Float,Float>)
    .def("__eq__", &eq<bool,Float,double>)
    .def("__eq__", &eq<bool,Float,Float>)
    .def("__ne__", &ne<bool,Float,double>)
    .def("__ne__", &ne<bool,Float,Float>)
    .def("__lt__", &lt<bool,Float,double>)
    .def("__lt__", &lt<bool,Float,Float>)
    .def("__gt__", &gt<bool,Float,double>)
    .def("__gt__", &gt<bool,Float,Float>)
    .def("__le__", &le<bool,Float,double>)
    .def("__le__", &le<bool,Float,Float>)
    .def("__ge__", &ge<bool,Float,double>)
    .def("__ge__", &ge<bool,Float,Float>)
    .def(self_ns::str(self))
  ;
  
}

template void export_float<Float>();