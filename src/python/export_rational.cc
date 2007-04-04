/***************************************************************************
 *            python/export_rational.cc
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

using namespace Ariadne;
using namespace Ariadne::Numeric;

#include <boost/python.hpp>
using namespace boost::python;


void export_rational() {
  class_<Rational>("Rational")
    .def(init<int,int>())
    .def(init<Integer,Integer>())
    .def(init<int>())
    .def(init<Integer>())
    .def(init<double>())
    .def(init<Float>())
    .def(init<Rational>())
    .def("__neg__", &neg<Rational,Rational>)
    .def("__add__", &add<Rational,Rational,int>)
    .def("__add__", &add<Rational,Rational,Integer>)
    .def("__add__", &add<Rational,Rational,double>)
    .def("__add__", &add<Rational,Rational,Float,Rational,Rational>)
    .def("__add__", &add<Rational,Rational,Rational>)
    .def("__radd__", &add<Rational,Rational,int>)
    .def("__radd__", &add<Rational,Rational,Integer>)
    .def("__radd__", &add<Rational,Rational,double>)
    .def("__radd__", &add<Rational,Rational,Float,Rational,Rational>)
    .def("__sub__", &sub<Rational,Rational,int>)
    .def("__sub__", &sub<Rational,Rational,Integer>)
    .def("__sub__", &sub<Rational,Rational,double>)
    .def("__sub__", &sub<Rational,Rational,Float,Rational,Rational>)
    .def("__sub__", &sub<Rational,Rational,Rational>)
    .def("__rsub__", &rsub<Rational,Rational,int>)
    .def("__rsub__", &rsub<Rational,Rational,Integer>)
    .def("__rsub__", &rsub<Rational,Rational,double>)
    .def("__rsub__", &rsub<Rational,Rational,Float,Rational,Rational>)
    .def("__mul__", &mul<Rational,Rational,int>)
    .def("__mul__", &mul<Rational,Rational,Integer>)
    .def("__mul__", &mul<Rational,Rational,double>)
    .def("__mul__", &mul<Rational,Rational,Float,Rational,Rational>)
    .def("__mul__", &mul<Rational,Rational,Rational>)
    .def("__rmul__", &mul<Rational,Rational,int>)
    .def("__rmul__", &mul<Rational,Rational,Integer>)
    .def("__rmul__", &mul<Rational,Rational,double>)
    .def("__rmul__", &mul<Rational,Rational,Float,Rational,Rational>)
    .def("__div__", &div<Rational,Rational,int>)
    .def("__div__", &div<Rational,Rational,Integer>)
    .def("__div__", &div<Rational,Rational,double>)
    .def("__div__", &div<Rational,Rational,Float,Rational,Rational>)
    .def("__div__", &div<Rational,Rational,Rational>)
    .def("__rdiv__", &rdiv<Rational,Rational,int>)
    .def("__rdiv__", &rdiv<Rational,Rational,Integer>)
    .def("__rdiv__", &rdiv<Rational,Rational,double>)
    .def("__rdiv__", &rdiv<Rational,Rational,Float,Rational,Rational>)
    .def("__eq__", &eq<bool,Rational,Rational>)
    .def("__eq__", &eq<bool,Rational,double>)
    .def("__eq__", &eq<bool,Rational,Integer>)
    .def("__eq__", &eq<bool,Rational,double>)
    .def("__ne__", &ne<bool,Rational,Rational>)
    .def("__ne__", &ne<bool,Rational,double>)
    .def("__ne__", &ne<bool,Rational,Integer>)
    .def("__ne__", &ne<bool,Rational,double>)
    .def("__lt__", &lt<bool,Rational,Rational>)
    .def("__lt__", &lt<bool,Rational,double>)
    .def("__lt__", &lt<bool,Rational,Integer>)
    .def("__lt__", &lt<bool,Rational,double>)
    .def("__gt__", &gt<bool,Rational,double>)
    .def("__gt__", &gt<bool,Rational,Integer>)
    .def("__gt__", &gt<bool,Rational,double>)
    .def("__le__", &le<bool,Rational,Rational>)
    .def("__le__", &le<bool,Rational,double>)
    .def("__le__", &le<bool,Rational,Integer>)
    .def("__le__", &le<bool,Rational,double>)
    .def("__ge__", &ge<bool,Rational,double>)
    .def("__ge__", &ge<bool,Rational,Integer>)
    .def("__ge__", &ge<bool,Rational,double>)
    .def("numerator", &Rational::numerator)
    .def("denominator", &Rational::denominator)
    .def(self_ns::str(self))
  ;
 
  
}