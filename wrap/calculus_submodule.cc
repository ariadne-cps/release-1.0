/***************************************************************************
 *            calculus_submodule.cc
 *
 *  Copyright 2008  Pieter Collins
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#include "taylor_variable.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;

#include "utilities.h"

namespace Ariadne {

TaylorVariable*
make_taylor_variable(const uint& as, const uint& d, const Interval& e, const boost::python::object& obj) 
{
  array<Float> data;
  read_array(data,obj);
  //assert(data.size()==compute_polynomial_data_size(1u,as,d));
  const Float* ptr=data.begin();
  return new TaylorVariable(as,d,e,ptr);
}

boost::python::list
make_taylor_variables(const Vector<Interval>& x) 
{
  boost::python::list result;
  for(uint i=0; i!=x.size(); ++i) {
    result.append(TaylorVariable::variable(x.size(),midpoint(x[i]),i));
  }
  return result;
}


template<class C, class I, class X> inline 
void set_item(C& c, const I& i, const X& x) {
  c[i]=x;
}


template<class C, class I> inline 
typename C::ValueType 
get_item(const C& c, const I& i) {
  return c[i];
}

template<> std::string __repr__(const TaylorVariable& tv) {
  std::stringstream ss;
  ss << "TV(" << tv.expansion() << "," << tv.error() << ")";
  return ss.str();
} 

} // namespace Ariadne

void export_taylor_variable() 
{
  typedef double D;
  typedef Float R;
  typedef Interval I;
  typedef MultiIndex A;
  typedef Vector<Float> V;
  typedef TaylorVariable T;


  class_<T> taylor_variable_class("TaylorVariable");
  taylor_variable_class.def("__init__", make_constructor(&make_taylor_variable) );
  taylor_variable_class.def( init< uint >());
  taylor_variable_class.def( init< uint, uint >());
  taylor_variable_class.def("error", (const I&(T::*)()const) &T::error, return_value_policy<copy_const_reference>());
  taylor_variable_class.def("__getitem__", &get_item<T,A>);
  taylor_variable_class.def("__setitem__",&set_item<T,A,D>);
  taylor_variable_class.def("__setitem__",&set_item<T,A,R>);
  taylor_variable_class.def(-self);
  taylor_variable_class.def(self+self);
  taylor_variable_class.def(self-self);
  taylor_variable_class.def(self*self);
  taylor_variable_class.def(self/self);
  taylor_variable_class.def(self+R());
  taylor_variable_class.def(self-R());
  taylor_variable_class.def(self+=R());
  taylor_variable_class.def(self-=R());
  taylor_variable_class.def(self*=R());
  taylor_variable_class.def(self/=R());
  taylor_variable_class.def(self_ns::str(self));
  taylor_variable_class.def("__repr__",&__repr__<T>);

  taylor_variable_class.def("variables",&make_taylor_variables);
  taylor_variable_class.staticmethod("variables");

  def("constant",(T(*)(uint, const R&))&T::constant);
  def("variable",(T(*)(uint, const R&, uint))&T::variable);
  

  def("max",(T(*)(const T&,const T&))&max);
  def("min",(T(*)(const T&,const T&))&min);
  def("abs",(T(*)(const T&))&abs);
  def("neg",(T(*)(const T&))&neg);
  def("rec",(T(*)(const T&))&rec);
  def("pow",(T(*)(const T&, int))&pow);

  def("sqrt", (T(*)(const T&))&sqrt);
  def("exp", (T(*)(const T&))&exp);
  def("log", (T(*)(const T&))&log);
  /*
  def("sin", (T(*)(const T&))&sin);
  def("cos", (T(*)(const T&))&cos);
  def("tan", (T(*)(const T&))&tan);
  def("asin", (T(*)(const T&))&asin);
  def("acos", (T(*)(const T&))&acos);
  def("atan", (T(*)(const T&))&atan);
  */
}

void calculus_submodule() 
{
  export_taylor_variable();
}


