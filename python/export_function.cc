/***************************************************************************
 *            python/export_function.cc
 *
 *  21 October 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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

#include "function.h"

#include <boost/python.hpp>

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

#include "real_typedef.h"

void export_function() {
  def("exp_approx", Ariadne::exp_approx<Real>, "approximate exponential function (maximum error e)" );
  def("cos_approx", Ariadne::exp_approx<Real>, "approximate sine function (maximum error e)" );
  def("sin_approx", Ariadne::exp_approx<Real>, "approximate cosine function (maximum error e)" );
}