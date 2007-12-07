/***************************************************************************
 *            python/export_vector_field_evolver.cc
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

#include "python/float.h"

#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/grid_set.h"
#include "geometry/list_set.h"

#include "system/vector_field.h"
#include "system/affine_vector_field.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/integrator_interface.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_vector_field_evolver() 
{
  typedef typename Evaluation::VectorFieldEvolver<R>::basic_set_type BS;

  class_< VectorFieldEvolver<R> > evolver_class("VectorFieldEvolver",init<const EvolutionParameters<R>&,const IntegratorInterface<BS>&>());
    evolver_class.def(init<const EvolutionParameters<R>&>());
    evolver_class.def("integrate",(SetInterface<R>*(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const SetInterface<R>&,const Rational&)const)
                      (&VectorFieldEvolver<R>::integrate),return_value_policy<manage_new_object>());
    evolver_class.def("reach",(SetInterface<R>*(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const SetInterface<R>&,const Rational&)const)
                      (&VectorFieldEvolver<R>::reach),return_value_policy<manage_new_object>());
    evolver_class.def("lower_reach",(SetInterface<R>*(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const SetInterface<R>&)const)
                      (&VectorFieldEvolver<R>::lower_reach),return_value_policy<manage_new_object>());
    evolver_class.def("chainreach",(SetInterface<R>*(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const SetInterface<R>&,const SetInterface<R>&)const)
                      (&VectorFieldEvolver<R>::chainreach),return_value_policy<manage_new_object>());
    evolver_class.def("viable",(SetInterface<R>*(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const SetInterface<R>&)const)
                      (&VectorFieldEvolver<R>::viable),return_value_policy<manage_new_object>());
    evolver_class.def("verify",(tribool(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const SetInterface<R>&,const SetInterface<R>&)const)
                      (&VectorFieldEvolver<R>::verify));
  

}

template void export_vector_field_evolver<FloatPy>();