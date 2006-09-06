/***************************************************************************
 *            discrete_time_system.tpl
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it,  Pieter.Collins@cwi.nl
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
 
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/simplex.h"
#include "../geometry/polyhedron.h"

#include "../system/discrete_time_system.h"

namespace Ariadne {
  namespace System {

    template<typename R>
    DiscreteTimeSystem<R>::~DiscreteTimeSystem() 
    {
    }
  
    template<typename R>
    Geometry::Point<R> 
    DiscreteTimeSystem<R>::operator() (const Geometry::Point<R>& x,
                                       const Geometry::Point<R>& u,
                                       const Geometry::Point<R>& v) const
    {
      throw std::invalid_argument(this->name()+"::operator() (Point,Point,Point) not implemented."); 
    }
    
    template<typename R>
    Geometry::Rectangle<R> 
    DiscreteTimeSystem<R>::operator() (const Geometry::Rectangle<R>& x,
                                       const Geometry::Rectangle<R>& u,
                                       const Geometry::Rectangle<R>& v) const
    {
      throw std::invalid_argument(this->name()+"::operator() (Rectangle,Rectangle,Rectangle) not implemented."); 
    }
    
    template<typename R>
    Geometry::Zonotope<R> 
    DiscreteTimeSystem<R>::operator() (const Geometry::Zonotope<R>& x,
                                       const Geometry::Zonotope<R>& u,
                                       const Geometry::Zonotope<R>& v) const
    {
      throw std::invalid_argument(this->name()+"::operator() (Zonotope,Zonotope,Zonotope) not implemented."); 
    }
    
    
    template<typename R>
    LinearAlgebra::Matrix<R> 
    DiscreteTimeSystem<R>::derivative(const Geometry::Point<R>& x,
                                      const Geometry::Point<R>& u,
                                      const Geometry::Point<R>& v) const
    {
      throw std::invalid_argument(this->name()+"::derivative(Point,Point,Point) not implemented."); 
    }

    template<typename R>
    LinearAlgebra::Matrix<R> 
    DiscreteTimeSystem<R>::derivative(const Geometry::Rectangle<R>& x,
                                      const Geometry::Rectangle<R>& u,
                                      const Geometry::Rectangle<R>& v) const
    {
      throw std::invalid_argument(this->name()+"::derivative(Rectangle,Rectangle,Rectangle) not implemented."); 
    }