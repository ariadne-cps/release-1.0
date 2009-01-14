/***************************************************************************
 *            solver_interface.h
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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
 
/*! \file solver_interface.h
 *  \brief Interface class for solving equations.
 */

#ifndef ARIADNE_SOLVER_INTERFACE_H
#define ARIADNE_SOLVER_INTERFACE_H

#include <exception>
#include <stdexcept>
#include <string>

namespace Ariadne {

    
class Interval;
template<class X> class Vector;
class FunctionInterface;

class EvaluationException : public std::runtime_error
{
  public:
    EvaluationException(const char* what) : std::runtime_error(what) { }
};

/*! \ingroup EvaluatorInterfaces \ingroup Solvers
 *  \brief %Interface for solving (nonlinear) equations. 
 */
class SolverInterface
{
  public:
    /*! \brief Virtual destructor. */
    virtual ~SolverInterface() { };
    
    /*! \brief The maximum permissible error of the solution. */
    virtual double maximum_error() const = 0;
    /*! \brief Set the maximum error. */
    virtual void set_maximum_error(double max_error) = 0;
    
    /*! \brief The maximum number of steps allowed before the method must quit. */
    virtual uint maximum_number_of_steps() const = 0;
    /*! \brief Set the maximum number of steps. */
    virtual void set_maximum_number_of_steps(uint max_steps) = 0;
    
    /*! \brief Solve \f$f(x)=0\f$, starting in the interval point \a pt. */
    virtual Vector<Interval> solve(const FunctionInterface& f,const Vector<Interval>& pt) = 0;
    /*! \brief Solve \f$f(x)=x\f$, starting in the interval point \a pt. */
    virtual Vector<Interval> fixed_point(const FunctionInterface& f,const Vector<Interval>& pt) = 0;
};
    

} // namespace Ariadne

#endif /* ARIADNE_SOLVER_INTERFACE_H */
