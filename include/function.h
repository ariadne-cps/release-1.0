/***************************************************************************
 *            function.h
 *
 *  Copyright 2008-9  Pieter Collins
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

/*! \file function.h
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_H
#define ARIADNE_FUNCTION_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function_interface.h"

#include "macros.h"
#include "pointer.h"
#include "container.h"

#include "vector.h"
#include "matrix.h"
#include "polynomial.h"
#include "affine.h"
#include "taylor_model.h"
#include "differential.h"
#include "real.h"

namespace Ariadne {

class Real;
template<class X> class Vector;
template<class X> class Polynomial;
template<class LHS, class RHS> class Assignment;
template<class T> class ExtendedVariable;
template<class T> class Variable;
template<class T> class Constant;
template<class T> class DottedVariable;
template<class T> class PrimedVariable;
template<class T> class Expression;
template<class T> class Space;
template<class T> class ConstantSetComparator;

typedef ExtendedVariable<Real> ExtendedRealVariable;
typedef Variable<Real> RealVariable;
typedef Constant<Real> RealConstant;
typedef Constant<Real> RealParameter;
typedef DottedVariable<Real> DottedRealVariable;
typedef PrimedVariable<Real> PrimedRealVariable;
typedef Expression<Real> RealExpression;
typedef Space<Real> RealSpace;

typedef Assignment<RealVariable,RealExpression> RealAssignment;
typedef Assignment<ExtendedRealVariable,RealExpression> ExtendedRealAssignment;
typedef Assignment<DottedRealVariable,RealExpression> DottedRealAssignment;
typedef Assignment<PrimedRealVariable,RealExpression> PrimedRealAssignment;
typedef std::set<Constant<Real>,ConstantSetComparator<Real> > RealParameterSet;

class ScalarFunction;
class VectorFunction;

//! A scalar function \f$f:\R^n\rightarrow\R\f$.
class ScalarFunction
{
    typedef uint Nat;
    typedef uint SizeType;
    typedef std::ostream OStream;
  public:
    static ScalarFunction constant(Nat n, double c);
    static ScalarFunction constant(Nat n, Real c);
    static ScalarFunction variable(Nat n, uint i);
    static ScalarFunction coordinate(Nat n, uint i);

    explicit ScalarFunction(Nat n=0u);
    ScalarFunction(const Expression<Real>& e, const Space<Real>& s);
    ScalarFunction(const Expression<Real>& e, const List< Variable<Real> >& s);
    ScalarFunction(const Expression<tribool>& e, const List< Variable<Real> >& s);
    ScalarFunction(const Polynomial<Real>& p);

    ScalarFunction(ScalarFunctionInterface* fptr) : _ptr(fptr) { }
    const ScalarFunctionInterface* pointer() const { return this->_ptr.operator->(); }
    operator const ScalarFunctionInterface& () const { return *this->_ptr; }

    ScalarFunction& operator=(const ScalarFunction& f) {
        this->_check_type(f.pointer()); this->_ptr=f._ptr; return *this; }
    ScalarFunction& operator=(int n) { (*this)=constant(0, n); return *this; }

    Nat argument_size() const { return this->_ptr->argument_size(); }
    Float evaluate(const Vector<Float>& x) const { return this->_ptr->evaluate(x); }
    Interval evaluate(const Vector<Interval>& x) const { return this->_ptr->evaluate(x); }
    TaylorModel evaluate(const Vector<TaylorModel>& x) const { return this->_ptr->evaluate(x); }

	/*! \brief Substitute the constant \a con, if present, on all expressions. */
	void substitute(const Constant<Real>& con);
	/*! \brief Get the parameters (i.e. the Constant<Real> whose name starts with a letter) from the function. */
	RealParameterSet parameters() const;

    Differential<Float> evaluate(const Vector< Differential<Float> >& x) const { return this->_ptr->evaluate(x); }
    Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const { return this->_ptr->evaluate(x); }

    Float operator() (const Vector<Float>& x) const { return this->_ptr->evaluate(x); }
    Interval operator() (const Vector<Interval>& x) const { return this->_ptr->evaluate(x); }

    Vector<Float> gradient(const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).gradient(); }
    Vector<Interval> gradient(const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).gradient(); }

    std::ostream& write(std::ostream& os) const { return this->_ptr->write(os); }

  public:
    ScalarFunction derivative(Nat j) const;
    Polynomial<Real> polynomial() const;
  public:
    friend ScalarFunction compose(const ScalarFunction&, const VectorFunction&);
    friend Real evaluate(const ScalarFunction&, const Vector<Real>&);

    friend ScalarFunction operator+(const ScalarFunction&);
    friend ScalarFunction operator-(const ScalarFunction&);
    friend ScalarFunction operator+(const ScalarFunction&, const ScalarFunction&);
    friend ScalarFunction operator-(const ScalarFunction&, const ScalarFunction&);
    friend ScalarFunction operator*(const ScalarFunction&, const ScalarFunction&);
    friend ScalarFunction operator/(const ScalarFunction&, const ScalarFunction&);
  public:
    const ScalarFunctionInterface* _raw_pointer() const { return this->_ptr.operator->(); }
  protected:
    virtual void _check_type(const ScalarFunctionInterface*) { }
  private:
    shared_ptr<ScalarFunctionInterface> _ptr;
};

inline Float evaluate_approx(const ScalarFunction& f, const Vector<Float>& x) { return f(x); }
inline Interval evaluate(const ScalarFunction& f, const Vector<Interval>& x) { return f(x); }
inline TaylorModel compose(const ScalarFunction& f, const Vector<TaylorModel>& x) { return f.evaluate(x); }
inline Vector<Float> gradient_approx(const ScalarFunction& f, const Vector<Float>& x) { return f.gradient(x); }
inline Vector<Interval> gradient(const ScalarFunction& f, const Vector<Interval>& x) { return f.gradient(x); }
inline std::ostream& operator<<(std::ostream& os, const ScalarFunction& f) { return f.write(os); }

ScalarFunction embed(const ScalarFunction&, uint i);

ScalarFunction operator+(const ScalarFunction&);
ScalarFunction operator-(const ScalarFunction&);
ScalarFunction operator+(const ScalarFunction&, const ScalarFunction&);
ScalarFunction operator-(const ScalarFunction&, const ScalarFunction&);
ScalarFunction operator*(const ScalarFunction&, const ScalarFunction&);
ScalarFunction operator/(const ScalarFunction&, const ScalarFunction&);
ScalarFunction operator+(const ScalarFunction&, const Real&);
ScalarFunction operator-(const ScalarFunction&, const Real&);
ScalarFunction operator*(const ScalarFunction&, const Real&);
ScalarFunction operator/(const ScalarFunction&, const Real&);
ScalarFunction operator+(const Real&, const ScalarFunction&);
ScalarFunction operator-(const Real&, const ScalarFunction&);
ScalarFunction operator*(const Real&, const ScalarFunction&);
ScalarFunction operator/(const Real&, const ScalarFunction&);

ScalarFunction pow(const ScalarFunction&, int);
ScalarFunction rec(const ScalarFunction&);
ScalarFunction sqr(const ScalarFunction&);
ScalarFunction sqrt(const ScalarFunction&);
ScalarFunction exp(const ScalarFunction&);
ScalarFunction log(const ScalarFunction&);
ScalarFunction sin(const ScalarFunction&);
ScalarFunction cos(const ScalarFunction&);
ScalarFunction tan(const ScalarFunction&);

//! A vector function \f$f:\R^n\rightarrow\R^m\f$.
class VectorFunction
{
    typedef uint Nat;
    typedef uint SizeType;
    typedef ushort SmoothnessType;
    typedef std::ostream OStream;
  public:
    static VectorFunction constant(const Vector<Real>& c, Nat as);
    static VectorFunction identity(Nat n);

    VectorFunction();
    VectorFunction(Nat rs, Nat as);
    VectorFunction(Nat rs, const ScalarFunction& sf);
    VectorFunction(const List<ScalarFunction>& lsf);
    VectorFunction(const Vector< Polynomial<Real> >& p);
    VectorFunction(VectorFunctionInterface*);
    const VectorFunctionInterface* pointer() const { return this->_ptr.operator->(); }

    VectorFunction(const List< Expression<Real> >& e, const Space<Real>& s);
    VectorFunction(const List< Expression<Real> >& e, const List< Variable<Real> >& v);
    VectorFunction(const List<ExtendedRealVariable>& rs, const Map<ExtendedRealVariable,RealExpression>& e, const List<RealVariable>& as);

    VectorFunction(const List<RealVariable>& rv,
                   const List<RealAssignment>& eq,
                   const List<RealVariable>& av);

    VectorFunction(const List<ExtendedRealVariable>& rv,
                   const List<ExtendedRealAssignment>& eq,
                   const List<RealVariable>& av);

    VectorFunction(const List<DottedRealVariable>& rv,
                   const List<DottedRealAssignment>& eq,
                   const List<RealVariable>& av);

    VectorFunction(const List<PrimedRealVariable>& rv,
                   const List<PrimedRealAssignment>& eq,
                   const List<RealVariable>& av);

    VectorFunction& operator=(const VectorFunction& f) {
        this->_check_type(f._raw_pointer()); this->_ptr=f._ptr; return *this; }

    Nat result_size() const { return this->_ptr->result_size(); }
    Nat argument_size() const { return this->_ptr->argument_size(); }

	/*! \brief Substitute the constant \a con, if present, on all expressions. */
	void substitute(const RealConstant& con);

	/*! \brief Get the parameters (i.e. Constant<Real> whose name starts with a letter) from the function. */
	RealParameterSet parameters() const;

    Vector<Float> evaluate(const Vector<Float>& x) const { return this->_ptr->evaluate(x); }
    Vector<Interval> evaluate(const Vector<Interval>& x) const { return this->_ptr->evaluate(x); }
    Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const { return this->_ptr->evaluate(x); }

    Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const { return this->_ptr->evaluate(x); }
    Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const { return this->_ptr->evaluate(x); }

    Matrix<Float> jacobian(const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).jacobian(); }
    Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).jacobian(); }

    std::ostream& write(std::ostream& os) const { return this->_ptr->write(os); }

    Vector<Float> operator()(const Vector<Float>& x) const { return this->evaluate(x); }
    Vector<Interval> operator()(const Vector<Interval>& x) const { return this->evaluate(x); };
  public:
    ScalarFunction operator[](Nat i) const;

    ScalarFunction get(Nat) const;
    void set(Nat,ScalarFunction);

    VectorFunction polynomial() const;
  public:
    friend VectorFunction join(const ScalarFunction&, const ScalarFunction&);
    friend VectorFunction join(const ScalarFunction&, const VectorFunction&);
    friend VectorFunction join(const VectorFunction&, const ScalarFunction&);
    friend VectorFunction join(const VectorFunction&, const VectorFunction&);

    friend VectorFunction compose(const VectorFunction&, const VectorFunction&);
    friend Vector<Real> evaluate(const VectorFunction&, const Vector<Real>&);

    friend VectorFunction operator+(const VectorFunction&, const VectorFunction&);
    friend VectorFunction operator-(const VectorFunction&, const VectorFunction&);
    friend VectorFunction operator*(const VectorFunction&, const ScalarFunction&);
    friend VectorFunction operator*(const ScalarFunction&, const VectorFunction&);
  public:
    const VectorFunctionInterface* _raw_pointer() const { return this->_ptr.operator->(); }
  protected:
    virtual void _check_type(const VectorFunctionInterface*) { }
  private:
    shared_ptr<VectorFunctionInterface> _ptr;
};

inline Vector<Float> evaluate_approx(const VectorFunction& f, const Vector<Float>& x) { return f(x); }
inline Vector<Interval> evaluate(const VectorFunction& f, const Vector<Interval>& x) { return f(x); }
inline Vector<TaylorModel> compose(const VectorFunction& f, const Vector<TaylorModel>& x) { return f.evaluate(x); }
inline Matrix<Float> jacobian_approx(const VectorFunction& f, const Vector<Float>& x);
inline Matrix<Interval> jacobian(const VectorFunction& f, const Vector<Interval>& x);

VectorFunction operator*(const ScalarFunction& sf, const Vector<Real>& e);
VectorFunction operator+(const VectorFunction& f1, const VectorFunction& f2);
VectorFunction operator-(const VectorFunction& f1, const VectorFunction& f2);
VectorFunction operator*(const VectorFunction& vf, const ScalarFunction& sf);
VectorFunction operator*(const ScalarFunction& sf, const VectorFunction& vf);

VectorFunction join(const ScalarFunction&, const ScalarFunction&);
VectorFunction join(const ScalarFunction&, const VectorFunction&);
VectorFunction join(const VectorFunction&, const ScalarFunction&);
VectorFunction join(const VectorFunction&, const VectorFunction&);

ScalarFunction compose(const ScalarFunction& f, const VectorFunction& g);
VectorFunction compose(const VectorFunction& f, const VectorFunction& g);
ScalarFunction lie_derivative(const ScalarFunction& g, const VectorFunction& f);

inline std::ostream& operator<<(std::ostream& os, const VectorFunction& f) { return f.write(os); }




class ScalarAffineFunction
    : public ScalarFunction
{
  public:
    //! \brief Construct the affine function \f$f(x)=\sum a_ix_i+b\f$.
    ScalarAffineFunction(const Vector<Real>& a, const Real& b);
  protected:
    virtual void _check_type(const ScalarFunctionInterface* ptr) const;
};



class VectorConstantFunction
    : public VectorFunction
{
  public:
    VectorConstantFunction(const Vector<Real>& c, uint as);
  protected:
    virtual void _check_type(const VectorFunctionInterface* ptr) const;
};


class VectorAffineFunction
    : public VectorFunction
{
  public:
    VectorAffineFunction(const Matrix<Real>& A, const Vector<Real>& b);
    const Matrix<Real> A() const;
    const Vector<Real> b() const;
  protected:
    virtual void _check_type(const VectorFunctionInterface* ptr) const;
};



class IdentityFunction
    : public VectorFunction
{
  public:
    //! \brief Construct the identity function in dimension \a n.
    IdentityFunction(uint n);
  protected:
    virtual void _check_type(const VectorFunctionInterface* ptr) const;
};


class ProjectionFunction
    : public VectorFunction
{
  public:
    //! \brief Construct the identity function in dimension \a n.
    ProjectionFunction(uint n);
    //! \brief Construct the projection functions \f$f_i(x)=x_{i+k}\f$ for \f$i=0,\ldots,m-1\f$. Precondition: \f$m+k\leq n\f$.
    ProjectionFunction(uint m, uint n, uint k);
    //! \brief Construct the projection function  with \f$f_i(x)=x_{p_i}\f$ for \f$i=0,\ldots,m-1\f$.
    ProjectionFunction(uint m, uint n, const array<uint>& p);
    //! \brief Construct the projection function with \f$f_i(x)=x_{p_i}\f$ for \f$i=0,\ldots,|p|-1\f$.
    ProjectionFunction(const array<uint>& p, uint n);

    const array<uint>& p() const;
    const uint p(uint i) const;
  protected:
    virtual void _check_type(const VectorFunctionInterface* ptr) const;
};



} // namespace Ariadne

#endif
