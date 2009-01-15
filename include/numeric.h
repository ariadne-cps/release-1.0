/***************************************************************************
 *            numeric.h
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
 
/*! \file numeric.h
 *  \brief Numerical classes.
 */
#ifndef ARIADNE_NUMERIC_H
#define ARIADNE_NUMERIC_H

#ifdef HAVE_GMPXX_H
#include <gmpxx.h>
#endif // HAVE_GMPXX_H

#include <cmath>
#include <limits>
#include <stdint.h>

#include "tribool.h"
#include "rounding.h"

// Simplifying typedefs for unsigned types
// These may be inclused in other headers,
// but repeating a typedef is not an error
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

namespace Ariadne {

#ifdef DOXYGEN
//! \brief Integers of arbitrary size. 
//! (Only available if the Gnu Multiple Precision library (GMP) is installed.)
class Integer { };
//! \brief Rationals numbers.
//! (Only available if the Gnu Multiple Precision library (GMP) is installed.)
class Rational { };
//! \brief Floating point numbers (double precision).
class Float { };
#endif // DOXYGEN


#ifdef HAVE_GMPXX_H
typedef mpq_class Rational;
Rational sqr(const Rational& q);
Rational pow(const Rational& q, int n);
Rational pow(const Rational& q, uint n);
#else
#endif // HAVE_GMPXX_H 


typedef double Float;


using std::min;
using std::max;
using std::abs;

uint16_t fac(uint16_t n);
uint32_t fac(uint32_t n);
uint64_t fac(uint64_t n);
uint16_t bin(uint16_t n, uint16_t k);
uint32_t bin(uint32_t n, uint32_t k);
uint64_t bin(uint64_t n, uint64_t k);


template<class X> X pi();

inline Float inf() { return std::numeric_limits<double>::max(); }
inline Float eps() { return std::numeric_limits<double>::epsilon(); }

inline Float down(Float x) { return x>0 ? x*(1-2e-16) : x*(1+2e-16); }
inline Float up(Float x) { return x>0 ? x*(1+2e-16) : x*(1-2e-16); }

inline Float neg(Float x) { return -x; }
inline Float rec(Float x) { return 1.0/x; }

inline Float add(Float x, Float y) { return x+y; }
inline Float sub(Float x, Float y) { return x-y; }
inline Float mul(Float x, Float y) { return x*y; }
inline Float div(Float x, Float y) { return x/y; }

inline Float sqr(Float x) { return x*x; }
inline Float pow(Float x, int n) { return std::pow(x,Float(n)); }
inline Float pow(Float x, uint n) { return std::pow(x,Float(n)); }

inline Float sqrt(Float x) { return std::sqrt(x); }
inline Float exp(Float x) { return std::exp(x); }
inline Float log(Float x) { return std::log(x); }

template<> inline Float pi<Float>() { return 3.1415926535897931; }
inline Float sin(Float x) { return std::sin(x); }
inline Float cos(Float x) { return std::cos(x); }
inline Float tan(Float x) { return std::tan(x); }
inline Float asin(Float x) { return std::asin(x); }
inline Float acos(Float x) { return std::acos(x); }
inline Float atan(Float x) { return std::atan(x); }

inline Float add_up(Float x, Float y) { return up(x+y); }
inline Float sub_up(Float x, Float y) { return up(x-y); }
inline Float mul_up(Float x, Float y) { return up(x*y); }
inline Float div_up(Float x, Float y) { return up(x/y); }
inline Float pow_up(Float x, int n) { return up(pow(x,n)); }

inline Float add_down(Float x, Float y) { return down(x+y); }
inline Float sub_down(Float x, Float y) { return down(x-y); }
inline Float mul_down(Float x, Float y) { return down(x*y); }
inline Float div_down(Float x, Float y) { return down(x/y); }
inline Float pow_down(Float x, int n) { return down(pow(x,n)); }

inline Float add_approx(Float x, Float y) { return x+y; }
inline Float sub_approx(Float x, Float y) { return x-y; }
inline Float mul_approx(Float x, Float y) { return x*y; }
inline Float div_approx(Float x, Float y) { return x/y; }
inline Float pow_approx(Float x, int n) { return pow(x,n); }

inline Float rad_up(Float x, Float y) { return up((y-x)/2); }
inline Float med_approx(Float x, Float y) { return (x+y)/2; }

inline Float add_rnd(Float x, Float y) { return x+y; }
inline Float sub_rnd(Float x, Float y) { return x-y; }
inline Float mul_rnd(Float x, Float y) { return x*y; }
inline Float div_rnd(Float x, Float y) { return x/y; }

inline Float add_opp(Float x, Float y) { volatile double t=(-x)-y; return -t; }
inline Float sub_opp(Float x, Float y) { volatile double t=(-x)+y; return -t; }
inline Float mul_opp(Float x, Float y) { volatile double t=(-x)*y; return -t; }
inline Float div_opp(Float x, Float y) { volatile double t=(-x)/y; return -t; }


//! \brief Intervals supporting interval arithmetic.
class Interval {
  public:
    Interval() : l(0.0), u(0.0) { }
    Interval(uint m) : l(m), u(m) { }
    Interval(int n) : l(n), u(n) { }
    Interval(Float x) : l(x), u(x) { }
    Interval(const Interval& i) : l(i.l), u(i.u) { }
  
    Interval(Float lower, Float upper) : l(lower), u(upper) { }
#ifdef HAVE_GMPXX_H
    Interval(Rational q);
    Interval(Rational lower, Rational upper);
#endif // HAVE_GMPXX_H 

    const Float& lower() const { return l; }
    const Float& upper() const { return u; }
    const Float midpoint() const { return add_approx(l,u)/2; }
    const Float radius() const { return sub_up(u,l)/2; }
    const Float width() const { return sub_up(u,l); }

    bool empty() const { return l>u; }
    bool singleton() const { return l==u; }

    void set_lower(const Float& lower) { l=lower; }
    void set_upper(const Float& upper) { u=upper; }
  public:
    double l, u;
};

std::ostream& operator<<(std::ostream& os, const Interval& ivl);

inline Float midpoint(Interval i) { 
    return add_approx(i.l,i.u)/2; 
}

inline Float radius(Interval i) { 
    return sub_up(i.u,i.l)/2; 
}

inline Float width(Interval i) { 
    return sub_up(i.u,i.l); 
}

inline bool equal(Interval i1, Interval i2) { 
    //std::cerr<<"equal(i1,i2) with i1="<<i1<<"; i2="<<i2<<std::endl;
    return i1.l==i2.l && i1.u==i2.u;
}

inline bool empty(Interval i) { 
    return i.l>i.u;
}

inline bool bounded(Interval i) { 
    return i.l!=-inf() && i.u!=+inf();
}

inline Interval intersection(Interval i1, Interval i2) { 
    if(i1.l>i2.u || i1.u<i2.l) { return Interval(1,-1); }
    return Interval(max(i1.l,i2.l),min(i1.u,i2.u));
}

inline Interval hull(Interval i1, Interval i2) { 
    return Interval(min(i1.l,i2.l),max(i1.u,i2.u));
}

Interval trunc(Interval, uint eps);

inline Float med(Interval i) { return (i.l+i.u)/2; }
inline Float rad(Interval i) { return up((i.u-i.l)/2); }
inline Float diam(Interval i) { return up(i.u-i.l); }

inline Interval abs(Interval);
inline Interval neg(Interval);
inline Interval add(Interval, Interval);
inline Interval add(Interval, Float);
inline Interval sub(Interval, Interval);
inline Interval sub(Interval, Float);
inline Interval sub(Float, Interval);

Interval rec(Interval);
Interval mul(Interval, Interval);
Interval div(Interval, Interval);
Interval mul(Interval, Float);
Interval div(Interval, Float);
Interval div(Float, Interval);

Interval sqr(Interval);
Interval pow(Interval, uint);
Interval pow(Interval, int);

Interval sqrt(Interval);
Interval exp(Interval);
Interval log(Interval);

template<> Interval pi<Interval>();
Interval sin(Interval);
Interval cos(Interval);
Interval tan(Interval);
Interval asin(Interval);
Interval acos(Interval);
Interval atan(Interval);


inline Float mag(Interval i) { return max(abs(i.l),abs(i.u)); }
inline Float mig(Interval i) { return min(abs(i.l),abs(i.u)); }

inline bool contains(Interval i, Float x) { return i.l<=x && x<=i.u; }

inline bool subset(Interval i1, Interval i2) { return i1.l>=i2.l && i1.u<=i2.u; }
inline bool intersect(Interval i1, Interval i2) { return i1.l<=i2.u && i1.u>=i2.l; }
inline bool disjoint(Interval i1, Interval i2) { return i1.l>i2.u || i1.u<i2.l; }
inline bool overlap(Interval i1, Interval i2) { return i1.l<i2.u && i1.u>i2.l; }
inline bool inside(Interval i1, Interval i2) { return i1.l>i2.l && i1.u<i2.u; }
inline bool covers(Interval i1, Interval i2) { return i1.l<i2.l && i1.u>i2.u; }

Interval abs(Interval i) 
{
    if(i.l>=0) {
        return Interval(i.l,i.u);
    } else if(i.u<=0) {
        return Interval(-i.u,-i.l);
    } else {
        return Interval(0.0,max(-i.l,i.u));
    }
}

Interval neg(Interval i) 
{
    return Interval(-i.u,-i.l);
}

Interval add(Interval i1, Interval i2) 
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=i1.l;
    volatile double& i1u=i1.u;
    volatile double& i2l=i2.l;
    volatile double& i2u=i2.u;
    set_rounding_mode(downward);
    volatile double rl=i1l+i2l;
    set_rounding_mode(upward);
    volatile double ru=i1u+i2u;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

Interval add(Interval i1, Float x2) 
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=i1.l;
    volatile double& i1u=i1.u;
    volatile double& x2v=x2;
    set_rounding_mode(downward);
    volatile double rl=i1l+x2v;
    set_rounding_mode(upward);
    volatile double ru=i1u+x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

Interval sub(Interval i1, Interval i2) 
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=i1.l;
    volatile double& i1u=i1.u;
    volatile double& i2l=i2.l;
    volatile double& i2u=i2.u;
    set_rounding_mode(downward);
    volatile double rl=i1l-i2u;
    set_rounding_mode(upward);
    volatile double ru=i1u-i2l;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

Interval sub(Interval i1, Float x2) 
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& i1l=i1.l;
    volatile double& i1u=i1.u;
    volatile double& x2v=x2;
    set_rounding_mode(downward);
    volatile double rl=i1l-x2v;
    set_rounding_mode(upward);
    volatile double ru=i1u-x2v;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

Interval sub(Float x1, Interval i2) 
{
    rounding_mode_t rnd=get_rounding_mode();
    volatile double& x1v=x1;
    volatile double& i2l=i2.l;
    volatile double& i2u=i2.u;
    set_rounding_mode(downward);
    volatile double rl=x1v-i2u;
    set_rounding_mode(upward);
    volatile double ru=x1v-i2l;
    set_rounding_mode(rnd);
    return Interval(rl,ru);
}

// Standard equality operators
inline bool operator==(const Interval& i1, const Interval& i2) { return i1.l==i2.l && i1.u==i2.u; }
inline bool operator!=(const Interval& i1, const Interval& i2) { return i1.l!=i2.l || i1.u!=i2.u; }

// Boost-style tribool (in)equality operators
//inline tribool operator==(const Interval& i1, const Interval& i2) { 
//  if(i1.l>i2.u || i1.u<i2.l) { return false; } else if(i1.l==i2.u && i1.u==i2.l) { return true; } else { return indeterminate; } }
//inline tribool operator!=(const Interval& i1, const Interval& i2) { return !(i1==i2); }


inline Interval operator+(Interval i) { return Interval(i.l,i.u); }
inline Interval operator-(Interval i) { return Interval(-i.u,-i.l); }
inline Interval operator+(Interval i1, Interval i2) { return add(i1,i2); }
inline Interval operator-(Interval i1, Interval i2) { return sub(i1,i2); }
inline Interval operator*(Interval i1, Interval i2) { return mul(i1,i2); }
inline Interval operator/(Interval i1, Interval i2) { return div(i1,i2); };

inline Interval& operator+=(Interval& i1, Interval i2) { i1=add(i1,i2); return i1; }
inline Interval& operator-=(Interval& i1, Interval i2) { i1=sub(i1,i2); return i1; }
inline Interval& operator*=(Interval& i1, Interval i2) { i1=mul(i1,i2); return i1; }
inline Interval& operator/=(Interval& i1, Interval i2) { i1=div(i1,i2); return i1; }

inline Interval operator+(Interval i1, Float x2) { return add(i1,x2); }
inline Interval operator+(Float x1, Interval i2) { return add(i2,x1); }
inline Interval operator-(Interval i1, Float x2) { return sub(i1,x2); }
inline Interval operator-(Float x1, Interval i2) { return sub(x1,i2); }
inline Interval operator*(Interval i1, Float x2) { return mul(i1,x2); }
inline Interval operator*(Float x1, Interval i2) { return mul(i2,x1); }
inline Interval operator/(Interval i1, Float x2) { return div(i1,x2); }
inline Interval operator/(Float x1, Interval i2) { return div(x1,i2); }

inline Interval& operator+=(Interval& i1, Float x2) { i1=add(i1,x2); return i1; }
inline Interval& operator-=(Interval& i1, Float x2) { i1=sub(i1,x2); return i1; }
inline Interval& operator*=(Interval& i1, Float x2) { i1=mul(i1,x2); return i1; }
inline Interval& operator/=(Interval& i1, Float x2) { i1=div(i1,x2); return i1; }

inline tribool operator==(Interval i1, Float x2) { 
    if(i1.upper()<x2 || i1.lower()>x2) { return false; }
    else if(i1.lower()==x2 && i1.upper()==x2) { return true; }
    else { return indeterminate; }
}

inline tribool operator!=(Interval i1, Float x2) { 
    if(i1.upper()<x2 || i1.lower()>x2) { return true; }
    else if(i1.lower()==x2 && i1.upper()==x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator> (Interval i1, Float x2) { 
    if(i1.lower()> x2) { return true; }
    else if(i1.upper()<=x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator< (Interval i1, Float x2) { 
    if(i1.upper()< x2) { return true; }
    else if(i1.lower()>=x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator>=(Interval i1, Float x2) { 
    if(i1.lower()>=x2) { return true; }
    else if(i1.upper()< x2) { return false; }
    else { return indeterminate; }
}

inline tribool operator<=(Interval i1, Float x2) { 
    if(i1.upper()<=x2) { return true; }
    else if(i1.lower()> x2) { return false; }
    else { return indeterminate; }
}



inline tribool operator> (Interval i1, Interval i2) { 
    if(i1.lower()> i2.upper()) { return true; }
    else if(i1.upper()<=i2.lower()) { return false; }
    else { return indeterminate; }
}

inline tribool operator< (Interval i1, Interval i2) { 
    if(i1.upper()< i2.lower()) { return true; }
    else if(i1.lower()>=i2.upper()) { return false; }
    else { return indeterminate; }
}

inline tribool operator>=(Interval i1, Interval i2) { 
    if(i1.lower()>=i2.upper()) { return true; }
    else if(i1.upper()< i2.lower()) { return false; }
    else { return indeterminate; }
}

inline tribool operator<=(Interval i1, Interval i2) { 
    if(i1.upper()<=i2.lower()) { return true; }
    else if(i1.lower()> i2.upper()) { return false; }
    else { return indeterminate; }
}

template<class A> void serialize(A& a, Interval& ivl, const uint version) {
    a & ivl.l & ivl.u; }

std::ostream& operator<<(std::ostream&, const Interval&);
std::istream& operator>>(std::istream&, Interval&);

} // namespace Ariadne 

#endif