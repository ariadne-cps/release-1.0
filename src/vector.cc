/***************************************************************************
 *            vector.cc
 *
 *  Copyright 2008-11  Alberto Casagrande, Pieter Collins
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

#include "config.h"
#include "macros.h"
#include "numeric.h"
#include "real.h"
#include "vector.h"

template class boost::numeric::ublas::vector<Ariadne::Float>;
template class boost::numeric::ublas::vector<Ariadne::Interval>;

namespace Ariadne {


template<> Vector<Float>::Vector(size_t n, const double& t0, const double& t1, ...)
    : ublas::vector<Float>(n)
{
    assert(n>=2); va_list args; va_start(args,t1);
    (*this)[0]=t0; (*this)[1]=t1;
    for(size_t i=2; i!=n; ++i) { (*this)[i]=va_arg(args,double); }
    va_end(args);
}

template<> Vector<Real>::Vector(size_t n, const double& t0, const double& t1, ...)
    : ublas::vector<Real>(n)
{
    assert(n>=2); va_list args; va_start(args,t1);
    (*this)[0]=t0; (*this)[1]=t1;
    for(size_t i=2; i!=n; ++i) { (*this)[i]=va_arg(args,double); }
    va_end(args);
}

template<> Vector<Interval>::Vector(size_t n, const double& t0, const double& t1, ...)
    : ublas::vector<Interval>(n)
{
    assert(n>=1); va_list args; va_start(args,t1);
    (*this)[0]=Interval(t0,t1);
    for(size_t i=1; i!=n; ++i) {
        double l=va_arg(args,double);
        double u=va_arg(args,double);
        (*this)[i]=Interval(l,u);
    }
    va_end(args);
}

#ifdef HAVE_RATIONAL
template<> Vector<Rational>::Vector(size_t n, const double& t0, const double& t1, ...)
    : ublas::vector<Rational>(n)
{
    assert(n>=2); va_list args; va_start(args,t1);
    (*this)[0]=t0; (*this)[1]=t1;
    for(size_t i=2; i!=n; ++i) { (*this)[i]=va_arg(args,double); }
    va_end(args);
}
#endif // HAVE_RATIONAL

Vector<Interval> operator*(const Vector<Interval>& v, const Float& s)
{
  return v*Interval(s);
}

Vector<Interval> operator*(const Float& s, const Vector<Interval>& v)
{
  return v*s;
}

Vector<Interval> operator*(const Vector<Interval>& v, const int& s)
{
  return v*Interval(s);
}

Vector<Interval> operator*(const int& s, const Vector<Interval>& v)
{
  return v*s;
}

Vector<Interval> operator/(const Vector<Interval>& v, const Float& s)
{
  return v/Interval(s);
}

Vector<Interval> operator/(const Vector<Interval>& v, const int& s)
{
  return v/Interval(s);
}

Vector<Float> operator*(const Vector<Float>& v, const int& s)
{
  return v*Float(s);
}

Vector<Float> operator*(const int& s, const Vector<Float>& v)
{
  return v*s;
}

Vector<Float> operator/(const Vector<Float>& v, const int& s)
{
  return v/Float(s);
}


bool contains(const Vector<Interval>& v1, const Vector<Float>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!contains(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool subset(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!subset(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool intersect(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!intersect(v1[i],v2[i])) { return false; }
    }
    return true;
}


bool disjoint(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(disjoint(v1[i],v2[i])) { return true; }
    }
    return false;
}

bool overlap(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!overlap(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool covers(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!covers(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool inside(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!inside(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool empty(const Vector<Interval>& v)
{
    for(size_t i=0; i!=v.size(); ++i) {
        if(empty(v[i])) { return true; }
    }
    return false;
}


uint irmax(const Vector<Interval>& v) {
    uint imw(0);
    Float mw=v[0].width();
    for(uint i=1; i!=v.size(); ++i) {
        if(v[i].width()>mw) { imw=i; mw=v[i].width(); }
    }
    return imw;
}


Vector<Interval> split(const Vector<Interval>& v, uint k, tribool lr) {
    ARIADNE_ASSERT(k<v.size());
    Vector<Interval> r(v);
    Float c=v[k].midpoint();
    if(lr) {
        r[k].set_upper(c);
    } else if(!lr) {
        r[k].set_lower(c);
    } else {
        Float cl=(3*v[k].lower()+v[k].upper())/4;
        Float cu=(v[k].lower()+3*v[k].upper())/4;
        r[k].set_lower(cl);
        r[k].set_upper(cu);
    }
    return r;
}

std::pair< Vector<Interval>, Vector<Interval> > split(const Vector<Interval>& v, uint k) {
    ARIADNE_ASSERT(k<v.size());
    std::pair< Vector<Interval>, Vector<Interval> > r(v,v);
    Float c=v[k].midpoint();
    r.first[k].set_upper(c);
    r.second[k].set_lower(c);
    return r;
}

Vector<Interval> split(const Vector<Interval>& v, tribool lr) {
    return split(v,irmax(v),lr);
}

std::pair< Vector<Interval>, Vector<Interval> > split(const Vector<Interval>& v) {
    return split(v,irmax(v));
}



Vector<Float> midpoint(const Vector<Interval>& v)
{
    Vector<Float> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].midpoint();
    }
    return r;
}

Vector<Float> lower(const Vector<Interval>& v)
{
    Vector<Float> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].lower();
    }
    return r;
}

Vector<Float> upper(const Vector<Interval>& v)
{
    Vector<Float> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].upper();
    }
    return r;
}

Vector<Interval> hull(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    Vector<Interval> r(v1.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        r[i]=hull(v1[i],v2[i]);
    }
    return r;
}

Vector<Interval> intersection(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    Vector<Interval> r(v1.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        r[i]=intersection(v1[i],v2[i]);
    }
    return r;
}

Float radius(const Vector<Interval>& v)
{
    Float r=0;
    for(size_t i=0; i!=v.size(); ++i) {
        r=Ariadne::max(r,v[i].radius());
    }
    return r;
}

Float volume(const Vector<Interval>& v)
{
    Float r=1.0;
    for(size_t i=0; i!=v.size(); ++i) {
        r*=diam(v[i]);
    }
    return r;
}

Vector<Interval> intersection(const Vector<Interval>& v1, const Vector<Interval>& v2);
Vector<Float> midpoint(const Vector<Interval>& v);


} // namespace Ariadne
