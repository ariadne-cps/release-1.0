/***************************************************************************
 *            grid.cc
 *
 *  Copyright  2010  Pieter Collins, Luca Geretti
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

#include <iostream>
#include <fstream>
#include <iomanip>
#include "array.h"
#include "numeric.h"
#include "vector.h"
#include "grid.h"
#include "box.h"
#include "macros.h"
#include "exceptions.h"
#include "stlio.h"

namespace Ariadne {

typedef size_t size_type;

/****************************************Grid**********************************************/

Grid::~Grid()
{
}
 
Grid::Grid()
    : _data()
{
}
 
Grid::Grid(const Grid& gr)
    : _data(gr._data)
{
}
 
Grid::Grid(uint d)
    : _data()
{
    Vector<Float> origin(d,Float(0));
    Vector<Float> lengths(d,Float(1));
    this->_create(origin,lengths);
}
 
Grid::Grid(uint d, Float l)
    : _data()
{
    Vector<Float> origin(d,Float(0));
    Vector<Float> lengths(d,l);
    this->_create(origin,lengths);
}

Grid::Grid(const Vector<Float>& lengths)
    : _data()
{
    Vector<Float> origin(lengths.size(),0);
    this->_create(origin,lengths);
}
 
Grid::Grid(const Vector<Float>& origin, const Vector<Float>& lengths)
    : _data()
{
    if(origin.size() != lengths.size()) {
        throw IncompatibleSizes(ARIADNE_PRETTY_FUNCTION);
    }
    this->_create(origin,lengths);
}

void Grid::_create(const Vector<Float>& origin, const Vector<Float>& lengths) 
{
    this->_data._origin=origin;
    this->_data._lengths=lengths;
}

uint Grid::dimension() const
{
    return this->_data._lengths.size();
}

const Vector<Float>& Grid::origin() const
{
    return this->_data._origin;
}

const Vector<Float>& Grid::lengths() const
{
    return this->_data._lengths;
}

void Grid::set_origin(const Vector<Float>& origin)
{
    if(origin.size() != this->_data._origin.size()) {
        throw IncompatibleSizes(ARIADNE_PRETTY_FUNCTION);
    }
    this->_data._origin=origin;
}

void Grid::set_origin_coordinate(uint i, Float o)
{
    if(i >= this->_data._origin.size()) {
        throw IncompatibleSizes(ARIADNE_PRETTY_FUNCTION);
    }
    this->_data._origin[i]=o;
}


void Grid::set_lengths(const Vector<Float>& lengths)
{
    if(lengths.size() != this->_data._lengths.size()) {
        throw IncompatibleSizes(ARIADNE_PRETTY_FUNCTION);
    }
    this->_data._lengths=lengths;
}

void Grid::set_length(uint i, Float l)
{
    if(i >= this->_data._lengths.size()) {
        throw IncompatibleSizes(ARIADNE_PRETTY_FUNCTION);
    }
    this->_data._lengths[i]=l;
}


Float Grid::coordinate(uint d, dyadic_type x) const 
{
    return add_approx(this->_data._origin[d],mul_approx(this->_data._lengths[d],x));
}

Float Grid::subdivision_coordinate(uint d, dyadic_type x) const 
{
    return add_approx(this->_data._origin[d],mul_approx(this->_data._lengths[d],x));
}

Float Grid::subdivision_coordinate(uint d, integer_type n) const 
{
    return add_approx(this->_data._origin[d],mul_approx(this->_data._lengths[d],n));
}

int Grid::subdivision_index(uint d, const real_type& x) const 
{
    Float half=0.5;
    int n=int(floor(add_approx(div_approx(sub_approx(x,this->_data._origin[d]),this->_data._lengths[d]),half)));
    Float sc=add_approx(this->_data._origin[d],mul_approx(this->_data._lengths[d],n));
    if(sc == x) { 
        return n; 
    } else {
        std::cerr << std::setprecision(20) << std::boolalpha
                  << "sc=" << sc << " x=" << x << " sc-x=" << Interval(sc-x) << "\n"
                  << "sc==x=" << (sc==x) << " sc!=x=" << (sc!=x)
                  << " sc<x=" << (sc<x) << " sc>x=" << (sc>x) << " sc<=x=" << (sc<=x) << " sc>=x=" << (sc>=x) << std::endl; 
        ARIADNE_THROW(InvalidGridPosition,std::setprecision(20)<<"Grid::subdivision_index(uint d,real_type x)","d="<<d<<", x="<<x<<", this->origin[d]="<<this->_data._origin[d]<<", this->lengths[d]="<<this->_data._lengths[d]<<" (closest value is "<<sc<<")");
    }
}
 
int Grid::subdivision_lower_index(uint d, const real_type& x) const 
{
    int n=int(floor(div_down(sub_down(x,this->_data._origin[d]),this->_data._lengths[d])));
    if(x>=add_approx(this->_data._origin[d],mul_approx(this->_data._lengths[d],(n+1)))) {
        return n+1;
    } else {
        return n;
    }
}
 
int Grid::subdivision_upper_index(uint d, const real_type& x) const 
{
    int n=int(ceil(div_up(sub_up(x,this->_data._origin[d]),this->_data._lengths[d])));
    if(x<=add_approx(this->_data._origin[d],mul_approx(this->_data._lengths[d],(n-1)))) {
        return n-1;
    } else {
        return n;
    }
}
 
bool Grid::operator==(const Grid& g) const
{
    return this->_data._origin==g._data._origin && this->_data._lengths==g._data._lengths;
}
 
bool Grid::operator!=(const Grid& g) const
{
    return !(*this==g);
}

array<double> Grid::index(const Vector<Float>& pt) const
{
    array<double> res(pt.size());
    for(size_t i=0; i!=res.size(); ++i) {
        res[i]=subdivision_index(i,pt[i]);
    }
    return res;
}

array<double> Grid::lower_index(const Vector<Interval>& bx) const {
    array<double> res(bx.size());
    for(size_t i=0; i!=res.size(); ++i) {
        res[i]=subdivision_lower_index(i,bx[i].lower());
    }
    return res;
}

array<double> Grid::upper_index(const Vector<Interval>& bx) const {
    array<double> res(bx.size());
    for(size_type i=0; i!=res.size(); ++i) {
        res[i]=subdivision_upper_index(i,bx[i].upper());
    }
    return res;
}

Vector<Float> Grid::point(const array<int>& a) const
{
    Vector<Float> res(a.size());
    for(size_type i=0; i!=res.size(); ++i) {
        res[i]=this->_data._origin[i]+this->_data._lengths[i]*a[i];
    }
    return res;
}

Vector<Float> Grid::point(const array<double>& a) const
{
    ARIADNE_ASSERT(a.size() == this->dimension());
    Vector<float> res(a.size());
    for(size_type i=0; i!=res.size(); ++i) {
        res[i]=this->_data._origin[i]+this->_data._lengths[i]*a[i];
    }
    return res;
}

Vector<Interval> Grid::cell(const array<int>& a) const
{
    ARIADNE_ASSERT(a.size() == this->dimension());
    Vector<Interval> res(a.size());
    for(size_type i=0; i!=res.size(); ++i) {
        res[i]=Interval(this->_data._origin[i]+this->_data._lengths[i]*a[i],
                        this->_data._origin[i]+this->_data._lengths[i]*(a[i]+1));
    }
    return res;
}

Vector<Interval> Grid::box(const array<double>& lower, const array<double>& upper) const
{
    Vector<Interval> res(lower.size());
    for(size_type i=0; i!=res.size(); ++i) {
        res[i]=Interval(this->subdivision_coordinate(i,lower[i]),
                        this->subdivision_coordinate(i,upper[i]));
    }
    return res;
}

Box Grid::primary_cell() const
{
    Box res(this->dimension());
    for(size_type i=0; i!=res.size(); ++i) {
        res[i]=Interval(_data._origin[i],
                        _data._origin[i]+this->_data._lengths[i]);
    }
    return res;
}

Grid project_down(const Grid& original_grid, const Vector<uint>& indices)
{
	uint original_space_size = original_grid.dimension();
	for (Vector<uint>::const_iterator dim_it = indices.begin(); dim_it != indices.end(); ++dim_it) {
		ARIADNE_ASSERT_MSG(*dim_it < original_space_size,"Index " << *dim_it << " is outside the original space size (" << original_space_size << ")");
	}

	Vector<double> original_grid_lengths = original_grid.lengths();
	Vector<double> projected_grid_lengths(indices.size());
	for (uint i=0; i < indices.size(); ++i) {
		projected_grid_lengths[i] = original_grid_lengths[indices[i]];
	}

	return Grid(projected_grid_lengths);
}


std::ostream& operator<<(std::ostream& os, const Grid& gr) 
{
    os << "Grid( ";
    os << "origin=" << gr.origin() << ", ";
    os << "lengths=" << gr.lengths() << " )";
    return os;
}

} // namespace Ariadne

