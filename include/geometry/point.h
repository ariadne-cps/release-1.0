/***************************************************************************
 *            point.h
 *
 *  Sun Jan 23 18:00:21 2005
 *  Copyright  2005  Alberto Casagrande
 *  Email casagrande@dimi.uniud.it
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

/*! \file point.h
 *  \brief A point in Euclidean space.
 */

#ifndef _ARIADNE_POINT_H
#define _ARIADNE_POINT_H

#include <iosfwd>
#include <stdexcept>

#include "../declarations.h"
#include "../base/array.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include <gmpxx.h>

namespace Ariadne { 
  namespace Numeric {
    typedef mpq_class Rational;
  } 
}

namespace Ariadne {
  namespace Geometry {
   
    template<typename> class Rectangle;
    template<typename> class Parallelotope;
    template<typename> class Zonotope;
    template<typename> class Sphere;
    template<typename> class Ellipsoid;

    template<template <typename> class BS1, template <typename> class BS2> 
    inline bool is_a(){ return false; }

    template<> inline bool is_a<Point, Point>(){ return true; }

    /* Forward declaration of friends. */
    template<typename R> std::ostream& operator<<(std::ostream&, const Point<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Point<R>&);

    /*! \brief A point in Euclidean space. */
    template <typename R>
    class Point {
     public:
      /*!\brief The type of denotable real number giving the point's values. */
      typedef R real_type;
      /*!\brief The type of denotable vector in the space the point lies in. */
      typedef LinearAlgebra::Vector<R> vector_type;

      /*! \brief Default constructor. */
      Point() : _vector(0) { }

      /*! \brief The origin in dimension \a dim. */
      explicit Point(size_type dim) : _vector(dim) {
        for(size_type i=0; i!=dimension(); ++i) {
          _vector[i]=real_type(0);
        }
      }

      /*! \brief Construct a point from a range of values. */
      template<class ForwardIterator>
      Point(ForwardIterator b, ForwardIterator e) : _vector(std::distance(b,e))
      {
        for(size_type i=0; i!=dimension(); ++i) {
          _vector[i]=*b;
          ++b;
        }
      }

      /*! \brief Construct a point from a position Vector. */
      explicit Point(const vector_type& position) : _vector(position) { }

      /*! \brief Construct a point from a string literal. */
      explicit Point(const std::string& s);

      /*! \brief Construct from a point which possibly lies in a different space. */
      template<typename R2>
      explicit Point(const Point<R2>& original) : _vector(original.position_vector()) { }
      
      /*! \brief Copy constructor. */
      Point(const Point<R>& original) : _vector(original._vector) { }

      /*! \brief Assignment operator. */
      Point<R>& operator=(const Point<R>& original) {
        if(this!=&original) { this->_vector=original._vector; }
        return *this; 
      }
      
      #ifndef RATIONAL_REAL
      /*! \brief Convert to a rational point. */
      operator Point<Rational> () const;
      #endif 


      /*! \brief Checks equivalence between two states. */
      bool operator==(const Point<R>& A) const {
        /* Return false if states have different dimensions */
        if (this->dimension()!=A.dimension()) { return false; }
        for (size_type i=0; i<this->dimension(); i++) {
          if (this->_vector[i]!=A._vector[i]) { return false; }
        }
        return true;
      }
      
      /*! \brief Inequality operator. */
      bool operator!=(const Point<real_type>& A) const {
        return !( *this == A );
      }

      /*! \brief The dimension of the Euclidean space the state lies in. */
      size_type dimension() const {
        return this->_vector.size();
      }

      /*! \brief Subcripting operator. */
      real_type& operator[] (size_type index) {
        if(this->_vector.size() <= index) { 
          throw std::out_of_range("Out of the Vector's range.");
        }
        return  (this->_vector[index]);
      }

      /*! \brief Subcripting operator. */
      const real_type& operator[](size_type index) const {
        if(this->_vector.size() <= index) { 
            throw std::out_of_range("Out of the Vector's range.");
        }
        return  (this->_vector[index]);
      }

      /*! \brief The position Vector of the point. */
      const vector_type& position_vector() const {
        return this->_vector; 
      }
      
/*      real_type get(size_type index) const {
        if(this->_vector.size() <= index) { 
          throw std::out_of_range("Out of the Vector's range.");
        }
        return  (this->_vector[index]);
      }

      void set(size_type index, const real_type& r) {
        if(this->_vector.size() <= index) {
          throw std::out_of_range("Out of the Vector's range.");
        }
        this->_vector[index]=r;
      }
*/
      
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);

      friend std::ostream& operator<< <>(std::ostream& os, const Point<real_type>& state);
      friend std::istream& operator>> <> (std::istream& is, Point<real_type>& state);
     private:
      vector_type _vector;
    };


    template <typename R> 
    inline
    LinearAlgebra::Vector<R>
    operator-(const Point<R>& s1, const Point<R>& s2)
    {
      return s1.position_vector() - s2.position_vector();
    }

    template <typename R>
    inline
    Point<R> 
    operator+(const Point<R>& s, const LinearAlgebra::Vector<R>& v)
    {
      return Point<R>(s.position_vector() + v);
    }

    template <typename R>
    inline
    Point<R> 
    operator-(const Point<R>& s, const LinearAlgebra::Vector<R>& v)
    {
      return Point<R>(s.position_vector() - v);
    }

    template <typename R>
    inline
    LinearAlgebra::Vector<R>
    sub_approx(const Point<R>& s1, const Point<R>& s2)
    {
      return sub_approx(s1.position_vector(),s2.position_vector());
    }

    template <typename R>
    inline
    Point<R> 
    add_approx(const Point<R>& s, const LinearAlgebra::Vector<R>& v)
    {
      return Point<R>(add_approx(s.position_vector(),v));
    }

    template <typename R>
    inline
    Point<R> 
    sub_approx(const Point<R>& s, const LinearAlgebra::Vector<R>& v)
    {
      return Point<R>(sub_approx(s.position_vector(),v));
    }

    template <typename R>
    inline 
    Point<R> 
    project_on_dimensions(const Point<R> &A, const Base::array<bool>& dims) 
    {
      if (A.dimension()!=dims.size()) {
         throw "project_on_dimensions(const Point & ,...): the two parameters have different dimension";
      }
      size_type new_dim=0;

      for (size_type i=0; i< A.dimension(); i++) {
        if (dims[i]) {
          new_dim++;
        }
      }
      Point<R> new_point(new_dim);
      
      size_type new_i=0;
      for (size_t i=0; i<dims.size(); i++) {
        if (dims[i]) {
           new_point.set(new_i,A[i]);
           new_i++;
        } 
      }

      return new_point;
    }

    template <typename R>
    inline 
    Point<R> 
    project_on_dimensions(const Point<R> &A, 
                          const size_type& x, const size_type& y, const size_type& z) 
    {

      if ((A.dimension()<=x)||(A.dimension()<=y)||(A.dimension()<=z)) {
         throw "project_on_dimensions(const Point& ,...): one of the projection dimensions is greater than the Point dimension";
      }
      
      Point<R> new_point(3);
      
      new_point.set(0,A[x]);
      new_point.set(1,A[y]);
      new_point.set(2,A[z]);

      return new_point;
    }

    template <typename R>
    inline 
    Point<R> 
    project_on_dimensions(const Point<R> &A, 
                          const size_type &x, const size_type&y) 
    {
      if ((A.dimension()<=x)||(A.dimension()<=y)) {
         throw "project_on_dimensions(const Point& ,...): one of the projection dimensions is greater than the Point dimension";
      }
      Point<R> new_point(2);
      
      new_point.set(0,A[x]);
      new_point.set(1,A[y]);

      return new_point;
    }
    
    template<typename R> inline 
    std::ostream& operator<<(std::ostream& os, const Point<R>& pt) {
      return pt.write(os);
    }
    
    template<typename R> inline
    std::istream& operator>>(std::istream& is, Point<R>& pt) {
      return pt.read(is);
    }


  }
}

#endif /* _ARIADNE_POINT_H */
