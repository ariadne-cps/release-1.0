/***************************************************************************
 *            numerical_traits.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file numerical_traits.h
 *  \brief Traits classes to define properties of numerical types.
 */

#ifndef ARIADNE_NUMERICAL_TRAITS_H
#define ARIADNE_NUMERICAL_TRAITS_H

#include <string>
#include <gmpxx.h>

// FIXME: WE should not use GMP internals
// The typedef below  may be needed when using GMP 4.2 or above.
class __gmpq_value;

namespace Ariadne {
  namespace Numeric {
    class Integer;
    class Rational;
    class Float64;
    class FloatMP;
    template<class R> class Interval;
    template<class X, class V> class Differential;
      
    /* numerical traits */
    /*! \brief Tags a class representing a ring. */
    class ring_tag { };
    /*! \brief Tags a class representing a field. */
    class field_tag { };
    /*! \brief Tags a class based on integer numbers. */
    class integer_tag { };
    /*! \brief Tags a class based on rational numbers. */
    class rational_tag { };
    /*! \brief Tags a class based on floating-point numbers. */
    class float_tag { };
      
#ifndef DOXYGEN
    template<class T1, class T2=T1> struct traits { 
    };
    //template<class T1, class T2=T1> class traits { };

#else
    /*! \brief Typedef's describing the results of binary operations involving \a T1 and \a T2. */
    template<class T1, class T2> struct traits<T1,T2> { 
      /*!\brief The default type used to represent the result of a binary arithmetical operation. */
      typedef AT arithmetic_type;
    };
    
    /*! \brief Typedef's describing a numerical type. */
    template<class T> struct traits<T,T> { 
      /*!\brief The type of real number used for T. */
      typedef RT real_type;
      /*!\brief The type which can be assigned an element of \a T. */
      typedef T closure_type;
      /*!\brief The default type used to store the result of a binary arithmetical operation. */
      typedef AT arithmetic_type;
      /*!\brief The type used for interval arithmetic over \a T. Defaults to Interval<T> if \a T is
       * a number type, and T if \a T is an interval type.
       */
      typedef IT interval_type;
    };
#endif      

    template<> struct traits<int> {
      typedef integer_tag type;
      typedef int number_type;
      typedef int closure_type;
      typedef int arithmetic_type;
    };
  
    template<> struct traits<double> {
      typedef double number_type;
      typedef double closure_type;
      typedef double arithmetic_type;
      typedef Interval<double> interval_type;
    };
  
    template<> struct traits< mpf_class > { 
      typedef mpf_class number_type; 
      typedef mpf_class approximate_arithmetic_type; 
      typedef mpf_class closure_type; 
      typedef mpf_class arithmetic_type; 
      typedef Interval<FloatMP> interval_type; 
    };
    
    template<> struct traits< Float64 > { 
      typedef float_tag type;
      typedef Float64 number_type; 
      typedef double approximate_arithmetic_type; 
      typedef Float64 closure_type; 
      typedef Interval<Float64> arithmetic_type; 
      typedef Interval<Float64> interval_type;
    };
    
    template<> struct traits< FloatMP > { 
      typedef float_tag type;
      typedef FloatMP number_type; 
      typedef mpf_class approximate_arithmetic_type; 
      typedef FloatMP closure_type; 
      typedef Interval<FloatMP> arithmetic_type; 
      typedef Interval<FloatMP> interval_type; 
    };
    
    template<> struct traits< Rational > { 
      typedef rational_tag type;
      typedef Rational number_type; 
      typedef Rational approximate_arithmetic_type; 
      typedef Rational closure_type; 
      typedef Rational arithmetic_type; 
      typedef Interval<Rational> interval_type; 
    };

    // FIXME: WE should not use GMP internals
    // The following is needed for rational expressions in GMP 4.2.x
    template<class E> struct traits< __gmp_expr<mpq_t,E> > { 
      typedef Rational closure_type; 
    };

    // FIXME: WE should not use GMP internals
    // The following is needed for rational expressions in GMP 4.1.x
    template<class E> struct traits< __gmp_expr<__gmpq_value,E> > { 
      typedef Rational closure_type; 
    };
  


    template<class R> struct traits< Interval<R> > { 
      typedef R number_type; 
      typedef Interval<typename traits<R>::closure_type> closure_type; 
      typedef Interval<typename traits<R>::closure_type> arithmetic_type; 
      typedef Interval<R> interval_type; 
    };

    template<class X, class V> struct traits< Differential<X,V> > { 
      typedef typename traits<X>::number_type number_type; 
    };

    template<> struct traits< int, Rational > { 
      typedef Rational arithmetic_type; 
    };
    
    template<> struct traits< Rational, int > { 
      typedef Rational arithmetic_type; 
    };
    
    template<> struct traits< double, Rational > { 
      typedef Rational arithmetic_type; 
    };
    
    template<> struct traits< Rational, double > { 
      typedef Rational arithmetic_type; 
    };
    
    template<> struct traits< Float64, Rational > { 
      typedef Rational arithmetic_type; 
    };
    
    template<> struct traits< Rational, Float64 > { 
      typedef Rational arithmetic_type; 
    };
    
    template<> struct traits< FloatMP, Rational > { 
      typedef Rational arithmetic_type; 
    };
    
    template<> struct traits< Rational, FloatMP > { 
      typedef Rational arithmetic_type; 
    };
    


    template<class RF, class RI> struct traits< RF,Interval<RI> > { 
      typedef Interval<RI> arithmetic_type; 
    };
    
    template<class RF, class RI> struct traits< Interval<RI>, RF > { 
      typedef Interval<RI> arithmetic_type; 
    };    
    



    //! \name Numerical type description.
    //@{
    /*! \brief The name of class T. */
    template<class T> inline std::string name();

    //! \name Standard conversion operations. (Deprecated) 
    /*! \brief Approximate \a x by an element of Res. */
    template<class Res, class Arg> inline Res convert_to(const Arg& x) { return Res(x); }
    
    /*! \brief Approximate \a x by an element of Res with accuracy \a e. */
    template<class Res, class Arg, class Err> Res approximate(const Arg& x, const Err& e);
    //@}
  }   
}
  

#endif /* ARIADNE_NUMERICAL_TRAITS_H */
