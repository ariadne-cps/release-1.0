/***************************************************************************
 *            float64.h
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
 
/*! \file float64.h
 *  \brief Type definitions and conversion operators for 64-bit fixed precision floating point numbers.
 */

#ifndef _ARIADNE_FLOAT64_H
#define _ARIADNE_FLOAT64_H

#include <mpfr.h>

#include "../declarations.h"
#include "../numeric/numerical_traits.h"
#include "../numeric/function.h"
#include "../numeric/integer.h"
#include "../numeric/rational.h"

namespace Ariadne {
  namespace Numeric {

#ifdef DOXYGEN
    /*!\ingroup Numeric
     * \brief A 64-bit fixed-precision floating point number.
     *  
     * Standard operations are not exact, but must support interval arithmetic.
     *
     * Currently implemented by the built-in type double.
     */
  class Float64 { };
#else
    //typedef double Float64;
#endif
 
    class Float64 {
     public:
      Float64() : _value() { }
      Float64(const int& n) : _value(n) { }
      Float64(const uint& n) : _value(n) { }
      Float64(const double& x) : _value(x) { }
      Float64(const Float64& x) : _value(x._value) { }
      
      Float64& operator=(const int& n) { this->_value=n; return *this; }
      Float64& operator=(const double& x) { this->_value=x; return *this; }
      Float64& operator=(const Float64& x) { this->_value=x._value; return *this; }
      
      double get_d() const { return this->_value; }
      operator Rational () const { return mpq_class(this->_value); }
     private:
      double _value;
    };
    
    template<> class numerical_traits<Float64> {
     public:
      typedef ring_tag algebraic_category;
      typedef Rational field_extension_type;
      typedef Interval<Float64> arithmetic_type;
    };

    inline std::ostream& operator<<(std::ostream& os, const Float64& x) { return os << x.get_d(); }
    inline std::istream& operator>>(std::istream& is, Float64& x) { double d; is >> d; x=d; return is; }
    
    template<> inline std::string name<Numeric::Float64>() { return "Float64"; }
    template<> inline std::string name<Numeric::Interval<Numeric::Float64> >() { return "Interval<Float64>"; }
   
    inline bool operator==(const Float64& x1, const Float64& x2) { return x1.get_d()==x2.get_d(); }
    inline bool operator!=(const Float64& x1, const Float64& x2) { return x1.get_d()!=x2.get_d(); }
    inline bool operator<=(const Float64& x1, const Float64& x2) { return x1.get_d()<=x2.get_d(); }
    inline bool operator>=(const Float64& x1, const Float64& x2) { return x1.get_d()>=x2.get_d(); }
    inline bool operator< (const Float64& x1, const Float64& x2) { return x1.get_d()< x2.get_d(); }
    inline bool operator> (const Float64& x1, const Float64& x2) { return x1.get_d()> x2.get_d(); }
/*    
    bool operator==(const Float64& x1, const int& x2) { return x1.get_d()==x2; }
    bool operator!=(const Float64& x1, const int& x2) { return x1.get_d()!=x2; }
    bool operator<=(const Float64& x1, const int& x2) { return x1.get_d()<=x2; }
    bool operator>=(const Float64& x1, const int& x2) { return x1.get_d()>=x2; }
    bool operator< (const Float64& x1, const int& x2) { return x1.get_d()< x2; }
    bool operator> (const Float64& x1, const int& x2) { return x1.get_d()> x2; }
*/    
    int mpfr_hypot(mpfr_t y, const __mpfr_struct* x1, const __mpfr_struct* x2, mpfr_rnd_t r);

    typedef int mpfr_func(mpfr_t, const __mpfr_struct*, mp_rnd_t);
    typedef int mpfr_bin_func(mpfr_t, const __mpfr_struct*, const __mpfr_struct*, mp_rnd_t);
 
    inline
    void invoke_mpfr(Float64& y, const Float64& x, mpfr_func f, mp_rnd_t r)
    {
      mpfr_t xx;
      mpfr_init_set_d(xx, x.get_d(), r);
      f(xx, xx, r);
      y=mpfr_get_d(xx,r);
      mpfr_clear(xx);
    }
    
    inline 
    void invoke_mpfr(Float64& y, const Float64& x1, const Float64& x2, mpfr_bin_func f, mpfr_rnd_t r) 
    {
      mpfr_t xx1;
      mpfr_t xx2;
      mpfr_init_set_d(xx1, x1.get_d(), r);
      mpfr_init_set_d(xx2, x2.get_d(), r);
      f(xx1, xx1, xx2, r);
      y=mpfr_get_d(xx1, r);
      mpfr_clear(xx1);
      mpfr_clear(xx2);
    }

    inline
    Float64 invoke_mpfr(Float64 x, mpfr_func f, mp_rnd_t r)
    {
      Float64 y; invoke_mpfr(y, x, f, r); return y;
    }
     
    inline 
    Float64 invoke_mpfr(Float64 x1, Float64 x2, mpfr_bin_func f, mpfr_rnd_t r) 
    {
      Float64 y; invoke_mpfr(y, x1, x2, f, r); return y;
    }
    
    template<> inline Float64 min(const Float64& x1, const Float64& x2) { return (x1<=x2) ? x1 : x2; }
    template<> inline Float64 max(const Float64& x1, const Float64& x2) { return (x1>=x2) ? x1 : x2; }
    template<> inline Float64 neg(const Float64& x) { return Float64(-x.get_d()); }
    template<> inline Float64 abs(const Float64& x) { return (x>=0) ? x : neg(x); }
    
    template<> inline double conv_approx(const Float64& x) { return x.get_d(); }

    template<> inline Float64 conv_exact(const int& n) { return n; }
    template<> inline Float64 conv_approx(const int& n) { return conv_exact<Float64>(n); }
    template<> inline Float64 conv_down(const int& n) { return conv_exact<Float64>(n); }
    template<> inline Float64 conv_up(const int& n) { return conv_exact<Float64>(n); }

    template<> inline Float64 conv_exact(const double& x) { return x; }
    template<> inline Float64 conv_approx(const double& x) { return conv_exact<Float64>(x); }
    template<> inline Float64 conv_down(const double& x) { return conv_exact<Float64>(x); }
    template<> inline Float64 conv_up(const double& x) { return conv_exact<Float64>(x); }

    template<> inline Float64 min_exact(const Float64& x1,const Float64& x2) { 
      return min(x1,x2); }
    template<> inline Float64 min_approx(const Float64& x1,const Float64& x2) { 
      return min_exact(x1,x2); }
    template<> inline Float64 min_down(const Float64& x1,const Float64& x2) { 
      return min_exact(x1,x2); }
    template<> inline Float64 min_up(const Float64& x1,const Float64& x2) { 
      return min_exact(x1,x2); }
   
    template<> inline Float64 max_exact(const Float64& x1,const Float64& x2) { 
      return max(x1,x2); }
    template<> inline Float64 max_approx(const Float64& x1,const Float64& x2) { 
      return max_exact(x1,x2); }
    template<> inline Float64 max_down(const Float64& x1,const Float64& x2) { 
      return max_exact(x1,x2); }
    template<> inline Float64 max_up(const Float64& x1,const Float64& x2) { 
      return max_exact(x1,x2); }

    template<> inline Float64 neg_exact(const Float64& x) { return neg(x); }
    template<> inline Float64 neg_approx(const Float64& x) { return neg_exact(x); }
    template<> inline Float64 neg_down(const Float64& x) { return neg_exact(x); }
    template<> inline Float64 neg_up(const Float64& x) { return neg_exact(x); }
   
    template<> inline Float64 abs_exact(const Float64& x) { return abs(x); }
    template<> inline Float64 abs_approx(const Float64& x) { return abs_exact(x); }
    template<> inline Float64 abs_down(const Float64& x) { return abs_exact(x);  }
    template<> inline Float64 abs_up(const Float64& x) { return abs_exact(x);  }

    template<> inline Float64 add_down(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_add,GMP_RNDD); }
    template<> inline Float64 add_up(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_add,GMP_RNDU); }
    template<> inline Float64 add_approx(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_add,GMP_RNDN); }
    
    template<> inline Float64 sub_down(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_sub,GMP_RNDD); }
    template<> inline Float64 sub_up(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_sub,GMP_RNDU); }
    template<> inline Float64 sub_approx(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_sub,GMP_RNDN); }
    
    template<> inline Float64 mul_down(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_mul,GMP_RNDD); }
    template<> inline Float64 mul_up(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_mul,GMP_RNDU); }
    template<> inline Float64 mul_approx(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_mul,GMP_RNDN); }
    
    template<> inline Float64 div_down(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_div,GMP_RNDD); }
    template<> inline Float64 div_up(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_div,GMP_RNDU); }
    template<> inline Float64 div_approx(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_div,GMP_RNDN); }
    
    inline Float64 mul_approx(const int& n, const Float64& x) {
      return x.get_d()*n; }
    inline Float64 mul_approx(const Float64& x, const int& n) {
      return x.get_d()*n; }
    inline Float64 div_approx(const Float64& x, const int& n) {
      return x.get_d()/n; }
      
    inline Float64 mul_approx(const uint& n, const Float64& x) {
      return x.get_d()*n; }
    inline Float64 mul_approx(const Float64& x, const uint& n) {
      return x.get_d()*n; }
    inline Float64 div_approx(const Float64& x, const uint& n) {
      return x.get_d()/n; }
      
    inline Float64 mul_approx(const double& d, const Float64& x) {
      return x.get_d()*d; }
    inline Float64 mul_approx(const Float64& x, const double& d) {
      return x.get_d()*d; }
      
      
    template<> inline Float64 sqrt_down(const Float64& x) { 
      return invoke_mpfr(x,mpfr_sqrt,GMP_RNDD); }
    template<> inline Float64 sqrt_up(const Float64& x) { 
      return invoke_mpfr(x,mpfr_sqrt,GMP_RNDU); }
    template<> inline Float64 sqrt_approx(const Float64& x) { 
      return invoke_mpfr(x,mpfr_sqrt,GMP_RNDN); }

    template<> inline Float64 hypot_down(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_hypot,GMP_RNDD); }
    template<> inline Float64 hypot_up(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_hypot,GMP_RNDU); }
    template<> inline Float64 hypot_approx(const Float64& x1,const Float64& x2) {
      return invoke_mpfr(x1,x2,mpfr_hypot,GMP_RNDN); }
    
    template<> inline Integer int_down(const Float64& x);
    template<> inline Integer int_up(const Float64& x);

    template<> inline Float64 exp_approx(const Float64& x) {
      return invoke_mpfr(x,mpfr_exp,GMP_RNDN); }
    template<> inline Float64 exp_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_exp,GMP_RNDD); }
    template<> inline Float64 exp_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_exp,GMP_RNDU); };
    
    template<> inline Float64 log_approx(const Float64& x) {
      return invoke_mpfr(x,mpfr_log,GMP_RNDN); }
    template<> inline Float64 log_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_log,GMP_RNDD); }
    template<> inline Float64 log_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_log,GMP_RNDU); };

    template<> inline Float64 sin_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_sin,GMP_RNDD); }
    template<> inline Float64 sin_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_sin,GMP_RNDU); };

    template<> inline Float64 cos_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_cos,GMP_RNDD); }
    template<> inline Float64 cos_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_cos,GMP_RNDU); };

    template<> inline Float64 tan_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_tan,GMP_RNDD); }
    template<> inline Float64 tan_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_tan,GMP_RNDU); };

    template<> inline Float64 sinh_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_sinh,GMP_RNDD); }
    template<> inline Float64 sinh_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_sinh,GMP_RNDU); };

    template<> inline Float64 cosh_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_cosh,GMP_RNDD); }
    template<> inline Float64 cosh_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_cosh,GMP_RNDU); };

    template<> inline Float64 tanh_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_tanh,GMP_RNDD); }
    template<> inline Float64 tanh_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_tanh,GMP_RNDU); };

    template<> inline Float64 asin_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_asin,GMP_RNDD); }
    template<> inline Float64 asin_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_asin,GMP_RNDU); };

    template<> inline Float64 acos_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_acos,GMP_RNDD); }
    template<> inline Float64 acos_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_acos,GMP_RNDU); };

    template<> inline Float64 atan_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_atan,GMP_RNDD); }
    template<> inline Float64 atan_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_atan,GMP_RNDU); };

    template<> inline Float64 asinh_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_asinh,GMP_RNDD); }
    template<> inline Float64 asinh_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_asinh,GMP_RNDU); };

    template<> inline Float64 acosh_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_acosh,GMP_RNDD); }
    template<> inline Float64 acosh_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_acosh,GMP_RNDU); };

    template<> inline Float64 atanh_down(const Float64& x) {
      return invoke_mpfr(x,mpfr_atanh,GMP_RNDD); }
    template<> inline Float64 atanh_up(const Float64& x) {
      return invoke_mpfr(x,mpfr_atanh,GMP_RNDU); };
      
    
    inline Float64 operator-(const Float64& x) {
      return neg(x); }
      
    Interval<Float64> operator+(const Float64&, const Float64&);
    Interval<Float64> operator-(const Float64&, const Float64&);
    Interval<Float64> operator*(const Float64&, const Float64&);
    Interval<Float64> operator/(const Float64&, const Float64&);
      
  }

}

#endif /* _ARIADNE_FLOAT64_H */
