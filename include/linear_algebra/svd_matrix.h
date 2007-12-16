/***************************************************************************
 *            svd_matrix.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
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
 
/*! \file svd_matrix.h
 *  \brief Singular value decomposition.
 */

#ifndef ARIADNE_SVD_MATRIX_H
#define ARIADNE_SVD_MATRIX_H

#include <tblas/copy.hpp>
#include <tblas/geset.hpp>
#include <tlapack/gesvd.hpp>

#include "base/array.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*! \ingroup LinearAlgebra
     *  \brief A matrix stored in SVD product form. 
     */
    template<class Real>
    class SVDMatrix {
     public:
      /*! \brief Construct from an ordinary matrix. */
      SVDMatrix(const Matrix<Real>& A);
      
      /*! \brief The number of rows. */
      size_type number_of_rows() const { return _u.number_of_rows(); }
      /*! \brief The number of columns. */
      size_type number_of_columns() const { return _vt.number_of_rows(); }
      /*! \brief The \a i,\a j th element. */
      Real operator() (const size_type& i, const size_type& j) const;

      /*! \brief A vector containing the singular values in decreasing order. */
      const Vector<Real>& S() const { return _s; }
      /*! \brief The left orthogonal factor. */
      const Matrix<Real>& U() const { return _u; }
      /*! \brief The right orthogonal factor. */
      const Matrix<Real> V() const { return _vt.transpose(); }
      /*! \brief The transpose of the right orthogonal factor. */
      const Matrix<Real>& Vt() const { return _vt; }
      /*! \brief The diagonal matrix of singular values. */
      Matrix<Real> D() const;
      
      /*! \brief Convert to an ordinary matrix. */
      operator Matrix<Real> () const;
     private:
      Vector<Real> _s;
      Matrix<Real> _u;
      Matrix<Real> _vt;
    };
    
    template<class Real>
    inline 
    SVDMatrix<Real>::SVDMatrix(const Matrix<Real>& A) 
      : _s(TBLAS::min(A.number_of_rows(),A.number_of_columns())),
        _u(A.number_of_rows(),A.number_of_rows()),
        _vt(A.number_of_columns(),A.number_of_columns())
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();
      array<Real> work(A.begin(),A.begin()+m*n);
      TLAPACK::gesvd(TBLAS::RowMajor,m,n,work.begin(),n,_s.data().begin(),_u.data().begin(),m,_vt.data().begin(),n);
    }

    template<class Real>
    inline
    Matrix<Real>
    SVDMatrix<Real>::D() const
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();
      int nsv=TBLAS::min(m,n);
      Matrix<Real> result(m,n);
      TBLAS::geset(TBLAS::RowMajor,m,n,Real(0),Real(0),result.data().begin(),n);
      TBLAS::copy(nsv,_s.data().begin(),1,result.data().begin(),n+1);
      return result;
    }

    template<class Real>
    inline
    Real
    SVDMatrix<Real>::operator() (const size_type& i, const size_type& j) const
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();

      Real result=0;
      for(int k=0; k!=TBLAS::min(m,n); ++k) { 
        result += _u(i,k) * _s(k) * _vt(j,k); 
      } 
      return result; 
    }


    template<class Real>
    inline
    SVDMatrix<Real>::operator Matrix<Real> () const
    {
      int m=this->number_of_rows();
      int n=this->number_of_columns();

      Matrix<Real> result(m,n);
      for(int i=0; i!=m; ++i) {
        for(int j=0; j!=n; ++j) {
          result(i,j)=0;
          for(int k=0; k!=TBLAS::min(m,n); ++k) { 
            result(i,j) += _u(i,k) * _s(k) * _vt(j,k); 
          } 
        }
      }
      return result; 
    }

  }
}

#endif /* ARIADNE_SVD_MATRIX_H */
