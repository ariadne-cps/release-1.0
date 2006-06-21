/***************************************************************************
 *            zonotope.h
 *
 *  6 February 2006
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
 
/*! \file zonotope.h
 *  \brief Zonotopes (affine images of cuboids).
 */

#ifndef _ARIADNE_ZONOTOPE_H
#define _ARIADNE_ZONOTOPE_H

#include <iosfwd>

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../linear_algebra/transformation_system.h"

#include "../numeric/interval.h"

#include "../geometry/point.h"

namespace Ariadne {
  namespace Geometry {

    template<> 
    inline bool is_a<Zonotope,Zonotope>() { return true; }
    template<> 
    inline bool is_a<Zonotope,Polyhedron>() { return true; }

    /* Forward declaration of friends. */
    template<typename R> std::ostream& operator<<(std::ostream&, const Zonotope<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Zonotope<R>&);

    /*! \brief A zonotope of arbitrary dimension.
     * 
     * A zonotope is a set of the form \f$c+Ae\f$, where \f$||e||_{infty}\leq1\f$.
     * The intersection and membership tests are performed using algorithms from: <br>
     * Guibas, Leonidas J.; Nguyen, An; Zhang, Li, "Zonotopes as bounding volumes."  <i>Proceedings of the Fourteenth Annual ACM-SIAM Symposium on Discrete Algorithms</i> (Baltimore, MD, 2003),  803--812, ACM, New York, 2003.
     */
    template <typename R>
    class Zonotope {
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the rectangle. */
      typedef Point<R> state_type;
      /*! \brief The type of Vectors. */
      typedef Ariadne::LinearAlgebra::Vector<R> Vector_type;
      /*! \brief The type of Matrix giving principal directions. */
      typedef Ariadne::LinearAlgebra::Matrix<R> Matrix_type;

     private:
      /* Zonotope's centre. */
      state_type _centre;
      
      /* Zonotope's principal directions. */
      Matrix_type _generators;
     public:
      /*! \brief Default constructor constructs an empty zonotope of dimension \a n. */
      explicit Zonotope(size_type n = 0)
        : _centre(n),  _generators(n,0) { }
     
      /*! \brief Construct from centre and directions. */
      explicit Zonotope(const Vector_type& c, const Matrix_type& m)
        : _centre(c), _generators(m)
      {
        if (c.size()!=m.size1()) {
          throw std::domain_error("The centre and directions have different dimensions.");
        }
      }
      
      /*! \brief Construct from centre and directions. */
      explicit Zonotope(const state_type& c, const Matrix_type& m)
        : _centre(c), _generators(m)
      {
        if (c.dimension()!=m.size1()) {
          throw std::domain_error(
              "The the Matrix of principal directions does not have the same number of rows as the point dimension.");
        }

        this->minimize_generators();
      }
       
      /*! \brief Construct from a rectangle. */
      explicit Zonotope(const Rectangle<real_type>& r);

      /*! \brief Construct from a Parallelotope. */
      explicit Zonotope(const Parallelotope<real_type>& p);

      /*! \brief Construct from a string literal. */
      explicit Zonotope(const std::string& s);
      
      /*! \brief Copy constructor. */
      Zonotope(const Zonotope<R>& original)
        : _centre(original._centre),
          _generators(original._generators)
      { }
      
      /*! \brief Copy assignment operator. */
      Zonotope<R>& operator=(const Zonotope<R>& original) {
        if(this != &original) {
          this->_centre = original._centre;
          this->_generators = original._generators;
        }
        return *this;
      }
      
      /*! \brief The equality operator. */
      bool operator==(const Zonotope<real_type>& A) const;
      
      /*! \brief The inequality operator */
      bool operator!=(const Zonotope<real_type>& A) const {
        return !(*this == A);
      }

      /*! \brief The dimension of the Euclidean space the zonotope lies in. */
      dimension_type dimension() const {
        return this->_centre.dimension();
      }
      
      /*! \brief The centre of the zonotope. */
      const state_type &centre() const {
        return this->_centre;
      }
      
      /*! \brief The radius of the zonotope. */
      real_type radius() const {
        return this->bounding_box().radius();
      }
      
      /*! \brief The \a n th of principle direction. */
      Vector_type generator(size_type n) const {
        return column(this->_generators,n);
      }
      
      /*! \brief The Matrix of principle directions. */
      const Matrix_type &generators() const {
        return this->_generators;
      }
     
      /*! \brief The number of generators of the zonotope. */
      size_type number_of_generators() const {
        return this->_generators.size2();
      }
      
      /*! \brief True if the zonotope is empty. */
      bool empty() const {
        if ((this->_generators.size1()>0)&&
	    (this->_generators.size2()>0))
	    return false;

	return true;
      }
      
      /*! \brief True if the zonotope has empty interior. */
      bool empty_interior() const {
        using namespace Ariadne::LinearAlgebra;

        return !independent_rows(this->_generators);
      }
      
      /*! \brief A rectangle containing the given zonotope. */
      Rectangle<R> bounding_box() const;
      
      /*! \brief Returns the zonotope's vertices */
      std::vector< state_type > vertices() const;
      
      /*! \brief Tests if the zonotope contains \a point. */
      bool contains(const state_type& point) const;
     
      /*! \brief Tests if the zonotope contains a \a rectangle. */
      bool contains(const Rectangle<R>& rect) const; 
      
      /*! \brief Tests if the interior of the zonotope contains \a point. */
      bool interior_contains(const state_type& point) const;

      /*! \brief Subdivide into two smaller pieces. */
      ListSet<R,Geometry::Zonotope> divide() const;
      
      /*! \brief Subdivide into smaller pieces in each dimension. */
      ListSet<R,Geometry::Zonotope> subdivide() const;
      
      /*! \brief Convert to a polyhedron. */
      operator Polyhedron<R> () const;
      
      /*! \brief Construct a parallelopic over-approximation. */
      Parallelotope<R> over_approximating_parallelotope() const;
      
     private: 
/*      friend Zonotope<R> minkowski_sum <> (const Zonotope<R>& A, const Zonotope<R>& B);
      friend Zonotope<R> minkowski_difference <> (const Zonotope<R>& A, const Zonotope<R>& B);       */
      friend std::ostream& operator<< <> (std::ostream& os, const Zonotope<R>& r);
      friend std::istream& operator>> <> (std::istream& is, Zonotope<R>& r);
     private:
      // Minimize the generator Matrix
      void minimize_generators(void);
      
      // Order the generator Matrix by norm.
      void sort_generators(void);
      
      // The linear inequalities defining the zonotope.
      void compute_linear_inequalities(Matrix_type&, Vector_type&) const;

      // A possible vertex is the image of a vertex of the cube in 
      // generator space under the affine transformation
      state_type _possible_vertex(const size_type& vert_num) const;
      std::vector< state_type > _get_possible_vertices() const ;
    };
  
    
    /*! \brief Performs the Minkoswi sum of two zonotopes */
    template<typename R> 
    Zonotope<R> 
    minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B);

    /*! \brief Performs the Minkoswi difference of two zonotopes */
    template<typename R> 
    Zonotope<R> 
    minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B); 

    /*! \brief Performs the Minkoswi sum of a zonotope and a basic set */
    template<typename R, template <typename> class BS> 
    inline Zonotope<R> 
    minkowski_sum(const Zonotope<R>& A, const BS<R>& B) {
      return minkowski_sum(A, Zonotope<R>(B));
    }

    /*! \brief Performs the Minkoswi difference of a zonotope and a basic set */
    template<typename R, template <typename> class BS> 
    inline
    Zonotope<R> 
    minkowski_difference(const Zonotope<R>& A, const BS<R>& B)
    {
       return minkowski_difference(A, Zonotope<R>(B));
    }

    /*! \brief Performs the Minkoswi sum of a zonotope and a basic set */
    template<typename R, template <typename> class BS> 
    inline Zonotope<R> 
    minkowski_sum(const BS<R>& A, const Zonotope<R>& B) {
      return minkowski_sum(B, Zonotope<R>(A));
    }

    /*! \brief Performs the Minkoswi difference of a zonotope and a basic set */
    template<typename R, template <typename> class BS> 
    inline
    Zonotope<R> 
    minkowski_difference(const BS<R>& A, const Zonotope<R>& B)
    {
       return minkowski_difference(B, Zonotope<R>(A));
    }


    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool 
    disjoint(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return !minkowski_difference(A,B).contains(Point<R>(A.dimension()));
    }
    
    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool 
    disjoint(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return disjoint(Zonotope<R>(A),B);
    }

    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool 
    disjoint(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return disjoint(B,A);
    }

    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool 
    disjoint(const Parallelotope<R>& A, const Zonotope<R>& B) 
    {
      return disjoint(Zonotope<R>(A),B);
    }

    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool 
    disjoint(const Zonotope<R>& A, const Parallelotope<R>& B) 
    {
      return disjoint(B,A);
    }

    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool 
    interiors_intersect(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      //return !minkowski_sum(A,B).interior_contains(Point<R>(A.dimension()));
      
      return interiors_intersect(Polyhedron<R>(A),Polyhedron<R>(B));
      
    }
   
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool 
    interiors_intersect(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      Zonotope<R> z_A(A);
      return interiors_intersect(Zonotope<R>(A),B);
    }

    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool 
    interiors_intersect(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return interiors_intersect(B,A);
    }
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool 
    interiors_intersect(const Parallelotope<R>& A, const Zonotope<R>& B) 
    {
      return interiors_intersect(Zonotope<R>(A),B);
    }
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool 
    interiors_intersect(const Zonotope<R>& A, const Parallelotope<R>& B) 
    {
      return interiors_intersect(B,A);
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool 
    inner_subset(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool 
    inner_subset(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool 
    inner_subset(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool 
    inner_subset(const Parallelotope<R>& A, const Zonotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool 
    inner_subset(const Zonotope<R>& A, const Parallelotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }


    /*! \brief Tests inclusion of \a A in \a B. */
    template <typename R>
    inline
    bool 
    subset(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in \a B. */
    template <typename R>
    inline
    bool 
    subset(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool 
    subset(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool 
    subset(const Parallelotope<R>& A, const Zonotope<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool 
    subset(const Zonotope<R>& A, const Parallelotope<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }
    
    /*! \brief Scale a zonotope. */
    template<typename R>
    inline
    Geometry::Zonotope<R> 
    scale(const Geometry::Zonotope<R>& z, const R& scale_factor) {

      const Geometry::Point<R>& centre=z.centre();
      const LinearAlgebra::Matrix<R>& generators=z.generators();
      
      Geometry::Point<R> new_centre(z.dimension());

      for(size_type i=0; i!=z.dimension(); ++i) {
        new_centre[i]=scale_factor*centre[i];
      }

      return Geometry::Zonotope<R>(new_centre, scale_factor*generators);
    }

    
    template<typename R>
    Zonotope<R>
    operator+(const Rectangle<R>& r, const LinearAlgebra::TransformationSystem<R>& v);
    
    template<typename R>
    Zonotope<R>
    operator+(const Zonotope<R>& z, const LinearAlgebra::TransformationSystem<R>& v);
    

  }
}

#endif /* _ARIADNE_ZONOTOPE_H */
