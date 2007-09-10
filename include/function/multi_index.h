/***************************************************************************
 *            multi_index.h
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
 
/*! \file multi_index.h
 *  \brief Multi-indices of symmetric tensors.
 */

#ifndef ARIADNE_MULTI_INDEX_H
#define ARIADNE_MULTI_INDEX_H

#include <cassert>

#include "../base/iterator.h"
#include "../base/array.h"
#include "../base/stlio.h"

#include "../numeric/integer.h"

#include "sorted_index.h"

namespace Ariadne {
  namespace Function {
    
    class SortedIndex;
    class PositionIndex;

    /*! \ingroup LinearAlgebra
     *  \brief An index of a symmetric object. 
     */
    class MultiIndex {
     public:
      typedef Ariadne::Base::size_type size_type;
      typedef size_type value_type;
     
      /*! Construct a multi index of degree \a 0 with \a nv variables. */
      explicit MultiIndex(size_type nv);
      /*! Construct a multi index with \a nv variables from the array \a ary. */
      explicit MultiIndex(size_type nv, const size_type* ary);
      /*! Construct a multi index from the sorted index \a a. */
      MultiIndex(const SortedIndex& a);
      /*! Construct a multi index from the positional index \a a. */
      MultiIndex(const PositionIndex& a);
      /*! Copy constructor. */
      MultiIndex(const MultiIndex& a);
      /*! Copy assignment operator. */
      MultiIndex& operator=(const MultiIndex& a);

      /*! The degree of the multi-index, equal to the sum of the number of occurrences of the variables. */
      const size_type degree() const;
      /*! The number of variables. */
      const size_type number_of_variables() const;
      /*! The number of occurrences of the \a i th variable. */
      const size_type& operator[](const size_type& i) const;
     
      /*! Equality operator. */
      bool operator==(const MultiIndex& a) const;
      /*! Inequality operator. */
      bool operator!=(const MultiIndex& a) const;
      /*! Comparison operator. */
      bool operator<(const MultiIndex& a) const; // inline

      /*! Increment. No post-increment operator as we sometimes pass MultiIndex by reference. */
      MultiIndex& operator++(); // inline
      // No post-increment operator as we sometimes pass MultiIndex by reference.Post increment. 
      // MultiIndex operator++(int); 
      /*! Inplace sum. */
      MultiIndex& operator+=(const MultiIndex& a); // inline
      /*! Sum. */
      MultiIndex operator+(const MultiIndex& a) const; // inline
     
      /*! Convert to a normal tensor index, with elements ordered lowest to highest. */
      operator SortedIndex () const;
     
      /*! The position of the element in the array of tensor values. */
      size_type position() const;
      
      /*! The number of ordered index arrays with each element occurring the number of times specified by the multi index. */
      size_type number() const;
      
      /*! Set the value of the \a i th index to \a j. */
      void set_index(const size_type& i, const size_type j);
      /*! Increment the the \a i th index, thereby increasing the degree. */
      void increment_index(const size_type& i);
      /*! Decrement the the \a i th index, thereby decreasing the degree. */
      void decrement_index(const size_type& i);
      
      /*! Set the value of the \a i th index to \a j. */
      void set(const size_type& i, const size_type j);
      friend std::ostream& operator<<(std::ostream&, const MultiIndex&);
     private:
      size_type _degree;
      array<size_type> _occurrences;
    };
      


    
    inline MultiIndex::MultiIndex(size_type nv)
      : _degree(0u), _occurrences(nv,0u) 
    {
    }

    inline MultiIndex::MultiIndex(size_type nv, const size_type* ary)
      : _degree(0u), _occurrences(nv,0u) 
    {
      for(size_type i=0; i!=nv; ++i) { this->_occurrences[i]=ary[i]; this->_degree+=ary[i]; }
    }

    inline MultiIndex::MultiIndex(const SortedIndex& a)
      : _degree(0u), _occurrences(a.number_of_variables(),0u)
    { 
      for(size_type i=0; i!=a.degree(); ++i) { this->increment_index(a[i]); } 
    }

    inline MultiIndex::MultiIndex(const MultiIndex& a)
      : _degree(a._degree), _occurrences(a._occurrences) 
    {
    }

    inline MultiIndex& MultiIndex::operator=(const MultiIndex& a) { 
      this->_degree=a._degree; this->_occurrences=a._occurrences; return *this; }


    inline const size_type MultiIndex::degree() const { 
      return this->_degree; 
    }

    inline const size_type MultiIndex::number_of_variables() const { 
      return this->_occurrences.size(); 
    }

    inline const size_type& MultiIndex::operator[](const size_type& i) const {
      assert(i<this->number_of_variables()); return this->_occurrences[i]; 
    }
     

    inline bool MultiIndex::operator==(const MultiIndex& a) const { 
      return this->_occurrences==a._occurrences; 
    }
    
    inline bool MultiIndex::operator!=(const MultiIndex& a) const { 
      return !(*this==a); 
    }

    inline void MultiIndex::set_index(const size_type& i, const size_type j) { 
      this->_degree+=(j-this->_occurrences[i]); this->_occurrences[i]=j; 
    }

    inline void MultiIndex::increment_index(const size_type& i) { 
      ++this->_occurrences[i]; ++this->_degree; 
    }

    inline void MultiIndex::decrement_index(const size_type& i) { 
      if(!(this->_occurrences[i]>0)) { 
        throw std::runtime_error("MultiIndex::decrement_index: the number of occurence of the index must be positive"); 
      }
      --this->_occurrences[i]; 
      --this->_degree; 
    }
      
    inline void MultiIndex::set(const size_type& i, const size_type j) { 
      this->set_index(i,j); 
    }





    inline
    size_type MultiIndex::number() const
    {
      size_type result=Numeric::factorial(this->degree());
      for(size_type k=0; k!=this->number_of_variables(); ++k) {
        result/=Numeric::factorial((*this)[k]);
      }
      return result;
    }
      
    inline
    size_type MultiIndex::position() const
    {
      size_type deg=this->degree()-1;
      size_type nvar=this->number_of_variables();
      size_type result=Numeric::choose(deg+nvar,nvar);
      for(size_type k=0; k!=this->number_of_variables()-1; ++k) {
        --nvar;
        deg-=(*this)[k];
        result+=Numeric::choose(deg+nvar,nvar);
      }
      return result;
    }
      
    inline
    MultiIndex::operator SortedIndex () const 
    {
      SortedIndex result(this->number_of_variables(),this->degree()); 
      size_type k=0;
      for(size_type i=0; i!=this->number_of_variables(); ++i) {
        for(size_type j=k; j!=k+(*this)[i]; ++j) {
          result._entries[j]=i;
        }
        k=k+(*this)[i];
      }
      return result;
    }
      
    inline
    bool MultiIndex::operator<(const MultiIndex& a2) const {
      const MultiIndex& a1=*this;
 
      if(a1.number_of_variables()!=a2.number_of_variables()) {
        throw std::runtime_error("operator<(SortedIndex,SortedIndex): number of variables must match");
      }

      if(a1.degree()!=a2.degree()) {
        return a1.degree()<a2.degree();
      } else {
        for(size_type i=0; i!=a1.number_of_variables(); ++i) {
          if(a1[i]!=a2[i]) { 
            return a1[i]<a2[i];
          }
        }
        return false;
      }
    }


    inline
    MultiIndex& MultiIndex::operator++() 
    {
      MultiIndex& a=*this;
      const size_type& nv=a.number_of_variables();
      assert(nv>0);
      if(nv==1) {
        a.increment_index(0);
        return a;
      }
      if(a[nv-2]!=0) {
        a.decrement_index(nv-2);
        a.increment_index(nv-1);
        return a;
      } else {
        size_type li=a[nv-1];
        a.set_index(nv-1,0);
        for(size_type k=nv-1; k!=0; --k) {
          if(a[k-1]!=0) {
            a.decrement_index(k-1);
            a.set_index(k,li+1);
            return a;
          }
        }
        a.set_index(0,li+1);
      }
      return a;
    }
    
    inline
    MultiIndex& MultiIndex::operator+=(const MultiIndex& a2) {
      for(size_type i=0; i!=this->number_of_variables(); ++i) {
        this->set_index(i,(*this)[i]+a2[i]);
      }
      return *this;
    }
    
    inline
    MultiIndex MultiIndex::operator+(const MultiIndex& a2) const {
      MultiIndex result(a2.number_of_variables());
      const MultiIndex& a1=*this;
      for(size_type i=0; i!=a1.number_of_variables(); ++i) {
        result.set_index(i,a1[i]+a2[i]);
      }
      return result;
    }
    

    class MultiIndexIterator
      : public boost::iterator_facade<MultiIndexIterator,
                               MultiIndex,
                               boost::forward_traversal_tag,
                               const MultiIndex&,
                               const MultiIndex*>
    {
      MultiIndex _index;
     public:
      MultiIndexIterator(const MultiIndex& i) : _index(i) { }
      MultiIndexIterator(const size_type& nv, const size_type& d)
        : _index(MultiIndex(nv)) { this->_index.set_index(0,d); }
       
      void increment() {
        const size_type nv=this->_index.number_of_variables();
        if(this->_index[nv-2]!=0) {
          this->_index.decrement_index(nv-2);
          this->_index.increment_index(nv-1);
          return;
        } else {
          size_type li=this->_index[nv-1];
          this->_index.set_index(nv-1,0);
          for(size_type k=nv-1; k!=0; --k) {
            if(this->_index[k-1]!=0) {
              this->_index.decrement_index(k-1);
              this->_index.set_index(k,li+1);
              return;
            }
          }
          this->_index.set_index(0,li+1);
        }
      }
      
      const MultiIndex& dereference() const { return this->_index; }
      
      bool equal(const MultiIndexIterator& other) const { return this->_index==other._index; }
    };
              
      
            
      
      
    inline 
    std::ostream& operator<<(std::ostream& os, const MultiIndex& a) {
      return os << a._occurrences;
    }
    
  }
}
#endif /* ARIADNE_MULTI_INDEX_H */