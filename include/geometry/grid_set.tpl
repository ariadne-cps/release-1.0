/***************************************************************************
 *            grid_set.tpl
 *
 *  10 January 2005
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "../base/arithmetic.h"

#include "../geometry/parallelotope.h"
#include "../geometry/partition_tree_set.h"
#include "../geometry/grid_set.h"

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g)
      : _grid(g), _position(g.dimension())
    { }

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const LatticeRectangle& b)
      : _grid(g), _position(b)
    {
      assert(g.dimension()==b.dimension());
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const IndexArray& l, const IndexArray& u)
      : _grid(g), _position(l,u)
    {
      assert(g.dimension()==l.size());
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const Grid<R>& g, const Rectangle<R>& r)
      : _grid(g), _position(g.dimension())
    {
      assert(g.dimension()==r.dimension());
      for(dimension_type i=0; i!=dimension(); ++i) {
        /* TODO: Catch and rethrow exceptions */
        _position.set_lower_bound(i,g.subdivision_index(i,r.lower_bound(i)));
        _position.set_upper_bound(i,g.subdivision_index(i,r.upper_bound(i)));
      }
    }

    template<typename R>
    GridRectangle<R>::GridRectangle(const GridCell<R>& c)
      : _grid(c._grid), _position(c.position())
    {
    }


    template<typename R>
    GridCell<R>::GridCell(const Grid<R>& g, const LatticeCell& pos)
      : _grid(g), _position(pos)
    {
      assert(_position.dimension()==_grid.dimension());
    }

    template<typename R>
    GridCell<R>::GridCell(const Grid<R>& g, const IndexArray& pos)
      : _grid(g), _position(pos)
    {
      assert(_position.dimension()==_grid.dimension());
    }



    template<typename R>
    GridCell<R>::operator Rectangle<R>() const {
      Rectangle<R> result(dimension());

      for(dimension_type i=0; i!=dimension(); ++i) {
        result.set_lower_bound(i, _grid.subdivision_coordinate(i,_position.lower_bound(i)));
        result.set_upper_bound(i, _grid.subdivision_coordinate(i,_position.upper_bound(i)));
      }

      return result;
    }


    template<typename R>
    GridRectangle<R>::operator Rectangle<R>() const {
      Rectangle<R> result(dimension());

      for(size_type i=0; i!=dimension(); ++i) {
        result.set_lower_bound(i, _grid.subdivision_coordinate(i,_position.lower_bound(i)));
        result.set_upper_bound(i, _grid.subdivision_coordinate(i,_position.upper_bound(i)));
      }

      return result;
    }




    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const GridRectangle<R>& r) {
      return os << Rectangle<R>(r);
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const GridCell<R>& c) {
      return os << Rectangle<R>(c);
    }



    /*! \brief Tests intersection of interiors */
    template<typename R>
    bool
    interiors_intersect(const GridRectangle<R>& A, const GridRectangle<R>& B) {
      assert(A.grid() == B.grid());
      return interiors_intersect(A.position(),B.position());
    }

    /*! \brief Tests intersection of interiors */
    template<typename R>
    bool
    interiors_intersect(const GridRectangle<R>& A, const GridMaskSet<R>& B) {
      assert(A.grid() == B.grid());
      return interiors_intersect(A._position,B._lattice_set);
    }

    /*! \brief Tests intersection of interiors */
    template<typename R>
    bool
    interiors_intersect(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.grid()==B.grid());
      return interiors_intersect(A._lattice_set,B._lattice_set);
    }



    /*! \brief Tests if A is a subset of B. */
    template<typename R>
    bool
    subset(const GridRectangle<R>& A, const GridRectangle<R>& B)
    {
      assert(A.dimension() == B.dimension());
      if(A.grid()==B.grid()) {
        return subset(A.position(),B.position());
      }
      return subset(Rectangle<R>(A),Rectangle<R>(B));
    }

    /*! \brief Tests if A is a subset of B. */
    template<typename R>
    bool
    subset(const GridRectangle<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.grid()==B.grid());
      return subset(A.position(),B.lattice_set()); 
    }

    /*! \brief Tests if A is a subset of B. */
    template<typename R>
    bool
    subset(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.grid()==B.grid());
      return subset(A.lattice_set(),B.lattice_set());
    }



    /*! \brief Tests intersection of interiors */
    template<typename R>
    bool
    interiors_intersect(const Rectangle<R>& A, const GridMaskSet<R>& B) {
      assert(A.dimension()==B.dimension());
      GridRectangle<R> gr=over_approximation_of_intersection(A,B.bounding_box(),B.grid());
      return interiors_intersect(gr,B);
    }
    
    /*! \brief Tests if A is a subset of B. */
    template<typename R>
    bool
    subset(const Rectangle<R>& r, const GridMaskSet<R>& B)
    {
      assert(r.dimension() == B.dimension());
      GridMaskSet<R> A(B.grid(),B.bounds());
      A.adjoin(over_approximation(r,A.grid()));
      return subset(A,B);
    }



    /*! \brief The union of \a A and \a B. */
    template<typename R>
    GridMaskSet<R>
    join(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.grid()==B.grid() && A.bounds()==B.bounds());
      return GridMaskSet<R>(A.grid(), A.bounds(), A.mask() | B.mask());
    }

    /*! \brief The closure of the intersection of the interior of \a A with the interior of \a B. */
    template<typename R>
    GridMaskSet<R>
    regular_intersection(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.grid()==B.grid() && A.bounds()==B.bounds());
      return GridMaskSet<R>(A.grid(), A.bounds(), A.mask() & B.mask());
    }

    /*! \brief The closure of the intersection of the interior of \a A with the interior of \a B. */
    template<typename R>
    GridMaskSet<R>
    difference(const GridMaskSet<R>& A, const GridMaskSet<R>& B)
    {
      assert(A.grid()==B.grid() && A.bounds()==B.bounds());
      return GridMaskSet<R>(A.grid(), A.bounds(), A.mask() - B.mask());
    }




    template<typename R>
    GridMaskSet<R>::operator ListSet<R,Rectangle>() const
    {
      //std::cerr << "GridMaskSet<R>::operator ListSet<R,Rectangle>() const" << std::endl;
      ListSet<R,Rectangle> result(this->dimension());
      for(typename GridMaskSet::const_iterator riter=begin(); riter!=end(); ++riter) {
        Rectangle<R> r(*riter);
        result.push_back(r);
      }
      return result;
    }


    template<typename R>
    GridCellListSet<R>::operator ListSet<R,Rectangle>() const
    {
      //std::cerr << "GridCellListSet<R>::operator ListSet<R,Rectangle>() const" << std::endl;
      ListSet<R,Rectangle> result(dimension());
      for(size_type i=0; i!=size(); ++i) {
        result.push_back((*this)[i]);
      }
      return result;
    }

    template<typename R>
    GridRectangleListSet<R>::operator ListSet<R,Rectangle>() const
    {
      //std::cerr << "GridRectangleListSet<R>::operator ListSet<R,Rectangle>() const" << std::endl;
      ListSet<R,Rectangle> result(dimension());
      for(size_type i=0; i!=size(); ++i) {
        result.push_back((*this)[i]);
      }
      return result;
    }



  
    template<typename R>
    GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& g)
      : _grid_ptr(&g), _lattice_set(g.bounds()) 
    { 
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const LatticeRectangle& b)
      : _grid_ptr(&g), _lattice_set(b) 
    {
      const FiniteGrid<R>* finite_grid_ptr=dynamic_cast<const FiniteGrid<R>*>(_grid_ptr);
      if(finite_grid_ptr) {
        assert(subset(b,finite_grid_ptr->bounds()));
      }
    }
  
    template<typename R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const LatticeMaskSet& ms)
      : _grid_ptr(&g), _lattice_set(ms)
    {
      const FiniteGrid<R>* finite_grid_ptr=dynamic_cast<const FiniteGrid<R>*>(_grid_ptr);
      if(finite_grid_ptr) {
        assert(subset(ms.bounds(),finite_grid_ptr->bounds()));
      }
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const LatticeRectangle& b, const BooleanArray& m)
      : _grid_ptr(&g), _lattice_set(b,m)
    {
      const FiniteGrid<R>* finite_grid_ptr=dynamic_cast<const FiniteGrid<R>*>(_grid_ptr);
      if(finite_grid_ptr) {
        assert(subset(b,finite_grid_ptr->bounds()));
      }
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const GridMaskSet<R>& gms) 
      : _grid_ptr(&gms.grid()), _lattice_set(gms._lattice_set)
    {
    }
    
    template<typename R>
    GridMaskSet<R>::GridMaskSet(const GridCellListSet<R>& gcls)
      : _grid_ptr(&gcls.grid()),
        _lattice_set(dynamic_cast<const FiniteGrid<R>*>(_grid_ptr)->bounds())
    {
      _lattice_set.adjoin(gcls._lattice_set);
    }
    
    template<typename R>
    GridMaskSet<R>::GridMaskSet(const GridRectangleListSet<R>& grls)
      : _grid_ptr(&grls.grid()),
        _lattice_set(dynamic_cast<const FiniteGrid<R>*>(_grid_ptr)->bounds())
    {
      _lattice_set.adjoin(grls._lattice_set);
    }

    template<typename R>
    GridMaskSet<R>::GridMaskSet(const ListSet<R,Rectangle>& rls) 
      : _grid_ptr(new FiniteGrid<R>(rls)), 
        _lattice_set(dynamic_cast<const FiniteGrid<R>*>(_grid_ptr)->bounds())
    {
      //FIXME: Memory leak!    

      for(typename ListSet<R,Rectangle>::const_iterator riter=rls.begin(); 
          riter!=rls.end(); ++riter) 
      {
        GridRectangle<R> r(grid(),*riter);
        adjoin(r);
      }
    }    
    
  



    template<typename R>
    GridCellListSet<R>::GridCellListSet(const Grid<R>& g)
      : _grid_ptr(&g), _lattice_set(g.dimension())
    {
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridMaskSet<R>& gms)
      : _grid_ptr(&gms.grid()), _lattice_set(gms.dimension())
    {
      this->_lattice_set.adjoin(gms._lattice_set);
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridCellListSet<R>& gcls)
      : _grid_ptr(&gcls.grid()), _lattice_set(gcls._lattice_set)
    {
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const GridRectangleListSet<R>& grls)
      : _grid_ptr(&grls.grid()), _lattice_set(grls.dimension())
    {
      this->_lattice_set.adjoin(grls._lattice_set); 
    }

    template<typename R>
    GridCellListSet<R>::GridCellListSet(const ListSet<R,Rectangle>& rls)
      : _grid_ptr(0), _lattice_set(rls.dimension())
    {
      (*this)=GridCellListSet<R>(GridRectangleListSet<R>(rls));      
    }




    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g)
      : _grid_ptr(&g), _lattice_set(g.dimension())
    {
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const ListSet<R,Rectangle>& s)
      : _grid_ptr(new FiniteGrid<R>(s)), _lattice_set(s.dimension())
    {
      //std::cerr << "GridRectangleListSet<R>::GridRectangleListSet(const ListSet<R,Rectangle>& s)" << std::endl;
      typedef typename ListSet<R,Rectangle>::const_iterator list_set_const_iterator;
      for(list_set_const_iterator iter=s.begin(); iter!=s.end(); ++iter) {
        this->adjoin(GridRectangle<R>(grid(),*iter));
      }
    }

    /* FIXME: This constructor is only included since Boost Python doesn't find the conversion operator */
    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const PartitionTreeSet<R>& pts)
      : _grid_ptr(0),_lattice_set(pts.dimension()) 
    {
      //std::cerr << "GridRectangleListSet<R>::GridRectangleListSet(const PartitionTreeSet<R>& pts)" << std::endl;
      SizeArray sizes=pts.depths();
      for(dimension_type i=0; i!=pts.dimension(); ++i) {
        sizes[i]=pow(2,sizes[i]);
      }
      _grid_ptr=new FiniteGrid<R>(pts.bounding_box(),sizes);
      for(typename PartitionTreeSet<R>::const_iterator iter=pts.begin(); iter!=pts.end(); ++iter) {
        Rectangle<R> rect(*iter);
        GridRectangle<R> grect(this->grid(),rect);
        this->adjoin(grect);
      }
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const GridRectangleListSet<R>& s)
      : _grid_ptr(&s.grid()), _lattice_set(s._lattice_set)
    {
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const GridCellListSet<R>& s)
      : _grid_ptr(&g), _lattice_set(g.dimension())
    {
      //TODO: re-write this function
      assert(false);
      assert(g.dimension()==s.dimension());

      const FiniteGrid<R>* to=dynamic_cast<const FiniteGrid<R>*>(&g);
      const FiniteGrid<R>* from=dynamic_cast<const FiniteGrid<R>*>(&s.grid());

      array< std::vector<index_type> > tr = FiniteGrid<R>::index_translation(*from,*to);
      //translate_cell_coordinates(&_lattice_set, s._lattice_set, tr);
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const GridRectangleListSet<R>& s)
      : _grid_ptr(&g), _lattice_set(s._lattice_set)
    {
      //TODO: re-write this function
      assert(false);

      const FiniteGrid<R>* to=dynamic_cast<const FiniteGrid<R>*>(&g);
      const FiniteGrid<R>* from=dynamic_cast<const FiniteGrid<R>*>(&s.grid());

      array< std::vector<index_type> > tr = FiniteGrid<R>::index_translation(*from,*to);
      //translate_rectangle_coordinates(&_lattice_set, s._lattice_set, tr);
    }

    template<typename R>
    GridRectangleListSet<R>::GridRectangleListSet(const Grid<R>& g, const ListSet<R,Rectangle>& ls)
      : _grid_ptr(&g), _lattice_set(ls.dimension())
    {
      typedef typename ListSet<R,Rectangle>::const_iterator ListSet_const_iterator;

      for(ListSet_const_iterator riter=ls.begin(); riter!=ls.end(); ++riter) {
        _lattice_set.adjoin(LatticeRectangle(grid().lower_index(*riter),grid().upper_index(*riter)));
      }
    }

    template<typename R>
    GridRectangle<R>
    over_approximation(const Rectangle<R>& r, const Grid<R>& g) 
    {
      if(r.empty()) {
        return GridRectangle<R>(g);
      }
      IndexArray lower(r.dimension());
      IndexArray upper(r.dimension());
      for(size_type i=0; i!=r.dimension(); ++i) {
        lower[i]=g.subdivision_lower_index(i,r.lower_bound(i));
        upper[i]=g.subdivision_upper_index(i,r.upper_bound(i));
      }
      return GridRectangle<R>(g,lower,upper);
    }
    
    
    template<typename R>
    GridCellListSet<R>
    over_approximation(const Parallelotope<R>& p, const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      assert(g.dimension()==p.dimension());
      if(p.empty()) {
        return result; 
      }
      Rectangle<R> bb=p.bounding_box();
      GridRectangle<R> gbb=over_approximation(bb,g);
      LatticeRectangle block=gbb.position();

      for(LatticeRectangle::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
        GridCell<R> cell(g,*iter);
        if(!disjoint(Rectangle<R>(cell),p)) {
          result.adjoin(cell);
        }
      }
      return result;
    }

    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    over_approximation(const ListSet<R,BS>& ls, const FiniteGrid<R>& g) 
    {
      //std::cerr << "GridMaskSet<R>::over_approximation(const ListSet<R,BS>& ls, const FiniteGrid<R>& g) " << std::endl;
      GridMaskSet<R> result(g);
      
      assert(g.dimension()==ls.dimension());

      for(typename ListSet<R,BS>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
       result.adjoin(over_approximation(*iter,g));
      }
        
      return result;
    }

    template<typename R>
    GridRectangle<R>
    over_approximation_of_intersection(const Rectangle<R>& r1, 
                                       const Rectangle<R>& r2,
                                       const Grid<R>& g) 
    {
      return over_approximation(intersection(r1,r2),g);
    }

    template<typename R>
    GridCellListSet<R>
    over_approximation_of_intersection(const Parallelotope<R>& p, 
                                       const Rectangle<R>& r,
                                       const Grid<R>& g) 
    {
      GridCellListSet<R> result(g);
      assert(p.dimension()==r.dimension());
      assert(g.dimension()==p.dimension());
      Rectangle<R> bb=intersection(p.bounding_box(),r);

      if(!bb.empty()) {
        GridRectangle<R> gbb=over_approximation(bb,g);
        LatticeRectangle block=gbb.position();

        for(LatticeRectangle::const_iterator iter=block.begin(); iter!=block.end(); ++iter) {
          GridCell<R> cell(g,*iter);
          if(!disjoint(Rectangle<R>(cell),p)) {
            result.adjoin(cell);
          }
        }
      }
      return result;
    }

    template<typename R, template<typename> class BS>
    GridMaskSet<R>
    over_approximation_of_intersection(const ListSet<R,BS>& ls, const Rectangle<R>& r, const FiniteGrid<R>& g) 
    {
      //std::cerr << "GridMaskSet<R>::over_approximation(const ListSet<R,BS>& ls, const FiniteGrid<R>& g) " << std::endl;
      GridMaskSet<R> result(g);
      
      assert(g.dimension()==ls.dimension());

      for(typename ListSet<R,BS>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
       result.adjoin(over_approximation_of_intersection(*iter,r,g));
      }
        
      return result;
    }

  

    template <typename R>
    std::ostream& operator<<(std::ostream& os,
                             const GridRectangleListSet<R>& set)
    {
      os << "GridRectangleListSet<" << name<R>() << ">(\n  rectangle_lattice_set: [\n    ";
      for (size_type i=0; i!=set.size(); i++) {
        if(i!=0) {
          os << ",\n    ";
        }
        os << Rectangle<R>(set[i]);
      }
      os << "\n  ]\n";
      os << "  grid: " << set.grid() << "\n";
      os << "  integer_rectangle_lattice_set: [ ";
      os.flush();
      for(size_type i=0; i!=set.size(); ++i) {
        if(i!=0) { os << ", "; }
        os << set[i].position();
      }
      os << " ]\n";
      os << ")" << std::endl;
      return os;
    }



    template <typename R>
    std::ostream& operator<<(std::ostream& os,
                             const GridCellListSet<R>& set)
    {
      os << "GridCellListSet<" << name<R>() << ">(\n";
      /*
      os << "  rectangle_lattice_set: [\n    ";
      if (!set.empty() ) {
        os << set[0];
      }
      for (size_t i=1; i<set.size(); ++i) {
        os << ",\n    " << Rectangle<R>(set[i]);
      }
      os << "\n  ]\n";
      */
      os << "  grid: " << set.grid() << "\n";
      os << "  integer_cell_lattice_set: [ ";
      os.flush();
      for(size_type i=0; i!=set.size(); ++i) {
        if(i!=0) { os << ", "; }
        os <<  set[i].position();
      }
      os << " ]\n";
      os << ")" << std::endl;
      return os;
    }

    template <typename R>
    std::ostream& operator<<(std::ostream& os,
                             const GridMaskSet<R>& set)
    {
      os << "GridMaskSet<" << name<R>() << ">(\n";
      os << "  grid: " << set.grid() << "\n";
      os << "  mask: " << set.mask() << "\n";
      os << ")\n";
      return os;
    }

  }
}