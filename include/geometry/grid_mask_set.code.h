/***************************************************************************
 *            grid_mask_set.code.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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

#include "grid_mask_set.h"

#include <ostream>

#include "../base/stlio.h"

#include "../combinatoric/array_operations.h"

#include "../geometry/rectangle.h"
#include "../geometry/zonotope.h"
#include "../geometry/polytope.h"
#include "../geometry/polyhedron.h"

#include "../geometry/grid_cell.h"
#include "../geometry/grid_block.h"
#include "../geometry/grid_cell_list_set.h"

#include "../geometry/list_set.h"
#include "../geometry/partition_tree_set.h"

#include "../geometry/set_interface.h"

#include "../output/logging.h"





namespace Ariadne {


template<class R>
Geometry::GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& fg)
  : _grid_ptr(new Grid<R>(fg.grid())), _lattice_set(fg.lattice_block()) 
{ 
}


template<class R>
Geometry::GridMaskSet<R>::GridMaskSet(const FiniteGrid<R>& fg, const BooleanArray& m)
  : _grid_ptr(new Grid<R>(fg.grid())), _lattice_set(fg.lattice_block(),m)
{
}


template<class R>
Geometry::GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Rectangle<R>& bb)
  : _grid_ptr(new Grid<R>(g)), _lattice_set(g.index_block(bb)) 
{ 
}

template<class R>
Geometry::GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b)
  : _grid_ptr(new Grid<R>(g)), _lattice_set(b) 
{ 
}


template<class R>
Geometry::GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeBlock& b, const BooleanArray& m)
  : _grid_ptr(new Grid<R>(g)), _lattice_set(b,m)
{
}


template<class R>
Geometry::GridMaskSet<R>::GridMaskSet(const Grid<R>& g, const Combinatoric::LatticeMaskSet& ms)
  : _grid_ptr(new Grid<R>(g)), _lattice_set(ms)
{
}


template<class R>
Geometry::GridMaskSet<R>::GridMaskSet(const GridMaskSet<R>& gms) 
  : _grid_ptr(gms._grid_ptr), _lattice_set(gms._lattice_set)
{
}


template<class R>
Geometry::GridMaskSet<R>&
Geometry::GridMaskSet<R>::operator=(const GridMaskSet<R>& gms) 
{
  if(this!=&gms) {
    this->_grid_ptr=gms._grid_ptr;
    this->_lattice_set=gms._lattice_set;
  }
  return *this;
}


template<class R>
Geometry::GridMaskSet<R>::GridMaskSet(const GridCellListSet<R>& gcls)
  : _grid_ptr(gcls._grid_ptr),
    _lattice_set(gcls.lattice_set())
{
}



template<class R>
Geometry::FiniteGrid<R>
Geometry::GridMaskSet<R>::finite_grid() const
{
  return FiniteGrid<R>(this->grid(),this->block());
}


template<class R>
Geometry::Rectangle<R>
Geometry::GridMaskSet<R>::extent() const
{
  return this->bounds();
}


template<class R> inline
Geometry::GridBlock<R> 
Geometry::GridMaskSet<R>::bounds() const 
{
  return GridBlock<R>(*this->_grid_ptr,_lattice_set.block()); 
}


// FIXME: Memory leak

template<class R>
Geometry::GridMaskSet<R>*
Geometry::GridMaskSet<R>::clone() const
{
  return new GridMaskSet<R>(*this);
}


template<class R>
tribool
Geometry::GridMaskSet<R>::contains(const Point<R>& pt) const
{
  return !Geometry::disjoint(*this,Rectangle<R>(pt));
}


template<class R>
tribool
Geometry::GridMaskSet<R>::superset(const Rectangle<R>& r) const
{
  return Geometry::subset(r,*this);
}


template<class R>
tribool
Geometry::GridMaskSet<R>::intersects(const Rectangle<R>& r) const
{
  return !Geometry::disjoint(*this,r);
}


template<class R>
tribool
Geometry::GridMaskSet<R>::disjoint(const Rectangle<R>& r) const
{
  return Geometry::disjoint(*this,r);
}


template<class R>
tribool
Geometry::GridMaskSet<R>::subset(const Rectangle<R>& r) const
{
  return Geometry::subset(*this,r);
}


template<class R> 
Geometry::Rectangle<R> 
Geometry::GridMaskSet<R>::bounding_box() const 
{
  return GridBlock<R>(grid(),bounds()); 
}


template<class R>
void
Geometry::GridMaskSet<R>::clear()
{
  this->_lattice_set.clear();
}




template<class R>
void
Geometry::GridMaskSet<R>::_instantiate_geometry_operators()
{
  typedef Numeric::Interval<R> I;
  tribool tb;
  //Point<I>* ipt=0;
  Rectangle<R>* r=0;
  
  ListSet< Rectangle<R> >* rls=0;
  ListSet< Zonotope<R,R> >* zls=0;
  ListSet< Zonotope<I,R> >* ezls=0;
  ListSet< Zonotope<I,I> >* izls=0;
  
  SetInterface<R>* set=0;
  
  FiniteGrid<R>* fg=0;
  GridCell<R>* gc=0;
  GridBlock<R>* gb=0;
  GridCellListSet<R>* gcls=0;
  GridMaskSet<R>* gms=0;
  
  *gms=outer_approximation(*rls,*fg);
  *gms=outer_approximation(*zls,*fg);
  *gms=outer_approximation(*ezls,*fg);
  *gms=outer_approximation(*izls,*fg);
  
  *gms=outer_approximation(*set,*fg);
  *gms=inner_approximation(*set,*fg);
  *rls=lower_approximation(*set,*fg);
  
  tb=Geometry::subset(*r,*gms);
  tb=Geometry::subset(*gms,*r);
  tb=Geometry::disjoint(*r,*gms);
  tb=Geometry::disjoint(*gms,*r);
  
  tb=Geometry::overlap(*gb,*gms);
  tb=Geometry::overlap(*gms,*gb);
  tb=Geometry::overlap(*gcls,*gms);
  tb=Geometry::overlap(*gms,*gcls);
  tb=Geometry::overlap(*gms,*gms);
  
  tb=Geometry::subset(*gc,*gms);
  tb=Geometry::subset(*gb,*gms);
  tb=Geometry::subset(*gcls,*gms);
  tb=Geometry::subset(*gms,*gms);
  
  *gms=Geometry::regular_intersection(*gb,*gms);
  *gms=Geometry::regular_intersection(*gms,*gb);
  *gcls=Geometry::regular_intersection(*gcls,*gms);
  *gcls=Geometry::regular_intersection(*gms,*gcls);
  *gms=Geometry::regular_intersection(*gms,*gms);
  *gcls=Geometry::difference(*gcls,*gms);
  *gms=Geometry::difference(*gms,*gms);
  *gms=Geometry::join(*gms,*gms);
}




// Geometric predicates ---------------------------------------------------



template<class R>
tribool
Geometry::disjoint(const GridBlock<R>& gb, const GridMaskSet<R>& gms) {
  ARIADNE_CHECK_SAME_GRID(gb,gms,"tribool disjoint(GridBlock gb, GridMaskSet gms)");
  return disjoint(gb.lattice_set(),gms.lattice_set());
}


template<class R>
tribool
Geometry::disjoint(const GridMaskSet<R>& gms, const GridBlock<R>& gb) {
  ARIADNE_CHECK_SAME_GRID(gms,gb,"tribool disjoint(GridMaskSet gms, GridBlock gb)");
  return disjoint(gms.lattice_set(),gb.lattice_set());
}


template<class R>
tribool
Geometry::disjoint(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
{
  ARIADNE_CHECK_SAME_GRID(gms1,gms2,"tribool disjoint(GridMaskSet gms1, GridMaskSet gms2)");
  return disjoint(gms1.lattice_set(),gms2.lattice_set());
}


template<class R>
tribool
Geometry::disjoint(const Rectangle<R>& r, const GridMaskSet<R>& gms) 
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r,gms,"tribool disjoint(Rectangle r, GridMaskSet gms)");
  Rectangle<R> br=closed_intersection(r,Rectangle<R>(gms.bounding_box()));
  GridBlock<R> gb=outer_approximation(br,gms.grid());
  return !overlap(gb,gms);
}


template<class R>
tribool
Geometry::disjoint(const GridMaskSet<R>& gms, const Rectangle<R>& r) {
  return disjoint(r,gms);
}





template<class R>
tribool
Geometry::overlap(const GridBlock<R>& gb, const GridMaskSet<R>& gms) {
  ARIADNE_CHECK_SAME_GRID(gb,gms,"tribool overlap(GridBlock gb, GridMaskSet gms)");
  return overlap(gb.lattice_set(),gms.lattice_set());
}


template<class R>
tribool
Geometry::overlap(const GridMaskSet<R>& gms, const GridBlock<R>& gb) {
  ARIADNE_CHECK_SAME_GRID(gms,gb,"tribool overlap(GridMaskSet gms, GridBlock gb)");
  return overlap(gms.lattice_set(),gb.lattice_set());
}


template<class R>
tribool
Geometry::overlap(const GridCellListSet<R>& A, const GridMaskSet<R>& B) {
  ARIADNE_CHECK_SAME_GRID(A,B,"overlap(GridCellListSet<R>,GridMaskSet<R>)");
  return overlap(A.lattice_set(),B.lattice_set());
}


template<class R>
tribool
Geometry::overlap(const GridMaskSet<R>& A, const GridCellListSet<R>& B) {
  ARIADNE_CHECK_SAME_GRID(A,B,"overlap(GridMaskSet<R>,GridCellListSet<R>)");
  return overlap(A.lattice_set(),B.lattice_set());
}


template<class R>
tribool
Geometry::overlap(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
{
  ARIADNE_CHECK_SAME_GRID(gms1,gms2,"tribool overlap(GridMaskSet gms2, GridMaskSet gms2)");
  return overlap(gms1.lattice_set(),gms2.lattice_set());
}






template<class R>
tribool
Geometry::subset(const GridMaskSet<R>& gms, const GridBlock<R>& gb)
{
  ARIADNE_CHECK_SAME_GRID(gms,gb,"tribool subset(GridMaskSet gms, GridBlock gb)");
  return subset(gms.lattice_set(),gb.lattice_set());
}

template<class R>
tribool
Geometry::subset(const GridCell<R>& gc, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_SAME_GRID(gc,gms,"tribool subset(GridCell gc, GridMaskSet gms)");
  return subset(gc.lattice_set(),gms.lattice_set()); 
}

template<class R>
tribool
Geometry::subset(const GridBlock<R>& gb, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_SAME_GRID(gb,gms,"tribool subset(GridBlock gb, GridMaskSet gms)");
  return subset(gb.lattice_set(),gms.lattice_set()); 
}

template<class R>
tribool
Geometry::subset(const GridCellListSet<R>& gcls, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_SAME_GRID(gcls,gms,"tribool subset(GridCellListSet gcls, GridMaskSet gms)");
  return subset(gcls.lattice_set(),gms.lattice_set()); 
}

template<class R>
tribool
Geometry::subset(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
{
  ARIADNE_CHECK_SAME_GRID(gms1,gms2,"tribool subset(GridMaskSet gms1, GridMaskSet gms2)");
  return subset(gms1.lattice_set(),gms2.lattice_set());
}




template<class R>
tribool
Geometry::subset(const Rectangle<R>& r, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(r,gms,"tribool subset(Rectangle r, GridMaskSet gms)");
  if(!subset(r,gms.bounding_box())) {
    return false;
  }
  return subset(outer_approximation(r,gms.grid()),gms);
}


template<class R>
tribool
Geometry::subset(const GridMaskSet<R>& gms, const Rectangle<R>& r)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(gms,r,"tribool subset(GridMaskSet gms, Rectangle r)");
  if(!subset(gms,r.bounding_box())) {
    return false;
  }
  return subset(gms,outer_approximation(r,gms.grid()));
}





template<class R>
Geometry::GridMaskSet<R>
Geometry::join(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
{
  ARIADNE_CHECK_SAME_GRID(gms1,gms2,"GridMaskSet join(GridMaskSet gms1, GridMaskSet gms2)");
  if(gms1.block()==gms2.block()) { 
    return GridMaskSet<R>(gms1.grid(), join(gms1.lattice_set(),gms2.lattice_set()));
  } else {
    GridMaskSet<R> result(gms1.grid(),Combinatoric::rectangular_hull(gms1.block(),gms2.block()));
    result.adjoin(gms1);
    result.adjoin(gms2);
    return result;
  }
}


template<class R>
Geometry::GridMaskSet<R>
Geometry::regular_intersection(const GridMaskSet<R>& gms, const GridBlock<R>& gb)
{
  ARIADNE_CHECK_SAME_GRID(gms,gb,"GridMaskSet regular_intersection(GridMaskSet gms, GridBlock gb)");
  return GridMaskSet<R>(gms.grid(), regular_intersection(gms.lattice_set(),gb.lattice_set()));
}

template<class R>
Geometry::GridMaskSet<R>
Geometry::regular_intersection(const GridBlock<R>& gb, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_SAME_GRID(gb,gms,"GridMaskSet regular_intersection(GridBlock gb, GridMaskSet gms)");
  return GridMaskSet<R>(gb.grid(), regular_intersection(gb.lattice_set(),gms.lattice_set()));
}

template<class R>
Geometry::GridMaskSet<R>
Geometry::regular_intersection(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
{
  ARIADNE_CHECK_SAME_GRID(gms1,gms2,"GridMaskSet regular_intersection(GridMaskSet gms1, GridMaskSet gms2)");
  if(gms1.block()==gms2.block()) { 
    return GridMaskSet<R>(gms1.grid(), regular_intersection(gms1.lattice_set(),gms2.lattice_set()));
  } else {
    GridMaskSet<R> result(gms1.grid(),Combinatoric::rectangular_hull(gms1.block(),gms2.block()));
    result.adjoin(gms1);
    result.restrict(gms2);
    return result;
  }
}

template<class R>
Geometry::GridCellListSet<R>
Geometry::regular_intersection(const GridCellListSet<R>& gcls, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_SAME_GRID(gcls,gms,"GridCellListSet regular_intersection(GridCellListSet gcls, GridMaskSet gms)");
  return GridCellListSet<R>(gcls.grid(), regular_intersection(gcls.lattice_set(),gms.lattice_set()));
}

template<class R>
Geometry::GridCellListSet<R>
Geometry::regular_intersection(const GridMaskSet<R>& gms, const GridCellListSet<R>& B)
{
  ARIADNE_CHECK_SAME_GRID(gms,B,"GridCellListSet regular_intersection(GridMaskSet<R>,GridCellListSet<R>)");
  return GridCellListSet<R>(gms.grid(), regular_intersection(gms.lattice_set(),B.lattice_set()));
}


template<class R>
Geometry::GridMaskSet<R>
Geometry::difference(const GridMaskSet<R>& gms1, const GridMaskSet<R>& gms2)
{
  ARIADNE_CHECK_SAME_GRID(gms1,gms2,"GridMaskSet difference(GridMaskSet<R>,GridMaskSet<R>)");
  if(gms1.block()==gms2.block()) { 
    return GridMaskSet<R>(gms1.grid(), difference(gms1.lattice_set(),gms2.lattice_set()));
  } else {
    GridMaskSet<R> result(gms1.grid(),Combinatoric::rectangular_hull(gms1.block(),gms2.block()));
    result.adjoin(gms1);
    result.remove(gms2);
    return result;
  }
}

template<class R>
Geometry::GridCellListSet<R>
Geometry::difference(const GridCellListSet<R>& gcls, const GridMaskSet<R>& gms)
{
  ARIADNE_CHECK_SAME_GRID(gcls,gms,"difference(GridCellListSet gcls, GridMaskSet gms)");
  return GridCellListSet<R>(gcls.grid(),Combinatoric::difference(gcls.lattice_set(),gms.lattice_set()));
}


template<class R>
Geometry::GridMaskSet<R>::operator ListSet< Rectangle<R> >() const
{
  ARIADNE_CHECK_BOUNDED(*this,"GridMaskSet::operator ListSet<Rectangle>()");
  ListSet< Rectangle<R> > result(this->dimension());
  for(typename GridMaskSet::const_iterator riter=begin(); riter!=end(); ++riter) {
    Rectangle<R> r(*riter);
    result.push_back(r);
  }
  return result;
}




template<class R, class BS>
Geometry::GridMaskSet<R>
Geometry::outer_approximation(const ListSet<BS>& ls, const FiniteGrid<R>& fg) 
{
  GridMaskSet<R> result(fg);
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls,fg,"outer_approximation(ListSet<BS> ls, FiniteGrid<R> fg)\n");
  
  for(typename ListSet<BS>::const_iterator bs_iter=ls.begin(); bs_iter!=ls.end(); ++bs_iter) {
    result.adjoin_outer_approximation(*bs_iter);
  }
  return result;
}


template<class R>
Geometry::GridMaskSet<R>
Geometry::outer_approximation(const SetInterface<R>& set, const FiniteGrid<R>& fg) 
{
  ARIADNE_LOG(4,"GridMaskSet outer_approximation(SetInterface, FiniteGrid)\n");
  GridMaskSet<R> result(fg);
  ARIADNE_CHECK_EQUAL_DIMENSIONS(set,fg,"outer_approximation(PartitionTreeSet<R>,FiniteGrid<R>)");
  
  const Grid<R>& g=fg.grid();
  const Combinatoric::LatticeBlock& lb=fg.lattice_block();
  
  for(typename Combinatoric::LatticeBlock::const_iterator iter=lb.begin(); iter!=lb.end(); ++iter) {
    GridCell<R> gc(g,*iter);
    if(!bool(set.disjoint(Rectangle<R>(gc)))) {
      result.adjoin(gc);
    }
  }
  return result;
}




template<class R>
Geometry::GridMaskSet<R>
Geometry::inner_approximation(const SetInterface<R>& s, const FiniteGrid<R>& fg) 
{
  ARIADNE_LOG(4,"GridMaskSet inner_approximation(SetInterface s, FiniteGrid fg)\n") 
    GridMaskSet<R> result(fg);
  GridBlock<R> gb(fg.grid(),fg.lattice_block());
  Rectangle<R> r;
  ARIADNE_LOG(5,"  testing cells in "<<gb<<"\n")
    for(typename GridBlock<R>::const_iterator gb_iter=gb.begin();
        gb_iter!=gb.end(); ++gb_iter)
      {
        r=*gb_iter;
        if(s.superset(r)) {
          result.adjoin(*gb_iter);
        }
        ARIADNE_LOG(6,"s.superset("<<r<<")="<<s.superset(r)<<"\n");
      }
  return result;
}





template<class R>  
Geometry::ListSet< Geometry::Rectangle<R> >
Geometry::lower_approximation(const SetInterface<R>& s, const FiniteGrid<R>& fg) 
{
  ARIADNE_LOG(4,"ListSet<Rectangle> lower_approximation(SetInterface s, FiniteGrid fg)\n") 
    ListSet< Rectangle<R> > result;
  GridBlock<R> gb(fg.grid(),fg.lattice_block());
  Rectangle<R> r;
  Rectangle<R> nr;
  ARIADNE_LOG(5,"  testing cells in "<<gb<<"\n")
    for(typename GridBlock<R>::const_iterator gb_iter=gb.begin();
        gb_iter!=gb.end(); ++gb_iter)
      {
        r=*gb_iter;
        if(!bool(s.disjoint(r))) {
          nr=gb_iter->neighbourhood();
          if(s.intersects(nr)) {
            result.adjoin(nr);
          }
        }
        ARIADNE_LOG(6,"s.interscts("<<r<<")="<<s.intersects(r)<<"\n");
      }
  return result;
}







template<class R>
std::ostream& 
Geometry::GridMaskSet<R>::write(std::ostream& os) const 
{
  os << "GridMaskSet( " << std::flush;
  os << " grid=" << this->grid() << ",";
  os << " block=" << this->block() << ",";
  os << " extent=" << Rectangle<R>(this->bounds()) << ",";
  os << " size=" << this->size() << ",";
  os << " capacity=" << this->capacity() << ",";
  //      os << "  mask=" << this->mask() << std::endl;
  os << " )";
  return os;
}


}
                                                    