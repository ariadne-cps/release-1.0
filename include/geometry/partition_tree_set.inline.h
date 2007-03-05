/***************************************************************************
 *            partition_tree_set.inline.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
namespace Ariadne {
  namespace Geometry {

          
    template<class R>
    class PartitionTreeIterator 
      : public boost::iterator_facade<PartitionTreeIterator<R>,
                                      PartitionTreeCell<R>,
                                      boost::forward_traversal_tag,
                                      PartitionTreeCell<R> >
    {
      friend class PartitionTree<R>;
     public:
      PartitionTreeIterator(const Rectangle<R>& bb, 
                            const Combinatoric::SubdivisionSequence& ss, 
                            Combinatoric::BinaryTree::const_iterator i)
        : _bounding_box(bb), _subdivisions(ss), _base(i) { }
     private:
      bool equal(const PartitionTreeIterator<R>& other) const {
        return this->_base == other._base; }
      void increment() { ++_base; }
      PartitionTreeCell<R> dereference() const { 
        return PartitionTreeCell<R>(_bounding_box,_subdivisions,*_base); }
     private:
      const Rectangle<R> _bounding_box;
      const Combinatoric::SubdivisionSequence _subdivisions;
      Combinatoric::BinaryTree::const_iterator _base;
    };

    template<class R>
    class PartitionTreeSetIterator 
      : public boost::iterator_facade<PartitionTreeSetIterator<R>,
                                      PartitionTreeCell<R>,
                                      boost::forward_traversal_tag,
                                      PartitionTreeCell<R> >
    {
      friend class PartitionTree<R>;
     public:
     private:
      bool equal(const PartitionTreeSetIterator<R>& other) const {
        return this->_base == other._base; }
      void increment() { ++_base; }
      PartitionTreeCell<R> dereference() const { 
        return PartitionTreeCell<R>(_bounding_box,_subdivisions,*_base); }
     private:
      const Rectangle<R> _bounding_box;
      const Combinatoric::SubdivisionSequence _subdivisions;
      Combinatoric::MaskedBinaryTree::const_iterator _base;
    };

    template<class R>
    PartitionScheme<R>::PartitionScheme(const Rectangle<R>& bb, const Combinatoric::SubdivisionSequence& sc)
      : _unit_box(bb), _subdivisions(sc)
    {
    }


    template<class R>
    bool 
    PartitionScheme<R>::operator==(const PartitionScheme<R>& pg) const 
    {
      return _unit_box==pg._unit_box && _subdivisions==pg._subdivisions;
    }

    template<class R> inline
    bool 
    PartitionScheme<R>::operator!=(const PartitionScheme<R>& pg) const 
    {
      return !(*this==pg); 
    }
      
    template<class R> inline
    const Rectangle<R>& 
    PartitionScheme<R>::unit_box() const 
    {
      return _unit_box; 
    }

    template<class R> inline
    const Combinatoric::SubdivisionSequence& 
    PartitionScheme<R>::subdivisions() const 
    {
      return _subdivisions; 
    }

    template<class R> inline
    dimension_type 
    PartitionScheme<R>::dimension() const 
    {
      return _subdivisions.dimension(); 
    }



    template<class R> inline
    PartitionTreeCell<R>::PartitionTreeCell(const Rectangle<R>& r, const Combinatoric::SubdivisionTreeCell& c)
      : _unit_box(r), _subdivision_cell(c)
    {
      check_equal_dimensions(r,c,__PRETTY_FUNCTION__); 
    }

    template<class R> inline
    PartitionTreeCell<R>::PartitionTreeCell(const Rectangle<R>& r, 
                                            const Combinatoric::SubdivisionSequence& s, 
                                            const Combinatoric::BinaryWord& w) 
      : _unit_box(r), _subdivision_cell(s,w) 
    { 
      //check_dimension_size(r,w,"PartitionTreeCell<R>::PartitionTreeCell(Rectangle<R>,SubdivisionSequence,BinaryWord"); 
    }

    template<class R> inline
    const Rectangle<R>& 
    PartitionTreeCell<R>::unit_box() const 
    {
      return this->_unit_box; 
    }

    template<class R> inline
    const Combinatoric::SubdivisionTreeCell& 
    PartitionTreeCell<R>::subdivision_cell() const 
    {
      return this->_subdivision_cell; 
    }

    template<class R> inline
    dimension_type 
    PartitionTreeCell<R>::dimension() const 
    {
      return this->_subdivision_cell.dimension(); 
    }

    template<class R> inline
    tribool 
    PartitionTreeCell<R>::empty() const 
    {
      return false; 
    }

    template<class R> inline
    tribool 
    PartitionTreeCell<R>::bounded() const 
    {
      return true; 
    }





    template<class R> inline
    PartitionTree<R>::PartitionTree(const Rectangle<R>& r, 
                                    const Combinatoric::SubdivisionSequence& s, 
                                    const Combinatoric::BinaryTree& t)
      : _unit_box(r), _subdivision_tree(s,t) 
    {
    }

    template<class R> inline
    PartitionTree<R>::PartitionTree(const PartitionScheme<R>& ps, const Combinatoric::BinaryTree& t)
        : _unit_box(ps.unit_box()), _subdivision_tree(ps.subdivisions(),t) 
    {
    }

    template<class R> inline
    const Rectangle<R>& 
    PartitionTree<R>::unit_box() const 
    {
      return _unit_box; 
    }

    template<class R> inline
    const Combinatoric::SubdivisionTree& 
    PartitionTree<R>::subdivision_tree() const 
    {
      return _subdivision_tree; 
    }

    template<class R> inline
    dimension_type 
    PartitionTree<R>::dimension() const 
    {
      return _subdivision_tree.dimension(); 
    }

    template<class R> inline
    const Combinatoric::SubdivisionSequence& 
    PartitionTree<R>::subdivisions() const 
    {
      return _subdivision_tree.subdivisions(); 
    }

    template<class R> inline
    const Combinatoric::BinaryTree& 
    PartitionTree<R>::binary_tree() const 
    {
      return _subdivision_tree.binary_tree(); 
    }

    template<class R> inline
    size_type 
    PartitionTree<R>::size() const 
    {
      return _subdivision_tree.size(); 
    }

    template<class R> inline
    PartitionScheme<R> 
    PartitionTree<R>::scheme() const 
    {
      return  PartitionScheme<R>(unit_box(),subdivisions()); 
    }
      
    template<class R> inline
    typename PartitionTree<R>::const_iterator 
    PartitionTree<R>::begin() const 
    {
      return const_iterator(_unit_box,_subdivision_tree.begin()); 
    }
    template<class R> inline
    typename PartitionTree<R>::const_iterator 
    PartitionTree<R>::end() const 
    {
      return const_iterator(_unit_box,_subdivision_tree.end()); 
    }





    template<class R> inline
    PartitionTreeSet<R>::PartitionTreeSet(const PartitionScheme<R>& g)
      : _unit_box(g.unit_box()), _subdivision_set(g.subdivisions()) 
    {
    }

    template<class R> inline
    PartitionTreeSet<R>::PartitionTreeSet(const PartitionScheme<R>& g, const Combinatoric::BinaryTree& t, const BooleanArray& m)
      : _unit_box(g.unit_box()), _subdivision_set(g.subdivisions(),t,m) 
    {
    }

    template<class R> inline
    PartitionTreeSet<R>::PartitionTreeSet(const PartitionTree<R>& t, const BooleanArray& m)
      : _unit_box(t.unit_box()), _subdivision_set(t.subdivisions(),t.binary_tree(),m)
    {
    }

    template<class R> inline
    PartitionTreeSet<R>::PartitionTreeSet(const Rectangle<R>& r, 
                                          const Combinatoric::SubdivisionSequence& s, 
                                          const Combinatoric::BinaryTree& t, 
                                          const BooleanArray& m)
      : _unit_box(r), _subdivision_set(s,t,m)
    {
    }


    template<class R> inline
    Rectangle<R> 
    PartitionTreeSet<R>::bounding_box() const 
    {
      return _unit_box; 
    }

    template<class R> inline
    const Rectangle<R>& 
    PartitionTreeSet<R>::unit_box() const 
    {
      return _unit_box; 
    }

    template<class R> inline
    const Combinatoric::SubdivisionTreeSet& 
    PartitionTreeSet<R>::subdivision_set() const 
    {
      return _subdivision_set; 
    }

    template<class R> inline
    dimension_type 
    PartitionTreeSet<R>::dimension() const 
    {
      return _subdivision_set.dimension(); 
    }

    template<class R> inline
    const Combinatoric::SubdivisionSequence& 
    PartitionTreeSet<R>::subdivisions() const 
    {
      return _subdivision_set.subdivisions(); 
    }

    template<class R> inline
    const Combinatoric::BinaryTree& 
    PartitionTreeSet<R>::binary_tree() const 
    {
      return _subdivision_set.binary_tree(); 
    }
      
    template<class R> inline
    const BooleanArray& 
    PartitionTreeSet<R>::mask() const 
    {
      return _subdivision_set.mask(); 
    }

    template<class R> inline
    size_type 
    PartitionTreeSet<R>::capacity() const 
    {
      return _subdivision_set.capacity(); 
    }

    template<class R> inline
    size_type 
    PartitionTreeSet<R>::size() const 
    {
      return _subdivision_set.size(); 
    }

    template<class R> inline
    SizeArray 
    PartitionTreeSet<R>::depths() const 
    {
      return _subdivision_set.depths(); 
    }

    template<class R> inline
    size_type 
    PartitionTreeSet<R>::depth() const 
    {
      return _subdivision_set.depth(); 
    }
      
    template<class R> inline
    PartitionScheme<R> 
    PartitionTreeSet<R>::scheme() const 
    {
      return PartitionScheme<R>(bounding_box(),subdivisions()); 
    }
      
    template<class R> inline
    PartitionTree<R> 
    PartitionTreeSet<R>::partition_tree() const 
    {
      return PartitionTree<R>(bounding_box(),subdivisions(),binary_tree()); 
    }
        
    template<class R> inline
    typename PartitionTreeSet<R>::const_iterator 
    PartitionTreeSet<R>::begin() const 
    {
      return const_iterator(_unit_box,_subdivision_set.begin()); 
    }

    template<class R> inline
    typename PartitionTreeSet<R>::const_iterator 
    PartitionTreeSet<R>::end() const 
    {
      return const_iterator(_unit_box,_subdivision_set.end()); 
    }
      
    template<class R> inline
    tribool 
    PartitionTreeSet<R>::empty() const 
    {
      return this->size()==0; 
    }

    template<class R> inline
    tribool 
    PartitionTreeSet<R>::bounded() const 
    {
      return true; 
    }
   





    template<class R> inline
    std::ostream& 
    operator<<(std::ostream& os, const PartitionScheme<R>& ps) 
    { 
      return ps.write(os);
    }
    
    template<class R> inline
    std::ostream& 
    operator<<(std::ostream& os, const PartitionTreeCell<R>& ptc)
    { 
      return ptc.write(os);
    }
    
    template<class R> inline
    std::ostream& 
    operator<<(std::ostream& os, const PartitionTree<R>& pt) 
    { 
      return pt.write(os);
    }
    
    template<class R> inline
    std::ostream& 
    operator<<(std::ostream& os, const PartitionTreeSet<R>& pts)
    { 
      return pts.write(os);
    }
    
    
  }
}