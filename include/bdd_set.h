/***************************************************************************
 *            bdd_set.h
 *
 *  Copyright  2011  Davide Bresolin
 *            davide.bresolin@univr.it
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file bdd_set.h
 *  \brief Implicit representation of %BDDTreeSets using Binary Decision Diagrams.
 */

#ifndef ARIADNE_BDD_SET_H
#define ARIADNE_BDD_SET_H

#include <iostream>
#include <string>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/shared_ptr.hpp>

#include <bdd.h>

#include "tribool.h"

#include "box.h"
#include "list_set.h"

#include "numeric.h"

#include "set_interface.h"
#include "taylor_set.h"
#include "function_set.h"

#include "grid.h"
#include "graphics_interface.h"

using namespace std;
using namespace Ariadne;

namespace Ariadne {

/*Some type definitions*/
typedef std::vector<bool> BooleanArray;
typedef array<int> IndexArray;
typedef array<unsigned int> SizeArray;

typedef unsigned short dimension_type;

/*Some pre-declarations*/
class Grid;
// class BDDCell;
class BDDTreeSet;
class BDDTreeConstIterator;
struct PathElement;

/*Declarations of classes in other files*/
template<class BS> class ListSet;

// std::ostream& operator<<(std::ostream& os, const BDDCell& cell);
// std::ostream& operator<<(std::ostream& os, const BDDTreeCursor& treecursor);
std::ostream& operator<<(std::ostream& os, const BDDTreeSet& set);
std::ostream& operator<<(std::ostream& os, const PathElement& pe);
// 
// // bool subset(const BDDCell& cell, const BDDTreeSet& set);
// // bool overlap(const BDDCell& cell, const BDDTreeSet& set);
bool subset(const BDDTreeSet& set1, const BDDTreeSet& set2);
bool superset(const BDDTreeSet& set1, const BDDTreeSet& set2);
bool disjoint(const BDDTreeSet& set1, const BDDTreeSet& set2);
bool overlap(const BDDTreeSet& set1, const BDDTreeSet& set2);
bool restricts(const BDDTreeSet& set1, const BDDTreeSet& set2);
// 
BDDTreeSet join(const BDDTreeSet& set1, const BDDTreeSet& set2);
BDDTreeSet intersection(const BDDTreeSet& set1, const BDDTreeSet& set2);
BDDTreeSet difference(const BDDTreeSet& set1, const BDDTreeSet& set2);
// 
// BDDTreeSet outer_approximation(const Box& box, const Grid& grid, const uint subdiv);
// BDDTreeSet outer_approximation(const CompactSetInterface& set, const Grid& grid, const uint subdiv);
// BDDTreeSet outer_approximation(const CompactSetInterface& set, const uint subdiv);
// template<class BS> BDDTreeSet outer_approximation(const ListSet<BS>& set, const uint subdiv);
// BDDTreeSet inner_approximation(const OpenSetInterface& set, const Grid& grid, const uint height, const uint subdiv);

// TO DO: Implement serialization
// template<class A> void serialize(A& archive, const BDDTreeSet& set, const uint version);

/*! \brief The BDDTreeSet class that represents a set of cells on a variable size grid.
 * The cells can be enabled or disabled (on/off), indicating whether they belong to the set or not.
 *
 * A Binary Decision Diagram (BDDs) is used to represent the set in a symbolic way. Every cell 
 * is identified by a path from the root of the BDD. Going down one level in the BDD correspond
 * to splitting the cell into two subcells on one dimension. The false branch is the subcell with
 * smaller values, while the true branch is the subcell with greater values.
 *
 * The size of the cells is determined by the PrimaryCell, that is not necessarily the root of the
 * BDD. In general, is an internal cell of the BDD.
 * 
 * TO DO: improve description.
 */
class BDDTreeSet : public DrawableInterface {
  protected:
    // The underlying grid on wicht the cells are built
    Grid _grid;
    // The height of the root cell, that is, the number of subdivision needed to obtain a cell of the base grid
    uint _root_cell_height;             
    // The coordinates of the root cell
    array<int> _root_cell_coordinates;  
    // The bdd encoding the enabled cells
    bdd _bdd;
    // The depth to which the cells must me minced
    int _mince_depth;
    
  public:
    /*! \brief A short name for the constant iterator */
    typedef BDDTreeConstIterator const_iterator;

    //@{
    //! \name Constructors

    /*! \brief Create a %BDDTreeSet based on zero dimensions.
     *  This constructor is needed to use the Boost Serialization library.
     */
    BDDTreeSet( );

    /*! \brief The copy constructor.
     */
    BDDTreeSet( const BDDTreeSet & set );

    /*! Create a %BDDTreeSet based on \a grid with the given \a root_cell_height and 
     *  \a root_cell_coordinates. Activate the cells defined by \a enabled_cells.
     */
    BDDTreeSet( const Grid& grid, const uint root_cell_height, 
                const array<int>& root_cell_coordinates, const bdd& enabled_cells);

    /*! A simple constructor that creates the [0, 1]*...*[0, 1] cell in the
     *  \a dimension - dimensional space, with root cell height = 0. 
     * If enable == true then the cell is enabled.
     */
    explicit BDDTreeSet( const uint dimension, const bool enable = false );

    /*! \brief Construct a set with cells based on \a grid. 
     * If enable == true then the primary cell is enabled.
     */
    explicit BDDTreeSet( const Grid& grid, const bool enable = false  );

    //@}

    //@{
    //! \name Cloning/Copying/Assignment

    /*! \brief Return a copy of the %BDDTreeSet.
     */
    BDDTreeSet* clone() const;

    //@}

    //@{
    //! \name Properties

    /*! \brief True if the set is empty. */
    bool empty() const;

    /*! \brief The number of activated cells in the set. */
    size_t size() const;

    /*! \brief The dimension of the set. */
    uint dimension() const;

    /*! \brief Returns a constant reference to the underlying grid. */
    const Grid& grid() const;
    
    /*! \brief The height of the root cell, that is, the number of subdivision needed to obtain a cell of the base grid.
     */
    uint root_cell_height() const;
    
    /*! \brief The depth of the set, that is, the longest path from a cell of the base grid to a leaf.
     */
    uint depth() const;

    /*! \brief The coordinate of the lower corner of the root cell
     */
    array<int> root_cell_coordinates() const;
    
    /*! \brief The depth to wich all cells are minced (returns -1 if the cells are not minced).
     */
    int mince_depth() const;
    
    /*! \brief The BDD representing the enabled cells.
     */
    const bdd& enabled_cells() const;

    /*! The measure (area, volume) of the set in Euclidean space. */
    double measure() const;

    /*! \brief Returns the %Box corresponding to the root cell of this \a BDDTreeSet
     */
    Box root_cell() const;

    /*! \brief Computes a bounding box for a %BDDTreeSet. */
    Box bounding_box() const;

    /*! \brief Allows to test if two BDDTreeSet are equal. The method returns true if
     * the grids are equal, the root cells have the same coordinate and height, and the BBDs are equal.
     */
    bool operator==(const BDDTreeSet& anotherBDDTreeSet) const;

    bool operator!=(const BDDTreeSet& anotherBDDTreeSet) const;

    //@}

    //@{
    //! \name Geometric Predicates
    
    /*! \brief Tests if a %BDDTreeSet \a set1 is a subset of \a set2. */
    friend bool subset( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    /*! \brief Tests if a %BDDTreeSet \a set1 is a superset of \a set2. */
    friend bool superset( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    /*! \brief Tests if two %BDDTreeSets are disjoint.
     */
    friend bool disjoint( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    /*! \brief Tests if two %BDDTreeSets overlap.
     */
    friend bool overlap( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    /*! \brief Tests if \a set1 restricts \a set2.
     */
    friend bool restricts( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    /*! \brief Tests if a %BDDTreeSet is a subset of a box. */
    tribool subset( const Box& box ) const ;

    /*! \brief Tests if a %BDDTreeSet is a superset of a box. */
    tribool superset( const Box& box ) const;

    /*! \brief Tests if (the closure of) a %BDDTreeSet is disjoint from a box. */
    tribool disjoint( const Box& box  ) const;

    /*! \brief Tests if a %BDDTreeSet overlaps a box. */
    tribool overlaps( const Box& box ) const;

    //@}
    
    //@{
    //! \name Subdivisions

    /*! \brief Subdivides the cells in such a way every cell is equal or smaller to the base cell 
     *  subdivided \a subdiv times in each dimension.
     */
    void mince( const uint subdiv );

    /*! \brief Recombines the subdivisions, for instance if all subcells of a cell are
     * enabled/disabled then they are put together.
     */
    void recombine();

    //@}

    //@{
    //! \name Geometric Operations

    /*! \brief Clears the set (makes empty set on same grid). */
    void clear( );

    /*! \brief Minimize the height of the set so that the root cell is the minimal cell covering the set.
     *  Returns the new root_cell_height.
     */
    int minimize_height();

    /*! \brief Increase the root cell height to \a new_height, and returns the new root_cell_height.
     *  If \a new_height is smaller than the current \a root_cell_height the set is not modified.
     */
    int increase_height(uint new_height);

    /*! \brief Increase the root cell until \a box is a subset of the root cell, and returns the new root_cell_height.
     */
    int increase_height(const Box& box);

    /*! \brief Join (make union of) two %BDDTreeSet. */
    friend BDDTreeSet join( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    /*! \brief The intersection of two %BDDTreeSet. 
     */
    friend BDDTreeSet intersection( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    /*! \brief The difference of two %BDDTreeSet. (Results in set1 minus set2) */
    friend BDDTreeSet difference( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    /*! \brief Adjoin (make inplace union with) another %BDDTreeSet. */
    void adjoin( const BDDTreeSet& set );

    /*! \brief Restrict to (make inplace intersection with) another %BDDTreeSet. */
    void restrict( const BDDTreeSet& set );

    /*! \brief Remove cells in another %BDDTreeSet. */
    void remove( const BDDTreeSet& set );

    //@}

    //@{
    //! \name Geometric Approximation

    /*! \brief Adjoin an over approximation of a box, computing to the given depth:
     *   \a subdiv -- defines, how many subdivisions in each dimension from the level of
     *   the primary cell we should make to get the proper cells for outer approximating \a set.
     *   \pre The box must have nonempty interior.
     */
    void adjoin_over_approximation( const Box& box, const uint subdiv );

    /*! \brief Adjoin an outer approximation of a given set, computing to the given depth.
     */
    void adjoin_outer_approximation( const CompactSetInterface& set, const uint subdiv );

    /*! \brief Adjoin a lower approximation to a given set, computing to the given height and depth:
     *   \a subdiv defines how many subdivisions in each dimension from the level of the
     *   base grid we should make to get the proper cells for outer approximating \a set.
     *   \a height defines how many merging in each dimension of the base grid are allowed to
     *   enlarge the root cell.
     *   A lower approximation comprises all cells intersecting a given set.
     */
    void adjoin_lower_approximation( const OvertSetInterface& set, const uint height, const uint subdiv );

    /*! \brief Adjoin a lower approximation to a given set restricted to the given bounding box,
     *   computing to the given depth: \a subdiv -- defines, how many subdivisions in each
     *   dimension from the level of the zero cell we should make to get the proper cells for outer
     *   approximating \a set. A lower approximation comprises all cells intersecting a given set.
     */
    void adjoin_lower_approximation( const OvertSetInterface& set, const Box& bounding_box, const uint subdiv );

    /*! \brief Adjoin a lower approximation to a given set, computing to the given depth
     *   \a subdiv -- defines, how many subdivisions in each dimension from the level of the
     *   zero cell we should make to get the proper cells for outer approximating \a set.
     *   A lower approximation comprises all cells intersecting a given set.
     */
    void adjoin_lower_approximation( const LocatedSetInterface& set, const uint subdiv );

    /*! \brief Adjoin an inner approximation to a given set, computing to the given height and depth:
     *   \a subdiv -- defines, how many subdivisions in each dimension from the level of the
     *   zero cell we should make to get the proper cells for outer approximating \a set.
     *   An inner approximation comprises all cells that are sub-cells of the given set.
     */
    void adjoin_inner_approximation( const OpenSetInterface& set, const uint height, const uint subdiv );

    /*! \brief Adjoin an inner approximation to a given set restricted to the given bounding box,
     *   computing to the given depth: \a subdiv -- defines, how many subdivisions in each
     *   dimension from the level of the zero cell we should make to get the proper cells for outer
     *   approximating \a set. An inner approximation comprises all cells that are sub-cells of
     *   the given set.
     */
    void adjoin_inner_approximation( const OpenSetInterface& set, const Box& bounding_box, const uint subdiv );

    /*! \brief Restrict to the cells that possibly overlaps with \a set. 
     *  The result is an outer approximation of the intersection with \a set.
     */
    void outer_restrict( const OvertSetInterface& set );

    /*! \brief Restrict to the cells that are definitely inside \a set.
     *  The result is an inner approximation of the intersection with \a set.
     */
    void inner_restrict( const OpenSetInterface& set );

    //@}

    //@{
    //! \name Iterators

    /*! \brief A constant iterator through the enabled leaf nodes of the subpaving. */
    const_iterator begin() const;

    /*! \brief A constant iterator to the end of the enabled leaf nodes of the subpaving. */
    const_iterator end() const;

    //@}

    //@{
    //! \name Conversions

    /*! \brief Convert to a list of ordinary boxes, unrelated to the grid. */
    operator ListSet<Box>() const;

    //@}

    //@}
    //! \name Input/Output

    /*! \brief Draw on a two-dimensional canvas. */
    void draw(CanvasInterface& canvas) const;

    /*! \brief Write to an output stream. */
    std::ostream& write(std::ostream& os) const;

    /*! \brief Import the content from the file \a filename.
     */
//    void import_from_file(const char*& filename);

    /*! \brief Export the tree to the file \a filename (without appending).
	 */
//    void export_to_file(const char*& filename);

    //@}
    
};

// definition of the path element for BDDTreeConstIterator
enum PEStatus {PE_NEW, PE_LEFT, PE_RIGHT};

struct PathElement {
    bdd obdd;
    PEStatus status;
    Box cell;
    uint split_coordinate;
    int root_var;
};

/*! \brief This class allows to iterate through the enabled leaf nodes of BDDTreeSet.
 * The return objects for this iterator are constant Boxes.
 */
class BDDTreeConstIterator 
    : public boost::iterator_facade< 
          BDDTreeConstIterator
        , Box const
        , boost::forward_traversal_tag 
      > 
{
  private:    
    /*! \brief The path from the root of the current cell */
    std::vector< PathElement > _path;
    // The depth to which the cells must me minced
    int _mince_depth;
        
    friend class boost::iterator_core_access;

    //@{
    //! \name Iterator Specific

    void increment();

    /*! \brief compare two iterators
     */
    bool equal( BDDTreeConstIterator const & other ) const;

    Box const& dereference() const;

    //@}

  public:
    //@{
    //! \name Constructors

    /*! \brief Default constructor constructs an invalid iterator.
     */
    BDDTreeConstIterator();

    /*! \brief The constructor that accepts the %BDDTreeSet \a set to iterate on
     */
    explicit BDDTreeConstIterator( const BDDTreeSet& set );

    /*! \brief The copy constructor.
     */
    BDDTreeConstIterator( const BDDTreeConstIterator& iter );

    //@}
};


template<class A> void serialize(A& archive, Ariadne::BDDTreeSet& set, const unsigned int version) {
    ARIADNE_NOT_IMPLEMENTED;
}


//! \brief Whether \a cons_set is disjoint from \a bdd_set.
tribool disjoint(const ConstraintSet& cons_set, const BDDTreeSet& bdd_set);
//! \brief Whether \a cons_set overlaps with \a bdd_set.
tribool overlaps(const ConstraintSet& cons_set, const BDDTreeSet& bdd_set);
//! \brief Whether \a cons_set covers \a bdd_set.
tribool covers(const ConstraintSet& cons_set, const BDDTreeSet& bdd_set);

//! \brief Evaluates \a bdd_set on \a cons_set in order to obtain (a superset of) the overlapping subset.
BDDTreeSet possibly_overlapping_subset(const BDDTreeSet& bdd_set, const ConstraintSet& cons_set);
//! \brief Applies \a cons_set to \a bdd_set in order to obtain the definitely covered subset.
BDDTreeSet definitely_covered_subset(const BDDTreeSet& bdd_set, const ConstraintSet& cons_set);

//! \brief Evaluates the codomain of \a func applied on the cells of \a bdd_set, each widened by \a eps.
//Box eps_codomain(const BDDTreeSet& bdd_set, const Vector<Float> eps, const VectorFunction& func);

//! \brief Projects \a bdd_set using the given \a indices.
BDDTreeSet project_down(const BDDTreeSet& bdd_set, const Vector<uint>& indices);

//! \brief Check whether \a covering_set covers \a covered_set with a tolerance of \a eps.
//! \details Since the cell boxes of \a covered_set, enlarged of \a eps, are checked against \a covering_set,
//! the two sets can feature different grids.
// tribool covers(const BDDTreeSet& covering_set, const BDDTreeSet& covered_set, const Vector<Float>& eps);

//! \brief Check whether \a covering_set covers \a covered_set with a tolerance of \a eps.
//! \details Since the cell boxes of \a covered_set are checked against an overapproximation (using \a accuracy) of the
//! epsilon-enlargement of \a covering_set, the two sets can feature different grids.
// tribool inside(const BDDTreeSet& covered_set, const BDDTreeSet& covering_set, const Vector<Float>& eps, int accuracy);

} // namespace Ariadne

#endif /* ARIADNE_BDD_SET_H */

