/***************************************************************************
 *            bdd_set.cc
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

#include <bdd.h>

#include "bdd_set.h"

namespace Ariadne {

typedef size_t size_type;

/********************************** Utility functions **********************************/

// Test if BuDDy has been initialized. If not, initialize it.
void _initialize_bddlib() {
    if(!bdd_isrunning()) {
        // Set the initial nodenum and cachesize of a medium sized load.
        // TO DO: make this a parameter for the user.
        bdd_init(100000, 10000);
    }
}

// Compute the minimum cell depth so that the root cell covers a given box
uint _minimum_primary_cell_depth(const Grid& grid, const Box& box) {
    // Grid and box must have the same dimension
    ARIADNE_ASSERT(grid.dimension() == box.dimension());
    // raise an error if dimension is zero
    ARIADNE_ASSERT(grid.dimension() != 0);
    // start with cell size equal to the primary cell and depth 0
    Box cell = grid.primary_cell();
    uint depth = 0;
    uint dim = grid.dimension();
    uint i = 0;
    uint p = 0;
    // increase depth until cell covers box
    while(!definitely(cell.covers(box))) {
        depth++;
        // compute the dimension to increase
        i = (i + dim - 1) % dim;
        // determine if this is an odd or even instance of the increase
        p = (depth / dim) % 2;
        if(p == 0) {    // positive instance: extend to the right
            cell[i] = cell[i] + Interval(0.0,cell[i].width());
        } else {        // negative instance: extend to the left
            cell[i] = cell[i] - Interval(0.0,cell[i].width());
        }
    }
    return i;
}

// Compute the size of the root cell from the grid and the depth of the primary cell
Box _compute_root_cell(const Grid& grid, uint depth) {
    // raise an error if dimension is zero
    ARIADNE_ASSERT(grid.dimension() != 0);
    // start with the primary cell
    Box cell = grid.primary_cell();
    uint dim = grid.dimension();
    uint i = 0;
    uint p = 0;
    // enlarge the cell until depth is reached
    for(uint d = 0; d < depth; ++d) {
        // compute the dimension to increase
        i = (i + dim - 1) % dim;
        // determine if this is an odd or even instance of the increase
        p = (depth / dim) % 2;
        if(p == 0) {    // positive instance: extend to the right
            cell[i] = cell[i] + Interval(0.0,cell[i].width());
        } else {        // negative instance: extend to the left
            cell[i] = cell[i] - Interval(0.0,cell[i].width());
        }
    }
    
    return cell;
}

/************************************* BDDTreeSet **************************************/

BDDTreeSet::BDDTreeSet( )
    : _grid(0)
    , _primary_cell_depth(0)
{
    _initialize_bddlib();
    _bdd = bdd_false();
}

BDDTreeSet::BDDTreeSet( const BDDTreeSet & set )
    : _grid(set.grid())
    , _primary_cell_depth(set.primary_cell_depth())
    , _bdd(set._bdd)
{
}

BDDTreeSet::BDDTreeSet( const uint dimension, const bool enable )
    : _grid(dimension)
    , _primary_cell_depth(0)
{
    _initialize_bddlib();
    if(enable) {
        _bdd = bdd_true();
    } else {
        _bdd = bdd_false();
    }
}

BDDTreeSet::BDDTreeSet( const Grid& grid, const bool enable )
    : _grid(grid)
    , _primary_cell_depth(0)
{
    _initialize_bddlib();
    if(enable) {
        _bdd = bdd_true();
    } else {
        _bdd = bdd_false();
    }    
}

BDDTreeSet::BDDTreeSet( const Grid& grid, const Box & domain )
    : _grid(grid)
{
    _initialize_bddlib();
    _bdd = bdd_false();     // The set is initially empty
    // set the primary cell depth so that the root cell covers the domain    
    _primary_cell_depth = _minimum_primary_cell_depth(grid, domain);
}

BDDTreeSet& BDDTreeSet::operator=( const BDDTreeSet & set ) {
    //Dereference the old bdd and substitute with the new one
    if(this != &set) {
        this->_grid = set.grid();
        this->_primary_cell_depth = set.primary_cell_depth();
        this->_bdd = set._bdd;
    }
    return *this;
}

BDDTreeSet* BDDTreeSet::clone() const {
    return new BDDTreeSet( *this );
}

BDDTreeSet::~BDDTreeSet() {
    
}

bool BDDTreeSet::empty() const {
    // A zero dimension set is always empty
    if(this->dimension() == 0) return true;
    // otherwise, the set is empty iff the BDD is the constant false
    return (this->_bdd == bddfalse);
}

size_t BDDTreeSet::size() const {
    ARIADNE_NOT_IMPLEMENTED;
}

uint BDDTreeSet::dimension() const {
    return this->_grid.dimension();
}

const Grid& BDDTreeSet::grid() const {
    return this->_grid;
}

uint BDDTreeSet::primary_cell_depth() const {
    return this->_primary_cell_depth;
}

double BDDTreeSet::measure() const {
    ARIADNE_NOT_IMPLEMENTED;    
}

Box BDDTreeSet::primary_cell() const {
    return this->_grid.primary_cell();
}

Box BDDTreeSet::bounding_box() const {
    return _compute_root_cell(this->grid(), this->primary_cell_depth());
}

/*
    bool operator==(const BDDTreeSet& anotherBDDTreeSet) const;

    bool subset( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    bool superset( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    bool disjoint( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    bool overlap( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    tribool subset( const Box& box ) const;

    tribool superset( const Box& box ) const;

    tribool disjoint( const Box& box  ) const;

    tribool overlaps( const Box& box ) const;

    void clear( );

    BDDTreeSet join( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    BDDTreeSet intersection( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    BDDTreeSet difference( const BDDTreeSet& set1, const BDDTreeSet& set2 );

    void adjoin( const BDDTreeSet& set );

    void restrict( const BDDTreeSet& set );

    void remove( const BDDTreeSet& set );

    void restrict_to_height( const uint height );

    void adjoin_over_approximation( const Box& box, const uint subdiv );

    void adjoin_outer_approximation( const CompactSetInterface& set, const uint subdiv );

    void adjoin_lower_approximation( const OvertSetInterface& set, const uint height, const uint subdiv );

    void adjoin_lower_approximation( const OvertSetInterface& set, const Box& bounding_box, const uint subdiv );

    void adjoin_lower_approximation( const LocatedSetInterface& set, const uint subdiv );

    void adjoin_inner_approximation( const OpenSetInterface& set, const uint height, const uint subdiv );

    void adjoin_inner_approximation( const OpenSetInterface& set, const Box& bounding_box, const uint subdiv );

    const_iterator begin() const;

    const_iterator end() const;

    BDDTreeSet& operator=( const BDDTreeSet &otherSubset);

    operator ListSet<Box>() const;

    void draw(CanvasInterface& canvas) const;

    std::ostream& write(std::ostream& os) const;

    void import_from_file(const char*& filename);

    void export_to_file(const char*& filename);

*/

} // namespace Ariadne

