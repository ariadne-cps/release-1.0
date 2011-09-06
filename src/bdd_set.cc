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

/********************************** Utility functions **********************************/

// Test if BuDDy has been initialized. If not, initialize it.
void _initialize_bddlib() {
    // TO DO: make the different parameters definable by the user.
    if(!bdd_isrunning()) {
        // Set the initial nodenum and cachesize of a medium sized load.
        ARIADNE_ASSERT_MSG(bdd_init(100000, 10000) == 0, "Error in the initialization of the bdd library.");
        // Set the initial number of variables.
        ARIADNE_ASSERT_MSG(bdd_setvarnum(32) == 0, "Error in the initialization of the bdd library.");
    }
}

// Compute the grid corresponding to a given root height, starting from the base grid
Grid _compute_root_grid(const Grid& grid, uint height) {
    // std::cout << "_compute_root_grid(" << grid << ", " << height << ")" << std::endl;
    // raise an error if dimension is zero
    ARIADNE_ASSERT(grid.dimension() != 0);
    // Increase the height step by step;
    Grid res(grid);
    uint dim = res.dimension();
    for(uint h = 0; h != height; ++h) {
        // determine which dimension to merge 
        uint i = (dim - 1) - (h % dim);
        // if the occurrence of the merge is even, shift the origin
        if((h / dim) % 2 == 1) {
            res.set_origin_coordinate(i, res.origin()[i]-res.lengths()[i]);
        }
        res.set_length(i, 2.0*res.lengths()[i]);        
        // std::cout << "h = " << h << ", grid = " << res << std::endl;
    }
    
    return res;
}

// Compute the root cell from the grid, the height and the coordinates of the root cell
Box _compute_root_cell(const Grid& grid, uint height, const array<int>& coordinates) {
    // std::cout << "_compute_root_cell(" << grid << ", " << height << ", " << coordinates << ")" << std::endl;
    // raise an error if dimension is zero
    ARIADNE_ASSERT(grid.dimension() != 0);
    // compute the root grid
    Grid root_grid = _compute_root_grid(grid, height);
    // std::cout << "root grid: " << root_grid << std::endl;
    // compute the root cell from the coordinates and return it
    return root_grid.cell(coordinates);
}

// Compute recursively the number of enabled cells
size_t _enabled_cells_number(const bdd& b, int level) {
    // std::cout << "_enabled_cells_number(" << b.id() << "," << level << ")" << std::endl;
    // Base cases: the size of the constant true/false is 1/0
    if(b == bddtrue) return 1;
    if(b == bddfalse) return 0;
    
    // Inductive cases
    int rootlevel = bdd_var2level(bdd_var(b));      // get the level of the root
    // std::cout << "rootlevel = " << rootlevel << std::endl;
    if(rootlevel == level) {                        // no jumps of variables
        bdd left = bdd_low(b);                      // get the left branch
        bdd right = bdd_high(b);                    // get the right branch
        return _enabled_cells_number(left, level+1) +
               _enabled_cells_number(right, level+1);
    }
    if(rootlevel > level) {     // the root variable is greater than the current level
        // hence, the root variable has been deleted by the bdd construction
        // because it pointed to two identical childs.
        return 2*_enabled_cells_number(b, level+1);
    }
    // this branch should never be reached
    ARIADNE_FAIL_MSG("rootlevel cannot be greater than current level.");
    return 0;        
}

// Minimize the height of the bdd b, given a root_cell_height and a root_cell_coordinates.
// root_var is the index of the variable corresponding to the current root level
// This to allow representations where the root level is not the variable 0, so that
// only one shift of variable (of index n) at the end is needed.
int _minimize_height(bdd& b, uint& root_cell_height, array<int>& root_cell_coordinates, int root_var) 
{
    // if the height is zero do not decrease height
    if(root_cell_height == 0) return -root_var;
    // if the b is TRUE or FALSE do not decrease height
    if(b == bddtrue || b == bddfalse) return -root_var;
    // if the variable labelling the root is different from root_var do not decrease height
    if(bdd_var(b) != root_var) return -root_var;
    // get the two childs of b
    bdd left = bdd_low(b);
    bdd right = bdd_high(b);
    // if both children are different from FALSE, do not decrease height
    if(left != bddfalse && right != bddfalse) return -root_var;
    // exactly one of the two children is FALSE, decrease height by 1 and repeat
    root_cell_height--;
    // determine which dimension to split
    uint dim = root_cell_coordinates.size();
    uint i = (dim - 1) - (root_cell_height % dim);
    // determine if this is an odd or even occurrence of the split
    uint p = (root_cell_height / dim) % 2;
    // compute the new coordinate for the i-th dimension,
    // then repeat recursively on the correct child
    if(right == bddfalse) {    // the right child is false, consider the left one
        root_cell_coordinates[i] = 2*root_cell_coordinates[i] - p;
        b = left;
        return _minimize_height(b, root_cell_height, root_cell_coordinates, root_var+1);
    } else {    // the left child is false, consider the right one
        root_cell_coordinates[i] = 2*root_cell_coordinates[i] - p + 1;
        b = right;
        return _minimize_height(b, root_cell_height, root_cell_coordinates, root_var+1);
    }
}

// Increase the height of the bdd b to new_height, given a root_cell_height and a root_cell_coordinates.
// Assumes that the variables have been shifted to the correct indexes prior execution.
// WARNING: incorrect results are given if the variables are not shifted.
int _increase_height(bdd& b, uint& root_cell_height, array<int>& root_cell_coordinates, uint new_height) 
{
    // std::cout << "_increase_height(" << b << ", " << root_cell_height << ", " << root_cell_coordinates
    //          << ", " << new_height << ")" << std::endl;
    // if the new height is equal or smaller than the current one do nothing
    if(root_cell_height >= new_height) return root_cell_height;

    // Increase the root cell height by one and recursively call _minimize_height
    // compute the new root_cell_coordinate
    // determine which dimension to merge 
    uint dim = root_cell_coordinates.size();
    uint i = (dim - 1) - (root_cell_height % dim);
    // If the current height is an "odd" one, the origin of the grid is shifted
    // shift the i-th coordinate accordingly
    if((root_cell_height / dim) % 2 == 1) {
        root_cell_coordinates[i]++;
    }
    // Increase the root_cell_height by 1
    root_cell_height++;
    // the root_var is the difference between new_height and the old one
    int root_var = new_height - root_cell_height;
    // Test if the current bdd will be the left or right child of the new one
    if((root_cell_coordinates[i] % 2) == 0) {   // left child
        b = bdd_nithvar(root_var) & b;
    } else {    // right child
        b = bdd_ithvar(root_var) & b;
    }
    // half the i-th coordinate
    root_cell_coordinates[i] = root_cell_coordinates[i] / 2;
    
    // recursive call
    return _increase_height(b, root_cell_height, root_cell_coordinates, new_height);
}
 
 // Shift variables by n (possibly negative)
void _shift_variables(bdd& b, int n) {
    // std::cout << "Shifting " << b << " by " << n << std::endl;
    // If n == 0, do nothing
    if(n == 0) return;
    
    // get the list of variables that occurs in b
    int* oldvars;
    int varnum;
    ARIADNE_ASSERT_MSG(bdd_scanset(bdd_support(b), oldvars, varnum) == 0,
        "Error in scanning the variable set.");
    // If varnum is zero, the variable set of b is empty: no shift needed
    if(varnum == 0) return;
    // Check consistency of the shift index
    // std::cout << "Variable support: [" << oldvars[0] << " .. " << oldvars[varnum-1] << "]" << std::endl;
    ARIADNE_ASSERT_MSG(oldvars[0] + n > 0, "Wrong shift index: negative variable number.")
    // Add new vars if necessary
    if(oldvars[varnum-1] + n >= bdd_varnum()) {
        ARIADNE_ASSERT_MSG(bdd_setvarnum(oldvars[varnum-1] + n + 1) == 0, 
            "Cannot add new BDD variables.");
    }
    // create the new variable list
    int* newvars = new int[varnum];
    for(int i = 0; i != varnum; ++i) {
        newvars[i] = oldvars[i] + n;
    }
    // create the bddPair and rename variables        
    bddPair* pPair = bdd_newpair();
    ARIADNE_ASSERT_MSG(bdd_setpairs(pPair, oldvars, newvars, varnum) == 0,
        "Error in creating the bddPair.");
    b = bdd_replace(b, pPair);
    // deallocate temporary data
    free(oldvars);
    free(newvars);
    bdd_freepair(pPair);
}

// Increase the height of the two sets until they have the same root cell,
// that is, the same root_cell_height and the same root_cell_coordinates
void _equalize_root_cell(BDDTreeSet& set1, BDDTreeSet& set2) {
    // check whether they have the same grid
    ARIADNE_ASSERT_MSG(set1.grid() == set2.grid(),
        "Cannot equalize BDDTreeSets based on different grids.");

    // Make the two sets of the same height
    set1.increase_height(set2.root_cell_height());
    set2.increase_height(set1.root_cell_height());

    // Compute the minimum cell height to equalize the root cells
    uint mch = set1.root_cell_height();
    array<int> rcc1 = set1.root_cell_coordinates();
    array<int> rcc2 = set2.root_cell_coordinates();
    uint dim = rcc1.size();
    // increase mch until the two root cell coordinates are the same
    while(rcc1 != rcc2) {
        // determine which dimension to merge 
        uint i = (dim - 1) - (mch % dim);
        // half the i-th coordinates
        rcc1[i] = rcc1[i] / 2;
        rcc2[i] = rcc2[i] / 2;
        // Increase the mch by 1
        mch++;        
    }
    // Equalize the two sets
    set1.increase_height(mch);
    set2.increase_height(mch);    
}

// Test whether the a BDDTreeSet based on root_cell and enabled_cells, where the root variable is root_var
// is a subset of Box. 
tribool _subset(const Box& root_cell, const bdd& enabled_cells, int root_var, uint splitting_coordinate, const Box& box) {
    // std::cout << "_subset(" << root_cell << ", " << enabled_cells << ", " 
    //           << root_var << ", " << splitting_coordinate << ", " << box << ")" << std::endl;
    // if the set is empty, return true
    if(enabled_cells == bddfalse) return true;
    // if the root cell is a subset of box, return true
    if(root_cell.subset(box)) return true;
    // if the bdd is the constant true, return false
    if(enabled_cells == bddtrue) return false;
    
    // Split the root cell and repeat recursively on the two subcells
    std::pair<Box,Box> subcells = root_cell.split(splitting_coordinate);  
    bdd right_branch, left_branch;
    // if the variable labelling the root of the bdd is not root_var, do not split the bdd
    if(bdd_var(enabled_cells) != root_var) {
        left_branch = enabled_cells;
        right_branch = enabled_cells;
    } else {
        left_branch = bdd_low(enabled_cells);
        right_branch = bdd_high(enabled_cells);
    }
    
    int dim = root_cell.dimension();
    if(!_subset(subcells.first, left_branch, root_var+1, (splitting_coordinate+1) % dim, box))
        return false;
    return _subset(subcells.second, right_branch, root_var+1, (splitting_coordinate+1) % dim, box);
}

/************************************* BDDTreeSet **************************************/

BDDTreeSet::BDDTreeSet( )
    : _grid(0)
    , _root_cell_height(0)
    , _root_cell_coordinates()
{
    _initialize_bddlib();
    this->_bdd = bdd_false();
}


BDDTreeSet::BDDTreeSet( const BDDTreeSet & set )
    : _grid(set.grid())
    , _root_cell_height(set.root_cell_height())
    , _root_cell_coordinates(set.root_cell_coordinates())
    , _bdd(set.enabled_cells())
{
}

BDDTreeSet::BDDTreeSet( const Grid& grid, const uint root_cell_height, 
                        const array<int>& root_cell_coordinates, 
                        const bdd& enabled_cells)
    : _grid(grid)
    , _root_cell_height(root_cell_height)
    , _root_cell_coordinates(root_cell_coordinates)
    , _bdd(enabled_cells)
{
    ARIADNE_ASSERT_MSG(root_cell_coordinates.size() == grid.dimension(),
        "root_cell_coordinates and grid must be of the same dimension.");
}


BDDTreeSet::BDDTreeSet( const uint dimension, const bool enable )
    : _grid(dimension)
    , _root_cell_height(0)
    , _root_cell_coordinates(dimension, 0)
{
    _initialize_bddlib();

    if(enable) {
        this->_bdd = bdd_true();
    } else {
        this->_bdd = bdd_false();
    }
}

BDDTreeSet::BDDTreeSet( const Grid& grid, const bool enable )
    : _grid(grid)
    , _root_cell_height(0)
    , _root_cell_coordinates(grid.dimension(), 0)
{
    _initialize_bddlib();

    if(enable) {
        this->_bdd = bdd_true();
    } else {
        this->_bdd = bdd_false();
    }
}

BDDTreeSet* BDDTreeSet::clone() const {
    return new BDDTreeSet( *this );
}

bool BDDTreeSet::empty() const {
    // A zero dimension set is always empty
    if(this->dimension() == 0) return true;
    // otherwise, the set is empty iff the BDD is the constant false
    return (this->_bdd == bddfalse);
}

size_t BDDTreeSet::size() const {
    return _enabled_cells_number(this->enabled_cells(), 0);
}

uint BDDTreeSet::dimension() const {
    return this->_grid.dimension();
}

const Grid& BDDTreeSet::grid() const {
    return this->_grid;
}

uint BDDTreeSet::root_cell_height() const {
    return this->_root_cell_height;
}

array<int> BDDTreeSet::root_cell_coordinates() const {
    return this->_root_cell_coordinates;
}

const bdd& BDDTreeSet::enabled_cells() const {
    return this->_bdd;
}

double BDDTreeSet::measure() const {
    ARIADNE_NOT_IMPLEMENTED;    
}

Box BDDTreeSet::root_cell() const {
    return _compute_root_cell(this->grid(), this->root_cell_height(), this->root_cell_coordinates());
}


Box BDDTreeSet::bounding_box() {
    // minimize the cell height and then return the root_cell
    this->minimize_height();
    return this->root_cell();
}


bool BDDTreeSet::operator==(const BDDTreeSet& set) const {
    if(this->grid() != set.grid())
        return false;
    if(this->root_cell_height() != set.root_cell_height())
        return false;
    if(this->root_cell_coordinates() != set.root_cell_coordinates())
        return false;
        
    return (this->enabled_cells() == set.enabled_cells());                
}

bool BDDTreeSet::operator!=(const BDDTreeSet& set) const {
    return !(this->operator==(set));                
}

bool subset( const BDDTreeSet& set1, const BDDTreeSet& set2 ) {
    // raise an error if the current set is zero dimensional
    ARIADNE_ASSERT_MSG(set1.dimension() != 0, "Cannot test for subset of zero-dimensional sets.");
    
    // set2 must be based on the same grid of the current set
    ARIADNE_ASSERT_MSG(set1.grid() == set2.grid(),
        "Subset inclusion can be tested only for BDDTreeSets based on the same grid.");

    // If set1 is empty the test is true for every set2
    if(set1.empty()) return true;
    // Make copies of set1 and set2 that can be modified.
    BDDTreeSet set3 = set1;
    BDDTreeSet set4 = set2;
    // Equalize the height of the two sets.
    set3.increase_height(set4.root_cell_height());
    set4.increase_height(set3.root_cell_height());
    // since set3 and set4 have the same root_cell_height, we have that set3 is a subset of set4 
    // iff (1) they have the same root_cell_coordinates, AND
    //     (2) set3.enabled_cells() implies set4.enabled_cells() as boolean functions.
    if(set3.root_cell_coordinates() != set4.root_cell_coordinates())
        return false;
    return (bdd_imp(set3.enabled_cells(), set4.enabled_cells()) == bddtrue);
}

bool superset( const BDDTreeSet& set1, const BDDTreeSet& set2 ) {
    return subset(set2, set1);
}

bool disjoint( const BDDTreeSet& set1, const BDDTreeSet& set2 ) {
    // two sets are disjoint iff they do not overlap
    return !overlap(set1, set2);
}

bool overlap( const BDDTreeSet& set1, const BDDTreeSet& set2 ) {
    // raise an error if the current set is zero dimensional
    ARIADNE_ASSERT_MSG(set1.dimension() != 0, "Cannot test for overlapping of zero-dimensional sets.");
    
    // set2 must be based on the same grid of the current set
    ARIADNE_ASSERT_MSG(set1.grid() == set2.grid(),
        "Disjointness can be tested only for BDDTreeSets based on the same grid.");

    // If set1 is empty the test is true iff set2 is empty as well
    if(set1.empty()) return set2.empty();
    // Make copies of set1 and set2 that can be modified.
    BDDTreeSet set3 = set1;
    BDDTreeSet set4 = set2;
    // Equalize the height of the two sets.
    set3.increase_height(set4.root_cell_height());
    set4.increase_height(set3.root_cell_height());
    // since set3 and set4 have the same root_cell_height, we have that set3 overlaps set4 
    // iff (1) they have the same root_cell_coordinates, AND
    //     (2) (set3.enabled_cells() AND set4.enabled_cells()) is not the constant false.
    if(set3.root_cell_coordinates() != set4.root_cell_coordinates())
        return false;
    return (bdd_and(set3.enabled_cells(), set4.enabled_cells()) != bddfalse);    
}


tribool BDDTreeSet::subset( const Box& box ) const {
    // raise an error if the current set is zero dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot test for subset of a zero-dimensional set.");
    
    // the box and the set must have the same dimension
    ARIADNE_ASSERT_MSG(this->dimension() == box.dimension(), "Box and set must have the same dimension.");

    // compute the first splitting coordinate
    uint dim = this->dimension();
    uint height = this->root_cell_height();
    uint splitting_coordinate;
    if(height == 0) {
        splitting_coordinate = 0;
    } else {
        splitting_coordinate = (dim - 1) - ((height-1) % dim);
    }
    return _subset(this->root_cell(), this->enabled_cells(), 0, splitting_coordinate, box);        
}

/*

tribool superset( const Box& box ) const {
    // raise an error if the current set is zero dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot test for superset of a zero-dimensional set.");
    
    // the box and the set must have the same dimension
    ARIADNE_ASSERT_MSG(this->dimension() == box.dimension(), "Box and set must have the same dimension.");
    
    // compute the first splitting coordinate
    uint dim = this->dimension();
    uint height = this->root_cell_height();
    uint splitting_coordinate = (dim - 1) - (height % dim);
    return _superset(this->root_cell(), this->enabled_cells(), 0, spliting_coordinate, box);        
}

tribool disjoint( const Box& box  ) const {
    // raise an error if the current set is zero dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot test for disjoint of a zero-dimensional set.");
    
    // the box and the set must have the same dimension
    ARIADNE_ASSERT_MSG(this->dimension() == box.dimension(), "Box and set must have the same dimension.");
    
    // compute the first splitting coordinate
    uint dim = this->dimension();
    uint height = this->root_cell_height();
    uint splitting_coordinate = (dim - 1) - (height % dim);
    return _disjoint(this->root_cell(), this->enabled_cells(), 0, spliting_coordinate, box);        
}

tribool overlaps( const Box& box ) const {
    // raise an error if the current set is zero dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot test for overlaps of a zero-dimensional set.");
    
    // the box and the set must have the same dimension
    ARIADNE_ASSERT_MSG(this->dimension() == box.dimension(), "Box and set must have the same dimension.");
    
    // compute the first splitting coordinate
    uint dim = this->dimension();
    uint height = this->root_cell_height();
    uint splitting_coordinate = (dim - 1) - (height % dim);
    return _overlaps(this->root_cell(), this->enabled_cells(), 0, spliting_coordinate, box);        
}


    void clear( );
*/

int BDDTreeSet::minimize_height() {
    // Raise an error if the set is zero-dimensional
    ARIADNE_ASSERT_MSG(this->grid().dimension() != 0, "Cannot minimize height of a zero-dimensional set.");
    
    // minimize the height of the bdd without shifting the variables
    // n is the shift index necessary to make the variable numbering consistent again
    int n = _minimize_height(this->_bdd, this->_root_cell_height, this->_root_cell_coordinates, 0);
    // shift variables
    _shift_variables(this->_bdd, n);
    
    return this->root_cell_height();
}

int BDDTreeSet::increase_height(uint new_height) {
    // Raise an error if the set is zero-dimensional
    ARIADNE_ASSERT_MSG(this->grid().dimension() != 0, "Cannot increase height of a zero-dimensional set.");

    // If the new height is smaller or equal to the current one, do nothing
    if(new_height <= this->root_cell_height()) 
        return this->root_cell_height();
        
    // shift variables
    int diff = new_height - this->root_cell_height();
    _shift_variables(this->_bdd, diff);
    
    // increase the height of the bdd
    _increase_height(this->_bdd, this->_root_cell_height, this->_root_cell_coordinates, new_height);

    // return the new root cell height    
    return this->root_cell_height();
}


/*
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
*/

void BDDTreeSet::draw(CanvasInterface& canvas) const {
    ARIADNE_NOT_IMPLEMENTED;
}

/*
    std::ostream& write(std::ostream& os) const;

    void import_from_file(const char*& filename);

    void export_to_file(const char*& filename);

*/

/********************************** Stream operators **********************************/

std::ostream& operator<<(std::ostream& os, const BDDTreeSet& set) {
    return os << "BDDTreeSet( Grid: " << set.grid() << 
                 ", root cell height: " << set.root_cell_height() << 
                 ", root cell coordinates: " << set.root_cell_coordinates() <<
                 ", bdd: " << set.enabled_cells() <<" )";
}

} // namespace Ariadne
