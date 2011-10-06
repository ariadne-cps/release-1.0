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
#include "set_checker.h"

namespace Ariadne {

// definition of the approximation semantics
enum ApprSemantics {AS_OUTER, AS_INNER};


/********************************** Utility functions **********************************/
// The silent Garbage Collector Handler
void _silent_gbchandler(int pre, bddGbcStat *s) {
    // do nothing
}

// An Error Handler for the BDD library that throws an ARIADNE Exception
void _bdd_error_handler(int errcode) {
    if(bdd_errstring(errcode) == NULL) {
        ARIADNE_FAIL_MSG("BDD unknown error " << errcode);
    } else {
        ARIADNE_FAIL_MSG("in BDD library: " << bdd_errstring(errcode));    
    }
}


// Test if BuDDy has been initialized. If not, initialize it.
void _initialize_bddlib() {
    // TO DO: make the different parameters definable by the user.
    if(!bdd_isrunning()) {
        // Set the initial nodenum and cachesize of a medium sized load.
        ARIADNE_ASSERT_MSG(bdd_init(100000, 10000) == 0, "Error in the initialization of the bdd library.");
        // Set the initial number of variables.
        ARIADNE_ASSERT_MSG(bdd_setvarnum(32) == 0, "Error in the initialization of the bdd library.");
        // Set the Garbage Collector Handler to the silent one (removes logging)
        bdd_gbc_hook(_silent_gbchandler);
        // Set the error code handler
        bdd_error_hook(_bdd_error_handler);
    }
}

inline void _ensure_bdd_variables(int numvar) {
    if(bdd_varnum() < numvar) {
        ARIADNE_ASSERT_MSG(bdd_setvarnum(numvar) == 0, "Cannot add new BDD variables.");
    }
}

// compute the coordinate to merge when increasing height by 1
inline uint _merging_coordinate(uint height, uint dim) {
    return (dim - 1) - (height % dim);
}

// compute the coordinate to split when decreasing height by 1
inline uint _splitting_coordinate(uint height, uint dim) {
    if(height == 0) return 0;
    return (dim - 1) - ((height-1) % dim);
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
        uint i = _merging_coordinate(h, dim);
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
size_t _enabled_cells_number(const bdd& b, int level, int mince_depth) {
    // std::cout << "_enabled_cells_number(" << b.id() << "," << level << ", " << mince_depth << ")" << std::endl;
    // Base cases: the size of the constant true/false is 1/0
    if(b == bddfalse) return 0;
    if(b == bddtrue) {
        // the number of cells depends on the mince_depth
        if(mince_depth <= 0) return 1;
        return (1 << mince_depth);
    } 
    // Inductive cases
    int rootlevel = bdd_var2level(bdd_var(b));      // get the level of the root
    // std::cout << "rootlevel = " << rootlevel << std::endl;
    if(rootlevel == level) {                        // no jumps of variables
        bdd left = bdd_low(b);                      // get the left branch
        bdd right = bdd_high(b);                    // get the right branch
        return _enabled_cells_number(left, level+1, mince_depth - 1) +
               _enabled_cells_number(right, level+1, mince_depth - 1);
    }
    if(rootlevel > level) {     // the root variable is greater than the current level
        // hence, the root variable has been deleted by the bdd construction
        // because it pointed to two identical childs.
        return 2*_enabled_cells_number(b, level+1, mince_depth - 1);
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
    // std::cout << "_minimize_height(" <<  b << ", " << root_cell_height << ", " 
    //           << root_cell_coordinates << ", " << root_var << ")" << std::endl;
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
    // determine which dimension to split
    uint dim = root_cell_coordinates.size();
    uint i = _splitting_coordinate(root_cell_height, dim);
    root_cell_height--;
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
    //            << ", " << new_height << ")" << std::endl;
    // if the new height is equal or smaller than the current one do nothing
    if(root_cell_height >= new_height) return root_cell_height;

    // Increase the root cell height by one and recursively call _increase_height
    // compute the new root_cell_coordinate
    // determine which dimension to merge 
    uint dim = root_cell_coordinates.size();
    uint i = _merging_coordinate(root_cell_height, dim);
    // If the current height is an "odd" one, the origin of the grid is shifted
    // shift the i-th coordinate accordingly
    if(((root_cell_height / dim) % 2) == 1) {
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
    // half the i-th coordinate: the correct expression would be
    //      floor(root_cell_coordinate[i] / 2)
    // NOTE THAT floor(-1 / 2) = floor (-0.5) = -1, while -1 / 2 = 0 (integer division)
    // so, if the coordinate is negative, decrease it before halving
    if(root_cell_coordinates[i] < 0) root_cell_coordinates[i]--;
    root_cell_coordinates[i] = root_cell_coordinates[i] / 2;
    
    // recursive call
    return _increase_height(b, root_cell_height, root_cell_coordinates, new_height);
}
 
 // Shift variables by n (possibly negative)
void _shift_variables(bdd& b, int n) {
    // std::cout << "Shifting " << b.id() << " by " << n << std::endl;
    // If n == 0, do nothing
    if(n == 0) return;
    
    // get the list of variables that occurs in b
    int* oldvars;
    int varnum;
    ARIADNE_ASSERT_MSG(bdd_scanset(bdd_support(b), oldvars, varnum) == 0,
        "Error in scanning the variable set.");
    // If varnum is zero, the variable set of b is empty: no shift needed,
    // only check if new vars are necessary
    if(varnum == 0) {
        _ensure_bdd_variables(n + 1);
        return;
    }
    // Check consistency of the shift index
    // std::cout << "Variable support: [" << oldvars[0] << " .. " << oldvars[varnum-1] << "]" << std::endl;
    ARIADNE_ASSERT_MSG(oldvars[0] + n >= 0, "Wrong shift index: negative variable number.")
    // Add new vars if necessary
    _ensure_bdd_variables(oldvars[varnum-1] + n + 1);
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
    // std::cout << "_equalize_root_cell( ... )" << std::endl;
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
        uint i = _merging_coordinate(mch, dim);
        // compute the parity of the merge
        uint p = (mch / dim) % 2;
        // half the i-th cooridnate
        rcc1[i] = (rcc1[i] + p) / 2;
        rcc2[i] = (rcc2[i] + p) / 2;
        // Increase the mch by 1
        mch++;        
    }
    // Equalize the two sets
    set1.increase_height(mch);
    set2.increase_height(mch);    
}

// Test whether the a BDDTreeSet based on root_cell and enabled_cells, where the root variable is root_var
// is a subset of set. 
tribool _subset(const Box& root_cell, const bdd& enabled_cells, int root_var, uint splitting_coordinate, const RegularSetInterface& set) {
    // std::cout << "_subset(" << root_cell << ", " << enabled_cells << ", " 
    //           << root_var << ", " << splitting_coordinate << ", " << set << ")" << std::endl;
    // if the set is empty, return true
    if(enabled_cells == bddfalse) return true;
    tribool test = set.covers(root_cell);
    // if the root cell is definitely a subset of box, return true
    if(definitely(test)) {
        // std::cout << "set definitely covers the root cell, return TRUE." << std::endl;
        return true;
    }
    // if the bdd is the constant true, return test
    if(enabled_cells == bddtrue) {
        // std::cout << "set possibly covers the root cell, return " << test << std::endl;
        return test;
    }
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
    test = _subset(subcells.first, left_branch, root_var+1, (splitting_coordinate+1) % dim, set);
    if(!possibly(test)) {
        // std::cout << "set does not cover the left branch, return FALSE." << std::endl;
        return false;
    }
    test = test && _subset(subcells.second, right_branch, root_var+1, (splitting_coordinate+1) % dim, set);
    // std::cout << "returning " << test << std::endl;
    return test;
}

// Test whether the a BDDTreeSet based on root_cell and enabled_cells, where the root variable is root_var
// is a superset of Box. 
// WARNING: This procedure assumes that box is not empty and a subset of root_cell.
bool _superset(const Box& root_cell, const bdd& enabled_cells, int root_var, uint splitting_coordinate, const Box& box) {
    // std::cout << "_superset(" << root_cell << ", " << enabled_cells << ", " 
    //            << root_var << ", " << splitting_coordinate << ", " << box << ")" << std::endl;
    // if the set is empty, return false
    if(enabled_cells == bddfalse) return false;
    // if the bdd is the constant true, return true
    if(enabled_cells == bddtrue) return true;
    
    // Split the root cell and the box, and repeat recursively on the two subcells
    std::pair<Box,Box> subcells = root_cell.split(splitting_coordinate);  
    Box left_box = box;
    left_box[splitting_coordinate] = intersection(box[splitting_coordinate], subcells.first[splitting_coordinate]); 
    // std::cout << "left_box = " << left_box << std::endl;
    Box right_box = box;
    right_box[splitting_coordinate] = intersection(box[splitting_coordinate], subcells.second[splitting_coordinate]); 
    // std::cout << "right_box = " << right_box << std::endl;
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
    // If the left box is not empty, repeat the test on the left branch
    if(!left_box[splitting_coordinate].empty()) {
        if(!_superset(subcells.first, left_branch, root_var+1, (splitting_coordinate+1) % dim, left_box)) {
            return false;
        }
    }
     // If the right box is not empty, repeat the test on the right branch
    if(!right_box[splitting_coordinate].empty()) {
        return (_superset(subcells.second, right_branch, root_var+1, (splitting_coordinate+1) % dim, right_box));
    }

    return true;
}

// Test whether the a BDDTreeSet based on root_cell and enabled_cells, where the root variable is root_var
// is a disjoint from set. 
tribool _disjoint(const Box& root_cell, const bdd& enabled_cells, int root_var, uint splitting_coordinate, const RegularSetInterface& set) {
    // std::cout << "_disjoint(" << root_cell << ", " << enabled_cells << ", " 
    //           << root_var << ", " << splitting_coordinate << ", " << set << ")" << std::endl;
    // if the set is empty, return true
    if(enabled_cells == bddfalse) return true;
    tribool test = set.disjoint(root_cell);
    // if the root cell is definitely disjoint from set, return true
    if(definitely(test)) return true;
    // the set is possibly disjoint from the root_cell
    // if the bdd is the constant true, return the result of the test
    if(enabled_cells == bddtrue) return test;
    
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
    // test the left branch
    test = _disjoint(subcells.first, left_branch, root_var+1, (splitting_coordinate+1) % dim, set);
    if(!possibly(test)) return false;
    return test && _disjoint(subcells.second, right_branch, root_var+1, (splitting_coordinate+1) % dim, set);
}

// recursive function that adjoin a box to a bdd with a maximum depth
bdd _adjoin_over_approximation(const Box& box, const bdd& enabled_cells, const Box& root_cell,
                               uint depth, uint splitting_coordinate, int root_var)
{
    // std::cout << "_adjoin_over_approximation(" << box << ", " << enabled_cells << ", "
    //           << root_cell << ", " << depth << ", " << splitting_coordinate << ", " << root_var 
    //           << ")" << std::endl;
    // if the bdd is the constant true, do nothing
    if(enabled_cells == bddtrue) return bddtrue;
    // if the root cell is a subset of box, return the constant true bdd
    if(root_cell.subset(box)) return bddtrue;
    // if the root cell is disjoint from the box, do nothing
    if(root_cell.disjoint(box)) return enabled_cells;
    // now we can assume that box overlaps the root cell
    // if the depth is zero, return the constant true bdd
    if(depth == 0) return bddtrue;
    
    // if depth > 0, split the root cell and the bdd, and continue recursively
    std::pair <Box, Box> split_cells = root_cell.split(splitting_coordinate);
    bdd right_branch, left_branch;
    // if the variable labelling the root of the bdd is not root_var, do not split the bdd
    if(enabled_cells == bddfalse || bdd_var(enabled_cells) != root_var) {
        // std::cout << "bdd_var is different from root_var, do no split the bdd." << std::endl;
        left_branch = enabled_cells;
        right_branch = enabled_cells;
    } else {
        // std::cout << "bdd_var is equal to root_var, split the bdd." << std::endl;
        left_branch = bdd_low(enabled_cells);
        right_branch = bdd_high(enabled_cells);
    }
    left_branch = _adjoin_over_approximation(box, left_branch, split_cells.first, depth - 1,
                                    (splitting_coordinate + 1) % box.dimension(), root_var + 1);
    right_branch = _adjoin_over_approximation(box, right_branch, split_cells.second, depth - 1,
                                    (splitting_coordinate + 1) % box.dimension(), root_var + 1);
    return (bdd_nithvar(root_var) & left_branch) | (bdd_ithvar(root_var) & right_branch);
}

// recursive function that adjoin an outer approximation of set to a bdd with a maximum depth
bdd _adjoin_outer_approximation(const CompactSetInterface& set, const bdd& enabled_cells, const Box& root_cell,
                               uint depth, uint splitting_coordinate, int root_var)
{
    // std::cout << "_adjoin_outer_approximation(" << set << ", " << enabled_cells << ", "
    //           << root_cell << ", " << depth << ", " << splitting_coordinate << ", " << root_var 
    //           << ")" << std::endl;
    // if the bdd is the constant true, do nothing
    if(enabled_cells == bddtrue) return bddtrue;
    // if the root cell is disjoint from the set, do nothing
    if(definitely(set.disjoint(root_cell))) {
        // std::cout << "set is disjoint from the current cell, skipping." << std::endl;
        return enabled_cells;
    }
    // std::cout << "set overlaps the current cell, go on." << std::endl;
    // now we can assume that set overlaps the root cell
    // if the depth is zero, return the constant true bdd
    if(depth == 0) {
        // std::cout << "depth is zero, mark the current cell." << std::endl;
        return bddtrue;
    }
    // if depth > 0, split the root cell and the bdd, and continue recursively
    std::pair <Box, Box> split_cells = root_cell.split(splitting_coordinate);
    bdd right_branch, left_branch;
    // if the variable labelling the root of the bdd is not root_var, do not split the bdd
    if(enabled_cells == bddfalse || bdd_var(enabled_cells) != root_var) {
        // std::cout << "bdd_var is different from root_var, do no split the bdd." << std::endl;
        left_branch = enabled_cells;
        right_branch = enabled_cells;
    } else {
        // std::cout << "bdd_var is equal to root_var, split the bdd." << std::endl;
        left_branch = bdd_low(enabled_cells);
        right_branch = bdd_high(enabled_cells);
    }
    left_branch = _adjoin_outer_approximation(set, left_branch, split_cells.first, depth - 1,
                                    (splitting_coordinate + 1) % set.dimension(), root_var + 1);
    right_branch = _adjoin_outer_approximation(set, right_branch, split_cells.second, depth - 1,
                                    (splitting_coordinate + 1) % set.dimension(), root_var + 1);
    return (bdd_nithvar(root_var) & left_branch) | (bdd_ithvar(root_var) & right_branch);
}

// recursive function that adjoin a lower approximation of a set to a bdd with a maximum depth
bdd _adjoin_lower_approximation(const OvertSetInterface& set, const bdd& enabled_cells, const Box& root_cell,
                               uint depth, uint splitting_coordinate, int root_var)
{
    // std::cout << "_adjoin_lower_approximation(" << set << ", " << enabled_cells << ", "
    //           << root_cell << ", " << depth << ", " << splitting_coordinate << ", " << root_var 
    //           << ")" << std::endl;
    // if the bdd is the constant true, do nothing
    if(enabled_cells == bddtrue) return bddtrue;
    tribool test = set.overlaps(root_cell);
    // if the root cell is definitely disjoint from the set, do nothing
    if(!possibly(test)) {
        // std::cout << "set is definitley disjoint from the current cell, skipping." << std::endl;
        return enabled_cells;
    }
    // std::cout << "set possibly overlaps the current cell, go on." << std::endl;
    // if the depth is zero, return the true bdd if the set definitely overlaps the root cell
    // otherwise, do not mark any new cell
    if(depth == 0) {
        if(definitely(test)) {
            // std::cout << "depth is zero and the set definitely overlaps the current cell, mark it." << std::endl;
            return bddtrue;
        } else {
            // std::cout << "depth is zero but the set possibly do not overlap the current cell, skip it." << std::endl;        
            return enabled_cells;
        }
    }
    // if depth > 0, split the root cell and the bdd, and continue recursively
    std::pair <Box, Box> split_cells = root_cell.split(splitting_coordinate);
    bdd right_branch, left_branch;
    // if the variable labelling the root of the bdd is not root_var, do not split the bdd
    if(enabled_cells == bddfalse || bdd_var(enabled_cells) != root_var) {
        // std::cout << "bdd_var is different from root_var, do no split the bdd." << std::endl;
        left_branch = enabled_cells;
        right_branch = enabled_cells;
    } else {
        // std::cout << "bdd_var is equal to root_var, split the bdd." << std::endl;
        left_branch = bdd_low(enabled_cells);
        right_branch = bdd_high(enabled_cells);
    }
    left_branch = _adjoin_lower_approximation(set, left_branch, split_cells.first, depth - 1,
                                    (splitting_coordinate + 1) % set.dimension(), root_var + 1);
    right_branch = _adjoin_lower_approximation(set, right_branch, split_cells.second, depth - 1,
                                    (splitting_coordinate + 1) % set.dimension(), root_var + 1);
    return (bdd_nithvar(root_var) & left_branch) | (bdd_ithvar(root_var) & right_branch);
}

// recursive function that adjoin an inner approximation of a set to a bdd with a maximum depth
bdd _adjoin_inner_approximation(const OpenSetInterface& set, const bdd& enabled_cells, const Box& root_cell,
                               uint depth, uint splitting_coordinate, int root_var)
{
    // std::cout << "_adjoin_outer_approximation(" << set << ", " << enabled_cells << ", "
    //            << root_cell << ", " << depth << ", " << splitting_coordinate << ", " << root_var 
    //            << ")" << std::endl;
    // if the bdd is the constant true, do nothing
    if(enabled_cells == bddtrue) return bddtrue;
    // if the set does not overlap the root cell, do nothing
    if(!possibly(set.overlaps(root_cell))) {
        // std::cout << "set definitely not overlap the current cell, skipping." << std::endl;
        return enabled_cells;
    }
    // if the set definitely covers the current cell, mark it
    if(definitely(set.covers(root_cell))) {
        // std::cout << "set definitely covers the current cell, mark it." << std::endl;
        return bddtrue;
    }    
    // std::cout << "set possibly covers the current cell, go on." << std::endl;
    // if the depth is zero, do not mark any new cell
    if(depth == 0) {
        // std::cout << "depth is zero, skip." << std::endl;
        return enabled_cells;
    }
    // if depth > 0, split the root cell and the bdd, and continue recursively
    std::pair <Box, Box> split_cells = root_cell.split(splitting_coordinate);
    bdd right_branch, left_branch;
    // if the variable labelling the root of the bdd is not root_var, do not split the bdd
    if(enabled_cells == bddfalse || bdd_var(enabled_cells) != root_var) {
        // std::cout << "bdd_var is different from root_var, do no split the bdd." << std::endl;
        left_branch = enabled_cells;
        right_branch = enabled_cells;
    } else {
        // std::cout << "bdd_var is equal to root_var, split the bdd." << std::endl;
        left_branch = bdd_low(enabled_cells);
        right_branch = bdd_high(enabled_cells);
    }
    left_branch = _adjoin_inner_approximation(set, left_branch, split_cells.first, depth - 1,
                                    (splitting_coordinate + 1) % set.dimension(), root_var + 1);
    right_branch = _adjoin_inner_approximation(set, right_branch, split_cells.second, depth - 1,
                                    (splitting_coordinate + 1) % set.dimension(), root_var + 1);
    return (bdd_nithvar(root_var) & left_branch) | (bdd_ithvar(root_var) & right_branch);
}

// recursive function that restrict the bdd to the cells that possibly overlaps with a set, for AS_OUTER semantics,
// or to the cells that are definitely inside the set, for AS_INNER semantics
bdd _approximate_restrict(const OpenSetInterface& set, const bdd& enabled_cells, const Box& root_cell,
                       const uint depth, const uint splitting_coordinate, const int root_var, const ApprSemantics semantics)
{
    // std::cout << "_approximate_restrict(" << set << ", " << enabled_cells << ", "
    //           << root_cell << ", " << depth << ", " << splitting_coordinate << ", " << root_var << ", "
    //           << (semantics == AS_INNER ? "AS_INNER" : "AS_OUTER") << ")" << std::endl;
    // if the bdd is the constant false, do nothing
    if(enabled_cells == bddfalse) return enabled_cells;
    // if the set definitely not overlap the root cell, return the empty bdd
    if(!possibly(set.overlaps(root_cell))) {
        // std::cout << "set definitely not overlap the current cell, return false." << std::endl;
        return bddfalse;
    }
    // if set covers the current cell, return the current bdd
    if(definitely(set.covers(root_cell))) {
        // std::cout << "the set covers the current cell, return the current bdd." << std::endl;
        return enabled_cells;
    }    
    // if the depth is zero, return the empty bdd for inner semantics, and the current bdd for outer semantics
    if(depth == 0) {
        if(semantics == AS_INNER) {
            // std::cout << "depth is 0 and the set possibly not cover the current cell, skip it." << std::endl;
            return bddfalse;
        } else {
            // std::cout << "depth is 0 and the set possibly overlaps the current cell, return the current bdd." << std::endl;
            return enabled_cells;        
        }
    }
    // std::cout << "set possibly overlaps the current cell, but possibly not cover it, go on." << std::endl;
    // split the root cell and the bdd, and continue recursively
    std::pair <Box, Box> split_cells = root_cell.split(splitting_coordinate);
    bdd right_branch, left_branch;
    // if the variable labelling the root of the bdd is not root_var, do not split the bdd
    if(enabled_cells == bddtrue || bdd_var(enabled_cells) != root_var) {
        // std::cout << "bdd_var is different from root_var, do no split the bdd." << std::endl;
        left_branch = enabled_cells;
        right_branch = enabled_cells;
    } else {
        // std::cout << "bdd_var is equal to root_var, split the bdd." << std::endl;
        left_branch = bdd_low(enabled_cells);
        right_branch = bdd_high(enabled_cells);
    }
    left_branch = _approximate_restrict(set, left_branch, split_cells.first, depth - 1,
                                    (splitting_coordinate + 1) % set.dimension(), root_var + 1, semantics);
    right_branch = _approximate_restrict(set, right_branch, split_cells.second, depth - 1,
                                    (splitting_coordinate + 1) % set.dimension(), root_var + 1, semantics);
    return (bdd_nithvar(root_var) & left_branch) | (bdd_ithvar(root_var) & right_branch);
}

// recursive function that restrict the bdd to the cells that possibly respect a property, for AS_OUTER semantics,
// or to the cells that definitely respect it, for AS_INNER semantics
bdd _approximate_restrict(const SetCheckerInterface& checker, const bdd& enabled_cells, const Box& root_cell,
                       const uint depth, const uint splitting_coordinate, const int root_var, const ApprSemantics semantics)
{
    // std::cout << "_approximate_restrict(checker, " << enabled_cells << ", "
    //           << root_cell << ", " << depth << ", " << splitting_coordinate << ", " << root_var << ", "
    //           << (semantics == AS_INNER ? "AS_INNER" : "AS_OUTER") << ")" << std::endl;
    // if the bdd is the constant false, do nothing
    if(enabled_cells == bddfalse) return enabled_cells;
    // get the value of the property for the root cell
    tribool prop = checker.check(root_cell);
    // if the current cell definitely not respect the property, return the empty bdd
    if(!possibly(prop)) {
        // std::cout << "the current cell do not respect the property, return false." << std::endl;
        return bddfalse;
    }
    // if the current cell definitely respects the property, return the current bdd
    if(definitely(prop)) {
        // std::cout << "the current cell definitely respects the property, return the current bdd." << std::endl;
        return enabled_cells;
    }    
    // if the depth is zero, return the empty bdd for inner semantics, and the current bdd for outer semantics
    if(depth == 0) {
        if(semantics == AS_INNER) {
            // std::cout << "depth is 0 and the current cell possibly do not respect the property, skip it." << std::endl;
            return bddfalse;
        } else {
            // std::cout << "depth is 0 and the current cell possibly respects the property, return the current bdd." << std::endl;
            return enabled_cells;        
        }
    }
    // std::cout << "the current cell possibly respects the property, go on." << std::endl;
    // split the root cell and the bdd, and continue recursively
    std::pair <Box, Box> split_cells = root_cell.split(splitting_coordinate);
    bdd right_branch, left_branch;
    // if the variable labelling the root of the bdd is not root_var, do not split the bdd
    if(enabled_cells == bddtrue || bdd_var(enabled_cells) != root_var) {
        // std::cout << "bdd_var is different from root_var, do no split the bdd." << std::endl;
        left_branch = enabled_cells;
        right_branch = enabled_cells;
    } else {
        // std::cout << "bdd_var is equal to root_var, split the bdd." << std::endl;
        left_branch = bdd_low(enabled_cells);
        right_branch = bdd_high(enabled_cells);
    }
    uint dim = root_cell.dimension();
    left_branch = _approximate_restrict(checker, left_branch, split_cells.first, depth - 1,
                                    (splitting_coordinate + 1) % dim, root_var + 1, semantics);
    right_branch = _approximate_restrict(checker, right_branch, split_cells.second, depth - 1,
                                    (splitting_coordinate + 1) % dim, root_var + 1, semantics);
    return (bdd_nithvar(root_var) & left_branch) | (bdd_ithvar(root_var) & right_branch);
}


// Function that increments a BDDTreeSet iterator
void _compute_next_cell(std::vector< PathElement >& path, int mince_depth) {
    // std::cout << "_compute_next_cell( path , " << mince_depth << ")" << std::endl;
    // if the path is empty, do nothing and return
    if(path.empty()) {
        return;
    }
    // get the last step of the path
    PathElement tail = path.back();
    // std::cout << "  tail = " << tail << std::endl;
    uint dim = tail.cell.dimension();
    uint i;
    // if the status is RIGHT, both child have been already explored, backtrack
    if(tail.status == PE_RIGHT) {
        // std::cout << "  status is PE_RIGHT, backtracking." << std::endl;        
        path.pop_back();
        _compute_next_cell(path, mince_depth + 1);
        return;
    }
    // if the status is LEFT, explore the right child
    if(tail.status == PE_LEFT) {
        // std::cout << "  status is PE_LEFT, explore right child." << std::endl;        
        // set the status to RIGHT
        tail.status = PE_RIGHT;
        path.pop_back();
        path.push_back(tail);
        if(tail.obdd != bddtrue) {
            // if the var labeling the bdd correspond to che current level, get the right child
            if(bdd_var(tail.obdd) == tail.root_var) {
                tail.obdd = bdd_high(tail.obdd);
            }
        }
        // get splitting coordinate
        i = tail.split_coordinate;
        // compute new cell by splitting coordinate i
        tail.cell[i].set_lower((tail.cell[i].lower() + tail.cell[i].upper())/2.0);
        // set status to new
        tail.status = PE_NEW;
        tail.split_coordinate = (i + 1) % dim;
        tail.root_var = tail.root_var + 1;
        // append the new step of the path and continue recursively
        path.push_back(tail);
        _compute_next_cell(path, mince_depth - 1);
        return;
    }
    // if the status is NEW, we are exploring the node for the first time
    // std::cout << "  status is PE_NEW ";        

    // if the bdd is the constant false, backtrack
    if(tail.obdd == bddfalse) {
        // std::cout << "and the current bdd is FALSE, backtracking" << std::endl;
        path.pop_back();
        _compute_next_cell(path, mince_depth + 1);
        return;
    }
    // if the bdd is the constant true, we are in a cell:
    // if mince_depth is zero or negative, fix the flag to true and return
    // otherwise, go on and split
    if(tail.obdd == bddtrue && mince_depth <= 0) {
        // std::cout << ", the current bdd is TRUE and depth is greater or equal do mince_depth, exiting" << std::endl;
        // set the status to RIGHT
        tail.status = PE_RIGHT;
        path.pop_back();
        path.push_back(tail);
        return;
    }

    // either the bdd is different from true and false, or we have to mince the cell: explore the left child
    // std::cout << "and the current bdd is neither TRUE nor FALSE, or we have to mince further, go on" << std::endl;

    // set the status to left
    tail.status = PE_LEFT;
    path.pop_back();
    path.push_back(tail);
    if(tail.obdd != bddtrue) {
        // if the var labeling the bdd correspond to che current level, get the right child
        if(bdd_var(tail.obdd) == tail.root_var) {
            tail.obdd = bdd_low(tail.obdd);
        }
    }
    // get splitting coordinate
    i = tail.split_coordinate;
    // compute new cell by splitting coordinate i
    tail.cell[i].set_upper((tail.cell[i].lower() + tail.cell[i].upper())/2.0);
    // set status to new
    tail.status = PE_NEW;
    tail.split_coordinate = (i + 1) % dim;
    tail.root_var = tail.root_var + 1;
    // append the new step of the path and continue recursively
    path.push_back(tail);
    _compute_next_cell(path, mince_depth - 1);
}

// Function that computes the bounding_box
Box _bounding_box(bdd enabled_cells, Box const& cell, uint dim, uint split, int root_var) {
    // std::cout << "_bounding_box(" << enabled_cells << ", " << cell << ", " << dim << ", " 
    //           << split << ", " << root_var << ")" << std::endl;
    // if the current bdd is FALSE, return the empty box
    if(enabled_cells == bddfalse) {
        // std::cout << "The current bdd is FALSE, return the empty box." << std::endl;
        return Box::empty_box(dim);
    }
    // if the current bdd is TRUE, return the current cell
    if(enabled_cells == bddtrue) {
        // std::cout << "The current bdd is TRUE, return the current cell." << std::endl;
        return cell;
    }
    // split the current cell
    std::pair<Box, Box> cells = cell.split(split);
    // if the root_var is different from the var labelling the bdd, the bounding_box
    // can be obtained by enlarging the bbox of one of the child on the split coordinate
    if(bdd_var(enabled_cells) != root_var) {
        // get the bounding box of the left child
        Box bbox = _bounding_box(enabled_cells, cells.first, dim, (split + 1) % dim, root_var+1);
        bbox[split] = cell[split];
        // std::cout << "The variable labeling the bdd is different from " << root_var 
        //           << ", return " << bbox << std::endl;
        return bbox;
    }
    // otherwise, get the bounding box of two childs and return the convex hull
    Box lbox = _bounding_box(bdd_low(enabled_cells), cells.first, dim, (split + 1) % dim, root_var+1);
    Box rbox = _bounding_box(bdd_high(enabled_cells), cells.second, dim, (split + 1) % dim, root_var+1);
    
    lbox = hull(lbox, rbox);
    // std::cout << "Return the convex hull: " << lbox << std::endl;
    return lbox;
}

// function that project down a box. 
// WARNING: the procedure does not check if the set of indices is consistent.
Box _unchecked_project_down(const Box& original_box, const Vector<uint>& indices)
{
    Box new_box(indices.size());
	for (uint i=0; i < indices.size(); ++i) {
		new_box[i] = original_box[indices[i]];
	}
	return new_box;
}


// Function that returns the new mince_depth after projecting down to indices
// WARNING: the procedure does not check if the set of indices is consistent.
int _project_mince_depth(const int old_mince_depth, uint old_dim, const Vector<uint>& indices) {
    // std::cout << "_project_mince_depth(" << old_mince_depth << ", " << old_dim << ", " << indices << ")" << std::endl;
    // the original set was not minced
    if(old_mince_depth == -1) return -1;
    // if the original set was minced, compute the new mince_depth
    uint new_dim = indices.size();
    int new_mince_depth = (old_mince_depth / old_dim) * new_dim;
    // if old_mince_depth is not a multiple of old_dim, increase new_mince_depth by new_dim
    if((old_mince_depth % old_dim) != 0) {
        new_mince_depth = new_mince_depth + new_dim;
    }
    return new_mince_depth;
}

inline bool _present(int* vars, int varnum, int var) {
    for(int i = 0; i < varnum; ++i) {
        if(vars[i] == var) {
            return true;
        }
        if(vars[i] > var) {
            return false;
        }
    }
    return false;
}

// Compute the term to duplicate variable firstvar with secondvar in old_bdd.
inline bdd _duplicate_var(const bdd& old_bdd, int firstvar, int secondvar, int lastvar, int mince_depth) {
    bdd res = bdd_biimp(bdd_ithvar(firstvar), bdd_ithvar(secondvar));
    // if the index of the variable is greater than the mince_depth,
    // add the term "OR FORALL (firstvar,...,lastvar).old_bdd
    if(firstvar > mince_depth) {
        bdd vars = bddtrue;
        for(int i = firstvar; i <= lastvar; ++i) {
            vars = vars & bdd_ithvar(i);
        }
        res = res | bdd_forall(old_bdd, vars);
    }
    return res;
}


// remove and duplicate variables in old_bdd following indices. 
// During the procedure indices is modified to avoid duplicated entries.
// WARNING: the procedure does not check if the set of indices is consistent.
bdd _project_down_variables(const bdd& old_bdd, uint old_dim, 
                            const Vector<uint>& indices, int mince_depth) 
{
    // std::cout << "_project_down_variables(" << old_bdd << ", "  
    //           << old_dim << ", " << indices << ")" << std::endl; 
    // the first free var is x_{old_depth}
    int* vars;
    int varnum;
    ARIADNE_ASSERT_MSG(bdd_scanset(bdd_support(old_bdd), vars, varnum) == 0,
        "Error in scanning the variable set.");
    int freevar = 0;
    int lastvar = -1;
    if(varnum > 0) {
        lastvar = vars[varnum-1];
        freevar = lastvar + 1;
    }
    // assure that enough variables are present (in the worst case all vars are duplicated).
    _ensure_bdd_variables(freevar + indices.size() * (freevar / old_dim + 1)); 
    array<bool> keep(old_dim, false);
    bddPair* renPairs = bdd_newpair();
    ARIADNE_ASSERT_MSG(renPairs != NULL, "Cannot allocate a new bddPair.");
    bdd new_bdd = old_bdd;
    bdd remove = bddtrue;
    bdd duplicate = bddtrue;
    uint new_dim = indices.size();
    // scan indices to determine which variables to keep and which one to duplicate    
    for(uint i = 0; i != indices.size(); ++i) {
        // std::cout << "indices[" << i << "] = " << indices[i] << std::endl;
        if(!keep[indices[i]]) {
            // std::cout << "the variable must be kept." << std::endl;
            // the variable indices[i] must be kept
            keep[indices[i]] = true;
            // add the correct pairs for the renaming
            int l = 0;
            for(int k = 0; (int)indices[i] + k <= lastvar ; k += old_dim, l += new_dim) {
                // std::cout << "adding the pair (" << indices[i] + k 
                //           << ", " << i + l << ")" << std::endl;
                ARIADNE_ASSERT_MSG(bdd_setpair(renPairs, indices[i] + k, i + l) == 0,
                         "Error in adding a variable pair for the rename.");                    
            }
        } else {
            // std::cout << "the variable is duplicated." << std::endl;
            // the variable indices[i] is duplicated in indices
            // add the pairs for the renaming and duplicate the var in the bdd
            uint l = 0;
            for(int k = 0; (int)indices[i] + k <= lastvar ; k += old_dim, l += new_dim) {
                // duplicate variable x = indices[i] + k with variable y = freevar
                // this is done by adding the following term to the BDD:
                // (FORALL (x..x_max).old_bdd) OR (x <=> y)
                // std::cout << "duplicating variable " << indices[i] + k << " with " 
                //           << freevar << std::endl;
                duplicate = duplicate & _duplicate_var(old_bdd, indices[i] + k, freevar, lastvar, mince_depth);
                // std::cout << "duplicate = " << duplicate << std::endl;
                // add the pair for the renaming x_depth / x_{i+k*dim}
                // std::cout << "adding the pair (" << freevar  
                //           << ", " << i + l << ") for renaming" << std::endl;
                ARIADNE_ASSERT_MSG(bdd_setpair(renPairs, freevar, i + l) == 0,
                        "Error in adding a variable pair for the rename.");
                // increase freevar to get the next free variable index
                freevar++;
            }
        }   // if(!keep[indices[i]])
    }   // for(uint i = 0; i != indices.size(), ++i)  
    // std::cout << "duplicate = " << duplicate << std::endl;
    // scan keep to obtain which variable must be deleted
    for(int i = 0; i < keep.size(); ++i) {
        // if the dimension i must be removed, add all corresponding variables to remove
        if(!keep[i]) {
            for(int k = 0; (int)i + k <= lastvar ; k += old_dim) {
                remove = remove & bdd_ithvar(i + k);
            }
        }
    }
    // std::cout << "remove = " << remove << std::endl;
    // remove variables by applying the existential quantification 
    new_bdd = bdd_appex(new_bdd, duplicate, bddop_and, remove);
    // std::cout << "new_bdd = " << new_bdd << std::endl;
    // rename the variable to obtain the correct order
    new_bdd = bdd_replace(new_bdd, renPairs);
    // std::cout << "new_bdd = " << new_bdd << std::endl;
    // free renPairs
    bdd_freepair(renPairs);
    
    return new_bdd;
}


/************************************* BDDTreeSet **************************************/

BDDTreeSet::BDDTreeSet( )
    : _grid(0)
    , _root_cell_height(0)
    , _root_cell_coordinates()
    , _mince_depth(-1)
{
    _initialize_bddlib();
    this->_bdd = bdd_false();
}


BDDTreeSet::BDDTreeSet( const BDDTreeSet & set )
    : _grid(set.grid())
    , _root_cell_height(set.root_cell_height())
    , _root_cell_coordinates(set.root_cell_coordinates())
    , _bdd(set.enabled_cells())
    , _mince_depth(set.mince_depth())
{
}

BDDTreeSet::BDDTreeSet( const Grid& grid, const uint root_cell_height, 
                        const array<int>& root_cell_coordinates, 
                        const bdd& enabled_cells, const int mince_depth)
    : _grid(grid)
    , _root_cell_height(root_cell_height)
    , _root_cell_coordinates(root_cell_coordinates)
    , _bdd(enabled_cells)
    , _mince_depth(mince_depth)
{
    ARIADNE_ASSERT_MSG(root_cell_coordinates.size() == grid.dimension(),
        "root_cell_coordinates and grid must be of the same dimension.");
}


BDDTreeSet::BDDTreeSet( const uint dimension, const bool enable )
    : _grid(dimension)
    , _root_cell_height(0)
    , _root_cell_coordinates(dimension, 0)
    , _mince_depth(-1)
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
    , _mince_depth(-1)
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
    // compute mince depth
    int mince_depth = -1;
    if(this->mince_depth() >= 0) {
        mince_depth = this->root_cell_height() + this->mince_depth();
    }
    return _enabled_cells_number(this->enabled_cells(), 0, mince_depth);
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

uint BDDTreeSet::depth() const {
    // get the index of the smallest var in the bdd
    int* vars;
    int varnum;
    ARIADNE_ASSERT_MSG(bdd_scanset(bdd_support(this->_bdd), vars, varnum) == 0,
        "Error in scanning the variable set.");
    if(varnum != 0) {
        // The depth of the tree is the index of the last variable of the support + 1
        // std::cout << "vars = " << vars[0] << " ... " << vars[varnum-1] << std::endl;
        varnum = vars[varnum-1] + 1;
    }
    // Remove the variables that are necessary to reach the zero level
    varnum = varnum - this->root_cell_height();
    // If varnum is negative, the depth is zero
    if(varnum < 0) varnum = 0;
    // free vars
    free(vars);    
    // the depth is the maximum between the mince_depth and the real depth of the tree
    return max(varnum, this->mince_depth());
}

array<int> BDDTreeSet::root_cell_coordinates() const {
    return this->_root_cell_coordinates;
}

int BDDTreeSet::mince_depth() const {
    return this->_mince_depth;
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


Box BDDTreeSet::bounding_box() const {
    if(this->empty()) return Box::empty_box(this->dimension());
    uint dim = this->dimension();
    // get the first splitting coordinate
    uint split = _splitting_coordinate(this->root_cell_height(), dim);
    return _bounding_box(this->enabled_cells(), this->root_cell(), dim, split, 0);
}


bool BDDTreeSet::operator==(const BDDTreeSet& set) const {
    if(this->grid() != set.grid())
        return false;
    if(this->root_cell_height() != set.root_cell_height())
        return false;
    if(this->root_cell_coordinates() != set.root_cell_coordinates())
        return false;
    if(this->mince_depth() != set.mince_depth())    
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

tribool BDDTreeSet::subset( const BDDTreeSet& other ) const {
    return Ariadne::subset(*this, other);
}

tribool BDDTreeSet::superset( const BDDTreeSet& other ) const{
    return Ariadne::superset(*this, other);
}

tribool BDDTreeSet::disjoint( const BDDTreeSet& other  ) const{
    return Ariadne::disjoint(*this, other);
}

tribool BDDTreeSet::overlaps( const BDDTreeSet& other ) const{
    return Ariadne::overlap(*this, other);
}

tribool BDDTreeSet::subset( const Box& box ) const {
    // raise an error if the current set is zero dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot test for subset of a zero-dimensional set.");
    
    // the box and the set must have the same dimension
    ARIADNE_ASSERT_MSG(this->dimension() == box.dimension(), "Box and set must have the same dimension.");

    // compute the first splitting coordinate
    uint dim = this->dimension();
    uint height = this->root_cell_height();
    uint splitting_coordinate = _splitting_coordinate(height, dim);
    return _subset(this->root_cell(), this->enabled_cells(), 0, splitting_coordinate, box);        
}

tribool BDDTreeSet::superset( const Box& box ) const {
    // raise an error if the current set is zero dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot test for superset of a zero-dimensional set.");

    // If box is empty, the test is true
    if(box.empty()) return true;
    
    // the box and the set must have the same dimension
    ARIADNE_ASSERT_MSG(this->dimension() == box.dimension(), "Box and set must have the same dimension.");
    
    // If the box is not a subset of the root cell, the test is false
    if(!box.subset(this->root_cell())) return false;
    
    // compute the first splitting coordinate
    uint dim = this->dimension();
    uint height = this->root_cell_height();
    uint splitting_coordinate = _splitting_coordinate(height, dim);
    return _superset(this->root_cell(), this->enabled_cells(), 0, splitting_coordinate, box);        
}


tribool BDDTreeSet::disjoint( const Box& box  ) const {
    // raise an error if the current set is zero dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot test for disjoint of a zero-dimensional set.");

    // If box is empty, the test is true
    if(box.empty()) return true;
    
    // the box and the set must have the same dimension
    ARIADNE_ASSERT_MSG(this->dimension() == box.dimension(), "Box and set must have the same dimension.");
    
    // compute the first splitting coordinate
    uint dim = this->dimension();
    uint height = this->root_cell_height();
    uint splitting_coordinate = _splitting_coordinate(height, dim);
    return _disjoint(this->root_cell(), this->enabled_cells(), 0, splitting_coordinate, box);        
}

tribool BDDTreeSet::overlaps( const Box& box ) const {
    return !this->disjoint(box);        
}

void BDDTreeSet::mince( const uint subdiv ) {
    this->_mince_depth = subdiv * this->dimension();
};

void BDDTreeSet::recombine() {
    this->_mince_depth = -1;
};


void BDDTreeSet::clear( ) {
    this->_root_cell_height = 0;
    this->_root_cell_coordinates.fill(0);
    this->_bdd = bddfalse;
    this->_mince_depth = -1;
}


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
    // std::cout << "increase_height(" <<  new_height << ")" << std::endl;
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

int BDDTreeSet::increase_height(const Box& box) {
    // std::cout << "increase_height(" <<  box << ")" << std::endl;

    // Raise an error if the set is zero-dimensional
    ARIADNE_ASSERT_MSG(this->grid().dimension() != 0, "Cannot increase height of a zero-dimensional set.");

    // Get the root gird
    uint height = this->root_cell_height();
    Grid root_grid = _compute_root_grid(this->grid(), height);
    // compute the root cell from the coordinates
    array<int> coordinates = this->root_cell_coordinates();
    Box root_cell = root_grid.cell(coordinates);

    uint dim = this->dimension();
    int i = 0;
    while(!definitely(root_cell.covers(box))) {
        // std::cout << "height = " << height << ", root_cell = " << root_cell << std::endl;        
        // determine which dimension to merge 
        i = _merging_coordinate(height, dim);
        // if the occurrence of the merge is even, shift the origin of the grid
        if((height / dim) % 2 == 1) {
            root_grid.set_origin_coordinate(i, root_grid.origin()[i]-root_grid.lengths()[i]);
            coordinates[i]++;
        }
        root_grid.set_length(i, 2.0*root_grid.lengths()[i]);
        // get the new root cell
        if(coordinates[i] < 0) coordinates[i]--;
        coordinates[i] = coordinates[i] / 2;
        root_cell = root_grid.cell(coordinates);
        // increase height
        height++;    
    }
    // std::cout << "final height = " << height << ", final root_cell = " << root_cell << std::endl;        
    return this->increase_height(height);
}


BDDTreeSet join( const BDDTreeSet& set1, const BDDTreeSet& set2 ) {
    BDDTreeSet res(set1);
    res.adjoin(set2);
    return res;
}

BDDTreeSet intersection( const BDDTreeSet& set1, const BDDTreeSet& set2 ) {
    BDDTreeSet res(set1);
    res.restrict(set2);
    return res;
}

BDDTreeSet difference( const BDDTreeSet& set1, const BDDTreeSet& set2 ) {
    BDDTreeSet res(set1);
    res.remove(set2);
    return res;
}

void BDDTreeSet::adjoin( const BDDTreeSet& set ) {
    // raise an error if the current set is zero-dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot adjoin a BDDTreeSet to a zero-dimensional one.");
    // the two sets must have the same grid
    ARIADNE_ASSERT_MSG(this->grid() == set.grid(), "Cannot adjoin a BDDTreeSet with a different grid.");
    
    // Make a copy of set that can be modified
    BDDTreeSet set2 = set;
    // Equalize the root cell of the two sets
    _equalize_root_cell(*this, set2);
    // Since the two sets have the same root cell, their union is the logical OR of the BDDs
    this->_bdd = bdd_or(this->_bdd, set2._bdd);
    // Minimize the height of the result
    this->minimize_height();
}


void BDDTreeSet::restrict( const BDDTreeSet& set ) {
    // raise an error if the current set is zero-dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot restrict a zero-dimensional BDDTreeSet.");
    // the two sets must have the same grid
    ARIADNE_ASSERT_MSG(this->grid() == set.grid(), "Cannot intersect a BDDTreeSet with a different grid.");
    
    // Make a copy of set that can be modified
    BDDTreeSet set2 = set;
    // Equalize the root cell of the two sets
    _equalize_root_cell(*this, set2);
    // Since the two sets have the same root cell, their intersection is the logical AND of the BDDs
    this->_bdd = bdd_and(this->_bdd, set2._bdd);
    // Minimize the height of the result
    this->minimize_height();
}

void BDDTreeSet::remove( const BDDTreeSet& set ) {
    // std::cout << "BDDTreeSet::remove( ... )" << std::endl;
    // raise an error if the current set is zero-dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot remove from a zero-dimensional BDDTreeSet.");
    // the two sets must have the same grid
    ARIADNE_ASSERT_MSG(this->grid() == set.grid(), "Cannot remove a BDDTreeSet with a different grid.");
    
    // Make a copy of set that can be modified
    BDDTreeSet set2 = set;
    // Equalize the root cell of the two sets
    _equalize_root_cell(*this, set2);
    // Since the two sets have the same root cell, their difference is the logical "difference" of the BDDs
    this->_bdd = bdd_apply(this->_bdd, set2._bdd, bddop_diff);
    // Minimize the height of the result
    this->minimize_height();
}

void BDDTreeSet::adjoin_over_approximation( const Box& box, const uint subdiv ) {
    // std::cout << "adjoin_over_approximation( " << box << ", " << subdiv << ")" << std::endl;
    // raise an error if the set is zero dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot adjoin a box to a zero-dimensional BDDTreeSet.");
    // raise an error if the dimensions of box and set are different
    ARIADNE_ASSERT_MSG(this->dimension() == box.dimension(), "Cannot adjoin a box with different dimension.");
    // raise an error if the box is unbounded
    ARIADNE_ASSERT_MSG(box.bounded(), "Cannot adjoin an unbounded box.");

    // do nothing if the box is empty
    if(box.empty()) return;
    
    // First step: increase the height of the BDDTreeSet until the box is a subset of the root cell
    uint height = this->increase_height(box);
    uint dim = this->dimension();
    uint depth = height + dim*subdiv;
    // add new variables if needed
    _ensure_bdd_variables(depth);
    // call to recursive worker procedure that computes the new bdd
    // determine which dimension to split first 
    uint i = _splitting_coordinate(height, dim);
    this->_bdd = _adjoin_over_approximation(box, this->enabled_cells(), this->root_cell(), depth, i, 0);
    // minimize the result
    this->minimize_height();
}   

void BDDTreeSet::adjoin_outer_approximation( const CompactSetInterface& set, const uint subdiv ) {
    // std::cout << "BDDTreeSet::adjoin_outer_approximation(" << set << ", " << subdiv << ")" << std::endl;
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot adjoin to a zero-dimensional bdd set.");
    // the set and the bdd set must have the same dimension
    ARIADNE_ASSERT_MSG(this->dimension() == set.dimension(), "Cannot adjoin a set with different dimension.");

    Box bbox = set.bounding_box();
    // std::cout << "bounding box: " << std::flush << bbox << std::endl;
    // do nothing if the set is empty
    if(bbox.empty()) return;
    
    // First step: increase the height of the BDDTreeSet until the set is a subset of the root cell
    uint height = this->increase_height(bbox);
    uint dim = this->dimension();
    uint depth = height + dim*subdiv;
    // add new variables if needed
    _ensure_bdd_variables(depth);

    // call to recursive worker procedure that computes the new bdd
    // std::cout << "root cell after height increase: " << this->root_cell() << std::endl;
    // ARIADNE_ASSERT_MSG(this->root_cell().covers(bbox), "Error in increasing the set height: the new root cell must cover the set.");
    // determine which dimension to split first 
    uint i = _splitting_coordinate(height, dim);
    this->_bdd = _adjoin_outer_approximation(set, this->enabled_cells(), this->root_cell(), depth, i, 0);
    // minimize the result
    this->minimize_height();    
}

void BDDTreeSet::adjoin_lower_approximation( const OvertSetInterface& set, const uint height, const uint subdiv ) {
    // std::cout << "BDDTreeSet::adjoin_lower_approximation(" << set << ", " << height << ", " << subdiv << ")" << std::endl;
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot adjoin to a zero-dimensional bdd set.");
    // the set and the bdd set must have the same dimension
    ARIADNE_ASSERT_MSG(this->dimension() == set.dimension(), "Cannot adjoin a set with different dimension.");
    
    // First step: increase the height of the BDDTreeSet to height * dim
    uint dim = this->dimension();
    uint set_height = this->increase_height(height * dim);
    uint depth = set_height + dim*subdiv;
    // add new variables if needed
    _ensure_bdd_variables(depth);

    // determine which dimension to split first 
    uint i = _splitting_coordinate(set_height, dim);
    // call to worker procedure that computes the new bdd
    this->_bdd = _adjoin_lower_approximation(set, this->enabled_cells(), this->root_cell(), depth, i, 0);
    // minimize the result
    this->minimize_height();            
}

void BDDTreeSet::adjoin_lower_approximation( const OvertSetInterface& set, const Box& bounding_box, const uint subdiv ) {
    // std::cout << "BDDTreeSet::adjoin_lower_approximation(" << set << ", " << bounding_box << ", " << subdiv << ")" << std::endl;
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot adjoin to a zero-dimensional bdd set.");
    // the set and the bdd set must have the same dimension
    ARIADNE_ASSERT_MSG(this->dimension() == set.dimension(), "Cannot adjoin a set with different dimension.");
    
    // First step: increase the height of the BDDTreeSet until it covers the bounding_box
    uint height = this->increase_height(bounding_box);
    uint dim = this->dimension();
    uint depth = height + dim*subdiv;
    // add new variables if needed
    _ensure_bdd_variables(depth);

    // call to recursive worker procedure that computes the new bdd

    // determine which dimension to split first 
    uint i = _splitting_coordinate(height, dim);
    // call to worker procedure that computes the new bdd
    this->_bdd = _adjoin_lower_approximation(set, this->enabled_cells(), this->root_cell(), depth, i, 0);
    // minimize the result
    this->minimize_height();    
}

void BDDTreeSet::adjoin_lower_approximation( const LocatedSetInterface& set, const uint subdiv ) {
    this->adjoin_lower_approximation(set, set.bounding_box(), subdiv);
}

void BDDTreeSet::adjoin_inner_approximation( const OpenSetInterface& set, const uint height, const uint subdiv ) {
    // std::cout << "BDDTreeSet::adjoin_inner_approximation(" << set << ", " << height << ", " << subdiv << ")" << std::endl;
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot adjoin to a zero-dimensional bdd set.");
    // the set and the bdd set must have the same dimension
    ARIADNE_ASSERT_MSG(this->dimension() == set.dimension(), "Cannot adjoin a set with different dimension.");
    
    // First step: increase the height of the BDDTreeSet to height * dim
    uint dim = this->dimension();
    uint set_height = this->increase_height(height * dim);
    uint depth = set_height + dim*subdiv;
    // add new variables if needed
    _ensure_bdd_variables(depth);

    // call to recursive worker procedure that computes the new bdd

    // determine which dimension to split first 
    uint i = _splitting_coordinate(set_height, dim);
    // call to worker procedure that computes the new bdd
    this->_bdd = _adjoin_inner_approximation(set, this->enabled_cells(), this->root_cell(), depth, i, 0);
    // minimize the result
    this->minimize_height();    
}

void BDDTreeSet::adjoin_inner_approximation( const OpenSetInterface& set, const Box& bounding_box, const uint subdiv ) {
    // std::cout << "BDDTreeSet::adjoin_inner_approximation(" << set << ", " << bounding_box << ", " << subdiv << ")" << std::endl;
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot adjoin to a zero-dimensional bdd set.");
    // the set and the bdd set must have the same dimension
    ARIADNE_ASSERT_MSG(this->dimension() == set.dimension(), "Cannot adjoin a set with different dimension.");
    
    // First step: increase the height of the BDDTreeSet until it covers the bounding_box
    uint height = this->increase_height(bounding_box);
    uint dim = this->dimension();
    uint depth = height + dim*subdiv;
    // add new variables if needed
    _ensure_bdd_variables(depth);

    // call to recursive worker procedure that computes the new bdd

    // determine which dimension to split first 
    uint i = _splitting_coordinate(height, dim);
    // call to worker procedure that computes the new bdd
    this->_bdd = _adjoin_inner_approximation(set, this->enabled_cells(), this->root_cell(), depth, i, 0);
    // minimize the result
    this->minimize_height();    
}

void BDDTreeSet::outer_restrict(const OpenSetInterface& set) {
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot compare with a zero-dimensional BDDTreeSet.");
    ARIADNE_ASSERT_MSG(this->dimension() == set.dimension(), "Cannot compare sets with different dimensions.");

    // Do nothing if the bdd set is empty
    if(this->empty()) return;

    // determine which dimension to split first 
    uint height = this->root_cell_height();
    uint dim = this->dimension();
    uint i = _splitting_coordinate(height, dim);
    uint depth = this->depth();
    // call to worker procedure that computes the new bdd
    this->_bdd = _approximate_restrict(set, this->enabled_cells(), this->root_cell(), height + depth, i, 0, AS_OUTER);
    // minimize the result
    this->minimize_height();    
}

void BDDTreeSet::inner_restrict(const OpenSetInterface& set) {
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot compare with a zero-dimensional BDDTreeSet.");
    ARIADNE_ASSERT_MSG(this->dimension() == set.dimension(), "Cannot compare sets with different dimensions.");

    // Do nothing if the bdd set is empty
    if(this->empty()) return;

    // determine which dimension to split first 
    uint height = this->root_cell_height();
    uint dim = this->dimension();
    uint i = _splitting_coordinate(height, dim);
    uint depth = this->depth();
    // call to worker procedure that computes the new bdd
    this->_bdd = _approximate_restrict(set, this->enabled_cells(), this->root_cell(), height + depth, i, 0, AS_INNER);
    // minimize the result
    this->minimize_height();    
}

void BDDTreeSet::outer_restrict(const SetCheckerInterface& checker, const uint accuracy) {
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot compare with a zero-dimensional BDDTreeSet.");

    // Do nothing if the bdd set is empty
    if(this->empty()) return;

    // determine which dimension to split first 
    uint height = this->root_cell_height();
    uint dim = this->dimension();
    uint i = _splitting_coordinate(height, dim);
    // call to worker procedure that computes the new bdd
    this->_bdd = _approximate_restrict(checker, this->enabled_cells(), this->root_cell(), 
                                       height + dim*accuracy, i, 0, AS_OUTER);
    // minimize the result
    this->minimize_height();    
}

void BDDTreeSet::inner_restrict(const SetCheckerInterface& checker, const uint accuracy) {
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot compare with a zero-dimensional BDDTreeSet.");

    // Do nothing if the bdd set is empty
    if(this->empty()) return;

    // determine which dimension to split first 
    uint height = this->root_cell_height();
    uint dim = this->dimension();
    uint i = _splitting_coordinate(height, dim);
    // call to worker procedure that computes the new bdd
    this->_bdd = _approximate_restrict(checker, this->enabled_cells(), this->root_cell(), 
                                       height + dim*accuracy, i, 0, AS_INNER);
    // minimize the result
    this->minimize_height();    
}

BDDTreeSet::const_iterator BDDTreeSet::begin() const {
    return BDDTreeSet::const_iterator(*this);
}

BDDTreeSet::const_iterator BDDTreeSet::end() const {
    return BDDTreeSet::const_iterator();
}

BDDTreeSet::operator ListSet<Box>() const {
    ARIADNE_ASSERT_MSG(this->dimension() != 0,
        "Cannot convert a zero-dimensional BDDTreeSet to a list of boxes.");
        
    ListSet<Box> result(this->dimension());

    for (BDDTreeSet::const_iterator it = this->begin(), end = this->end(); it != end; it++ ) {
        result.push_back((*it));
    }

    return result;
}

void BDDTreeSet::draw(CanvasInterface& canvas) const {
    for(BDDTreeSet::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        iter->draw(canvas);
    }
}


std::ostream& BDDTreeSet::write(std::ostream& os) const {
    return os << (*this);
}


/*
    void import_from_file(const char*& filename);

    void export_to_file(const char*& filename);

*/


/********************************** BDDTreeConstIterator **********************************/
int operator==(const PathElement& a, const PathElement& b) {
	return (a.obdd == b.obdd) && (a.status == b.status) && (a.cell == b.cell) &&
	       (a.root_var == a.root_var) && (a.split_coordinate == b.split_coordinate);
}


BDDTreeConstIterator::BDDTreeConstIterator() 
    : _path()
    , _mince_depth(-1)
{
}

BDDTreeConstIterator::BDDTreeConstIterator( const BDDTreeSet& set ) 
{
    ARIADNE_ASSERT_MSG(set.dimension() != 0, "Cannot create an iterator for a zero-dimensional BDDTreeSet.");
    PathElement root;
    root.cell = set.root_cell();
    root.status = PE_NEW;
    root.obdd = set.enabled_cells();
    root.root_var = 0;
    // compute first splitting coordinate
    uint height = set.root_cell_height();
    root.split_coordinate = 0;
    if(height != 0) {
        uint dim = set.dimension();
        root.split_coordinate = (dim - 1) - ((height-1) % dim);
    }
    this->_path.push_back(root);
    // compute mince_depth
    if(set.mince_depth() >= 0) {
        this->_mince_depth = set.root_cell_height() + set.mince_depth();
    } else {
        this->_mince_depth = -1;    
    }
    _compute_next_cell(this->_path, this->_mince_depth);
}

BDDTreeConstIterator::BDDTreeConstIterator( const BDDTreeConstIterator& iter )
    : _path(iter._path)
    , _mince_depth(iter._mince_depth)
{
}

void BDDTreeConstIterator::increment() {
    // increment only if the iterator is a valid one
    if(!this->_path.empty()) {
        // compute initial mince_depth
        int mince_depth = this->_mince_depth - (this->_path.size() - 1);
        _compute_next_cell(this->_path, mince_depth);
    }
}

bool BDDTreeConstIterator::equal( BDDTreeConstIterator const & other ) const {
    // two iterators are equal if they are both invalid, or if they have the same path and mince_depth
    return (this->_path == other._path) && ((this->_mince_depth == other._mince_depth) || this->_path.empty());
}


Box const& BDDTreeConstIterator::dereference() const {
    ARIADNE_ASSERT_MSG(!this->_path.empty(), "Cannot dereference an invalid BDDTreeConstIterator.");
    const PathElement& tail = this->_path.back();
    return tail.cell;
}

/********************************** Stream operators **********************************/

std::ostream& operator<<(std::ostream& os, const BDDTreeSet& set) {
    return os << "BDDTreeSet( Grid: " << set.grid() << 
                 ", root cell height: " << set.root_cell_height() << 
                 ", root cell coordinates: " << set.root_cell_coordinates() <<
                 ", bdd: " << set.enabled_cells() <<" )";
}

std::ostream& operator<<(std::ostream& os, const PathElement& pe) {
    return os << "( " << pe.cell << 
                 ", " << pe.obdd << 
                 ", " << pe.status <<
                 ", " << pe.split_coordinate << 
                 ", " << pe.root_var <<" )";
}

/********************************** Comparison and filtering operators **********************************/

tribool disjoint(const ConstraintSet& cons_set, const BDDTreeSet& bdd_set) {
    ARIADNE_ASSERT_MSG(bdd_set.dimension() != 0, "Cannot compare with a zero-dimensional BDDTreeSet.");
    ARIADNE_ASSERT_MSG(bdd_set.dimension() == cons_set.dimension(), "Cannot compare sets with different dimensions.");

    if(bdd_set.empty()) return true;
    if(cons_set.unconstrained()) return false;
    if(cons_set.empty()) return true;
        
    // compute the first splitting coordinate
    uint dim = bdd_set.dimension();
    uint height = bdd_set.root_cell_height();
    uint splitting_coordinate;
    if(height == 0) {
        splitting_coordinate = 0;
    } else {
        splitting_coordinate = (dim - 1) - ((height-1) % dim);
    }
    return _disjoint(bdd_set.root_cell(), bdd_set.enabled_cells(), 0, splitting_coordinate, cons_set);            
}

tribool overlaps(const ConstraintSet& cons_set, const BDDTreeSet& bdd_set) {
    return !disjoint(cons_set, bdd_set);            
}

tribool covers(const ConstraintSet& cons_set, const BDDTreeSet& bdd_set) {
    // std::cout << "covers(" << cons_set << ", " << bdd_set << ")" << std::endl;
    ARIADNE_ASSERT_MSG(bdd_set.dimension() != 0, "Cannot compare with a zero-dimensional BDDTreeSet.");
    ARIADNE_ASSERT_MSG(bdd_set.dimension() == cons_set.dimension(), "Cannot compare sets with different dimensions.");

    if(bdd_set.empty()) return true;
    if(cons_set.unconstrained()) return true;
    if(cons_set.empty()) return false;

    // compute the first splitting coordinate
    uint dim = bdd_set.dimension();
    uint height = bdd_set.root_cell_height();
    uint splitting_coordinate;
    if(height == 0) {
        splitting_coordinate = 0;
    } else {
        splitting_coordinate = (dim - 1) - ((height-1) % dim);
    }
    return _subset(bdd_set.root_cell(), bdd_set.enabled_cells(), 0, splitting_coordinate, cons_set);        
}

BDDTreeSet outer_intersection(const BDDTreeSet& bdd_set, const ConstraintSet& cons_set) {
    // make a copy of bdd_set
    BDDTreeSet result = bdd_set;
    if (!cons_set.unconstrained())
    	result.outer_restrict(cons_set);
    return result;
}

BDDTreeSet inner_intersection(const BDDTreeSet& bdd_set, const ConstraintSet& cons_set) {
    // make a copy of bdd_set
    BDDTreeSet result = bdd_set;
    if (!cons_set.unconstrained())
    	result.inner_restrict(cons_set);
    return result;
}

// Box eps_codomain(const BDDTreeSet& bdd_set, const Vector<Float> eps, const VectorFunction& func) {
//     ARIADNE_NOT_IMPLEMENTED;    
// }

BDDTreeSet project_down(const BDDTreeSet& bdd_set, const Vector<uint>& indices) {
    // std::cout << "project_down(" << bdd_set << ", " << indices << ")" << std::endl;
    // project down the grid
    Grid new_grid = project_down(bdd_set.grid(), indices);
    // std::cout << "new_grid = " << new_grid << std::endl;
    // To simplify the procedure, increase the height of the bdd_set so that the first splitting coordinate is 0
    BDDTreeSet set = bdd_set;
    uint height = set.root_cell_height();
    uint old_dim = set.dimension();
    uint split = _splitting_coordinate(height, old_dim);
    if(split > 0) {
        set.increase_height(height + split);
    }
    // std::cout << "height = " << set.root_cell_height() << std::endl;
    ARIADNE_ASSERT(_splitting_coordinate(set.root_cell_height(), set.dimension()) == 0);
    // compute the new height and root cell coordinates
    uint new_height = (set.root_cell_height() / old_dim) * indices.size();
    array<int> new_coordinates(new_grid.dimension(), 0);
    array<int> old_coordinates = bdd_set.root_cell_coordinates();
    for(uint i = 0; i != indices.size(); ++i) {
        new_coordinates[i] = old_coordinates[indices[i]];
    }
    // compute the new mince depth
    uint new_mince_depth = _project_mince_depth(set.mince_depth(), old_dim, indices);
    // remove, duplicate and rename variables
    bdd new_bdd = _project_down_variables(set.enabled_cells(), old_dim, indices, set.root_cell_height() + set.mince_depth());
    // create the new BDDTreeSet
    set = BDDTreeSet(new_grid, new_height, new_coordinates, new_bdd, new_mince_depth);
    set.minimize_height();
    return set;
}

// tribool covers(const BDDTreeSet& covering_set, const BDDTreeSet& covered_set, const Vector<Float>& eps) {
//     ARIADNE_NOT_IMPLEMENTED;    
// }
// 
// tribool inside(const BDDTreeSet& covered_set, const BDDTreeSet& covering_set, const Vector<Float>& eps, int accuracy) {
//     ARIADNE_NOT_IMPLEMENTED;    
// }

} // namespace Ariadne

