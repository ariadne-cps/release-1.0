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
    //           << ", " << new_height << ")" << std::endl;
    // if the new height is equal or smaller than the current one do nothing
    if(root_cell_height >= new_height) return root_cell_height;

    // Increase the root cell height by one and recursively call _minimize_height
    // compute the new root_cell_coordinate
    // determine which dimension to merge 
    uint dim = root_cell_coordinates.size();
    uint i = (dim - 1) - (root_cell_height % dim);
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
    ARIADNE_ASSERT_MSG(oldvars[0] + n >= 0, "Wrong shift index: negative variable number.")
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
// is a subset of Box. 
bool _subset(const Box& root_cell, const bdd& enabled_cells, int root_var, uint splitting_coordinate, const Box& box) {
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
// is a disjoint from Box. 
bool _disjoint(const Box& root_cell, const bdd& enabled_cells, int root_var, uint splitting_coordinate, const Box& box) {
    // std::cout << "_disjoint(" << root_cell << ", " << enabled_cells << ", " 
    //           << root_var << ", " << splitting_coordinate << ", " << box << ")" << std::endl;
    // if the set is empty, return true
    if(enabled_cells == bddfalse) return true;
    // if the root cell is disjoint from box, return true
    if(root_cell.disjoint(box)) return true;
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
    if(!_disjoint(subcells.first, left_branch, root_var+1, (splitting_coordinate+1) % dim, box))
        return false;
    return _disjoint(subcells.second, right_branch, root_var+1, (splitting_coordinate+1) % dim, box);
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



// Function that increments a BDDTreeSet iterator
void _compute_next_cell(std::vector< PathElement >& path) {
    // std::cout << "_compute_next_cell( ... )" << std::endl;
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
        _compute_next_cell(path);
        return;
    }
    // if the status is LEFT, explore the right child
    if(tail.status == PE_LEFT) {
        // std::cout << "  status is PE_LEFT, explore right child." << std::endl;        
        // set the status to RIGHT
        tail.status = PE_RIGHT;
        path.pop_back();
        path.push_back(tail);
        // get the variable labeling the bdd
        int var = bdd_var(tail.obdd);
        // if the var labeling the bdd correspond to che current level, get the right child
        if(var == tail.root_var) {
            tail.obdd = bdd_high(tail.obdd);
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
        _compute_next_cell(path);
        return;
    }
    // if the status is NEW, we are exploring the node for the first time
    // if the bdd is the constant true, we are in a cell: fix the flag to true and return
    // std::cout << "  status is PE_NEW ";        
    if(tail.obdd == bddtrue) {
        // std::cout << "and the current bdd is TRUE, exiting" << std::endl;
        // set the status to RIGHT
        tail.status = PE_RIGHT;
        path.pop_back();
        path.push_back(tail);
        return;
    }
    // if the bdd is the constant false, backtrack
    if(tail.obdd == bddfalse) {
        // std::cout << "and the current bdd is FALSE, backtracking" << std::endl;
        path.pop_back();
        _compute_next_cell(path);
        return;
    }
    // the bdd is different from true and false, explore the left child
    // set the status to left
    // std::cout << "and the current bdd is neither TRUE nor FALSE, go on" << std::endl;
    tail.status = PE_LEFT;
    path.pop_back();
    path.push_back(tail);
    // get the variable labeling the bdd
    int var = bdd_var(tail.obdd);
    // if the var labeling the bdd correspond to che current level, get the right child
    if(var == tail.root_var) {
        tail.obdd = bdd_low(tail.obdd);
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
    _compute_next_cell(path);
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


Box BDDTreeSet::bounding_box() const {
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
    uint splitting_coordinate;
    if(height == 0) {
        splitting_coordinate = 0;
    } else {
        splitting_coordinate = (dim - 1) - ((height-1) % dim);
    }
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
    uint splitting_coordinate;
    if(height == 0) {
        splitting_coordinate = 0;
    } else {
        splitting_coordinate = (dim - 1) - ((height-1) % dim);
    }
    return _disjoint(this->root_cell(), this->enabled_cells(), 0, splitting_coordinate, box);        
}

tribool BDDTreeSet::overlaps( const Box& box ) const {
    return !this->disjoint(box);        
}

void BDDTreeSet::clear( ) {
    this->_root_cell_height = 0;
    this->_root_cell_coordinates.fill(0);
    this->_bdd = bddfalse;
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

    Box root_cell = this->root_cell();
    uint height = this->root_cell_height();
    uint dim = this->dimension();
    int i = 0;
    while(!definitely(box.subset(root_cell))) {
        // std::cout << "height = " << height << ", root_cell = " << root_cell << std::endl;        
        // determine which dimension to merge 
        i = (dim - 1) - (height % dim);
        // if the occurrence of the merge is even, enlarge the cell to the right
        if((height / dim) % 2 == 1) {
            root_cell[i] = root_cell[i] - Interval(0.0, root_cell[i].width());
        } else {    // otherwise, to the left
            root_cell[i] = root_cell[i] + Interval(0.0, root_cell[i].width());
        }
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
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot adjoin a BDDTeeSet to a zero-dimensional one.");
    // the two sets must have the same grid
    ARIADNE_ASSERT_MSG(this->grid() == set.grid(), "Cannot adjoin a BDDTeeSet with a different grid.");
    
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
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot restrict a zero-dimensional BDDTeeSet.");
    // the two sets must have the same grid
    ARIADNE_ASSERT_MSG(this->grid() == set.grid(), "Cannot intersect a BDDTeeSet with a different grid.");
    
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
    // raise an error if the current set is zero-dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot remove from a zero-dimensional BDDTeeSet.");
    // the two sets must have the same grid
    ARIADNE_ASSERT_MSG(this->grid() == set.grid(), "Cannot remove a BDDTeeSet with a different grid.");
    
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
    // raise an error if the set is zero dimensional
    ARIADNE_ASSERT_MSG(this->dimension() != 0, "Cannot adjoin a box to a zero-dimensional BDDTreeSet.");
    // raise an error if the dimensions of box and set are different
    ARIADNE_ASSERT_MSG(this->dimension() == box.dimension(), "Cannot adjoin a box with different dimension.");
    // raise an error if the box is unbounded
    ARIADNE_ASSERT_MSG(box.bounded(), "Cannot adjoin an unbounded box.");

    // do nothing if the box is empty
    if(box.empty()) return;
    
    // First step: increase the height of the BDDTreeSet until the box is a subset of the root cell
    this->increase_height(box);
    
    // recursive call to worker procedure that computes the new bdd
    uint height = this->root_cell_height();
    uint dim = this->dimension();
    // determine which dimension to split first 
    uint i = 0;
    if(height > 0) i = (dim - 1) - ((height-1) % dim);    
    this->_bdd = _adjoin_over_approximation(box, this->enabled_cells(), this->root_cell(), height + dim*subdiv, i, 0);
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
    this->increase_height(bbox);

    // recursive call to worker procedure that computes the new bdd
    uint height = this->root_cell_height();
    uint dim = this->dimension();
    // determine which dimension to split first 
    uint i = 0;
    if(height > 0) i = (dim - 1) - ((height-1) % dim);    
    this->_bdd = _adjoin_outer_approximation(set, this->enabled_cells(), this->root_cell(), height + dim*subdiv, i, 0);
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

    // determine which dimension to split first 
    uint i = 0;
    if(set_height > 0) i = (dim - 1) - ((set_height-1) % dim);    
    // call to worker procedure that computes the new bdd
    this->_bdd = _adjoin_lower_approximation(set, this->enabled_cells(), this->root_cell(), set_height + dim*subdiv, i, 0);
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

    // determine which dimension to split first 
    uint i = 0;
    if(height > 0) i = (dim - 1) - ((height-1) % dim);    
    // call to worker procedure that computes the new bdd
    this->_bdd = _adjoin_lower_approximation(set, this->enabled_cells(), this->root_cell(), height + dim*subdiv, i, 0);
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

    // determine which dimension to split first 
    uint i = 0;
    if(set_height > 0) i = (dim - 1) - ((set_height-1) % dim);    
    // call to worker procedure that computes the new bdd
    this->_bdd = _adjoin_inner_approximation(set, this->enabled_cells(), this->root_cell(), set_height + dim*subdiv, i, 0);
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

    // determine which dimension to split first 
    uint i = 0;
    if(height > 0) i = (dim - 1) - ((height-1) % dim);    
    // call to worker procedure that computes the new bdd
    this->_bdd = _adjoin_inner_approximation(set, this->enabled_cells(), this->root_cell(), height + dim*subdiv, i, 0);
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
    _compute_next_cell(this->_path);
}

BDDTreeConstIterator::BDDTreeConstIterator( const BDDTreeConstIterator& iter )
    : _path(iter._path)
{
}

void BDDTreeConstIterator::increment() {
    // increment only if the iterator is a valid one
    if(!this->_path.empty()) {
        _compute_next_cell(this->_path);
    }
}

bool BDDTreeConstIterator::equal( BDDTreeConstIterator const & other ) const {
    return (this->_path == other._path);
}

Box const& BDDTreeConstIterator::dereference() const {
    ARIADNE_ASSERT_MSG(!this->_path.empty(), "Cannot dereference and invalid BDDTreeConstIterator.");
    PathElement tail = this->_path.back();
    Box *pCell = new Box(tail.cell);
    return (*pCell);
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


} // namespace Ariadne

