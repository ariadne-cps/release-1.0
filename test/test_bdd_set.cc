/***************************************************************************
 *            test_grid_set.cc
 *
 *
 *  Copyright  2008  Ivan S. Zapreev, Pieter Collins
 *            ivan.zapreev@gmail.com, pieter.collins@cwi.nl
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "config.h"

#include "bdd_set.h"
#include "logging.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

void test_constructors() {
    ARIADNE_PRINT_TEST_COMMENT("Testing constructors.");
    BDDTreeSet set1;
    ARIADNE_TEST_EQUAL(set1.dimension(),0);

    BDDTreeSet set2(3);
    ARIADNE_TEST_EQUAL(set2.grid(),Grid(3));
    ARIADNE_TEST_EQUAL(set2.root_cell_height(),0);
    ARIADNE_TEST_EQUAL(set2.root_cell_coordinates(),array<int>(3, 0,0,0));
    ARIADNE_TEST_EQUAL(set2.enabled_cells(),bddfalse);

    Grid grid(4, 1.25);
    BDDTreeSet set3(grid, true);
    ARIADNE_TEST_EQUAL(set3.grid(),grid);
    ARIADNE_TEST_EQUAL(set3.root_cell_height(),0);
    ARIADNE_TEST_EQUAL(set3.root_cell_coordinates(),array<int>(4, 0,0,0,0));
    ARIADNE_TEST_EQUAL(set3.enabled_cells(),bddtrue);
    
    // Test copy constructor
    set2 = set3;
    ARIADNE_TEST_EQUAL(set2,set3);

    // Test cloning operator
    BDDTreeSet* pset = set2.clone();
    ARIADNE_TEST_EQUAL(set2, *pset);

    // Test construction from a bdd
    bdd enabled_cells = bdd_ithvar(2);
    BDDTreeSet set4(grid, 1, array<int>(4, 1,0,0,0), enabled_cells);
    ARIADNE_TEST_EQUAL(set4.grid(),grid);
    ARIADNE_TEST_EQUAL(set4.root_cell_height(),1);
    ARIADNE_TEST_EQUAL(set4.root_cell_coordinates(),array<int>(4, 1,0,0,0));
    ARIADNE_TEST_EQUAL(set4.enabled_cells(),enabled_cells);
    
    // If the dimension of root_cell_coordinates and grid differs an exception must be thrown.
    ARIADNE_TEST_FAIL(set4 = BDDTreeSet(grid, 1, array<int>(3, 0,0,0), enabled_cells));
}

void test_properties() {
    ARIADNE_PRINT_TEST_COMMENT("Testing properties.");
    // Test empty
    BDDTreeSet set0;
    // The zero-dimensional set is always empty
    ARIADNE_TEST_ASSERT(set0.empty());
    BDDTreeSet set1(2, false);
    ARIADNE_TEST_ASSERT(set1.empty());
    bdd enabled_cells = bdd_nithvar(0) & (bdd_ithvar(2) | bdd_ithvar(3));
    BDDTreeSet set2(Grid(3), 4, array<int>(3, 1,0,0), enabled_cells);
    ARIADNE_TEST_ASSERT(!set2.empty());
    
    // Test size
    ARIADNE_TEST_EQUAL(set1.size(), 0);
    ARIADNE_TEST_EQUAL(set2.size(), 4);
    
    // Measure is not implemented yet
    ARIADNE_TEST_FAIL(set2.measure());
    
    // Test root_cell
    // raise an error if the set is zero dimensional
    ARIADNE_TEST_FAIL(set0.root_cell());
    ARIADNE_TEST_EQUAL(set1.root_cell(), Box(2, 0.0,1.0, 0.0,1.0));
    set0 = BDDTreeSet(Grid(2), 3, array<int>(2, 0,0), enabled_cells);
    ARIADNE_TEST_EQUAL(set0.root_cell(), Box(2, 0.0,2.0, -2.0,2.0));
    ARIADNE_TEST_EQUAL(set2.root_cell(), Box(3, 2.0,4.0, 0.0,2.0, -2.0,2.0));   
    
    // Test bounding box
    // Note that for set0 and set2 the results are different from the previous call to root_cell
    // because bounding_box() minimizes the root cell height to compute the smallest box
    ARIADNE_TEST_EQUAL(set1.bounding_box(), Box(2, 0.0,1.0, 0.0,1.0));
    ARIADNE_TEST_EQUAL(set0.bounding_box(), Box(2, 0.0,2.0, -2.0,0.0));       
    ARIADNE_TEST_EQUAL(set2.bounding_box(), Box(3, 2.0,4.0, 0.0,2.0, -2.0,0.0));       
}

void test_predicates() {
    ARIADNE_PRINT_TEST_COMMENT("Testing predicates.");
    // Test subset and superset
    BDDTreeSet set0;
    BDDTreeSet set1(Grid(2), false);
    // testing w.r.t. a zero-dimensional set should raise an error
    ARIADNE_TEST_FAIL(subset(set0, set1));
    ARIADNE_TEST_FAIL(superset(set0, set1));
    // testing sets with different grids should raise an error
    BDDTreeSet set2(Grid(2, 1.25), true);
    ARIADNE_TEST_FAIL(subset(set1, set2));
    ARIADNE_TEST_FAIL(superset(set1, set2));
    
    // check test with an empty set
    bdd enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) | bdd_ithvar(4));    
    set2 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), enabled_cells);
    ARIADNE_TEST_ASSERT(subset(set1,set2));
    ARIADNE_TEST_ASSERT(!superset(set1,set2));
    
    // Sets with different root cell coordinates are disjoint
    set1 = BDDTreeSet(Grid(2), true);
    ARIADNE_TEST_ASSERT(!subset(set1,set2));
    ARIADNE_TEST_ASSERT(!superset(set1,set2));
    
    // Sets with the same root coordinates
    set1 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), bdd_ithvar(0));
    ARIADNE_TEST_ASSERT(subset(set2,set1));
    ARIADNE_TEST_ASSERT(!superset(set2,set1));
    
    // Subset is a reflexive relation
    ARIADNE_TEST_ASSERT(subset(set2,set2));
    ARIADNE_TEST_ASSERT(superset(set2,set2));    
    
    // Test disjoint and overlap
    // testing w.r.t. a zero-dimensional set should raise an error
    ARIADNE_TEST_FAIL(disjoint(set0, set1));
    ARIADNE_TEST_FAIL(overlap(set0, set1));
     // testing sets with different grids should raise an error
    set2 = BDDTreeSet(Grid(2, 1.25), true);
    ARIADNE_TEST_FAIL(subset(set1, set2));
    ARIADNE_TEST_FAIL(superset(set1, set2));
    // check test with an empty set
    set1 = BDDTreeSet(Grid(2), false);
    enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) | bdd_ithvar(4));    
    set2 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), enabled_cells);
    ARIADNE_TEST_ASSERT(disjoint(set1,set2));
    ARIADNE_TEST_ASSERT(!overlap(set1,set2));
    ARIADNE_TEST_ASSERT(!disjoint(set1,set1));
    ARIADNE_TEST_ASSERT(overlap(set1,set1));
    
    // Sets with different root cell coordinates are disjoint
    set1 = BDDTreeSet(Grid(2), true);
    ARIADNE_TEST_ASSERT(!subset(set1,set2));
    ARIADNE_TEST_ASSERT(!superset(set1,set2));

    // Sets with the same root coordinates
    set1 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), bdd_ithvar(0));
    ARIADNE_TEST_ASSERT(!disjoint(set1,set2));
    ARIADNE_TEST_ASSERT(overlap(set1,set2));
    set1 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), bdd_ithvar(0) & bdd_ithvar(1));
    ARIADNE_TEST_ASSERT(disjoint(set1,set2));
    ARIADNE_TEST_ASSERT(!overlap(set1,set2));
    
    // A set always overlaps itself
    ARIADNE_TEST_ASSERT(!disjoint(set2,set2));
    ARIADNE_TEST_ASSERT(overlap(set2,set2));    
    
    // Test subset inclusion with a Box
    // testing w.r.t. a zero-dimensional set or box should raise an error
    ARIADNE_TEST_FAIL(set0.subset(Box(2, 0.0,1.0, 0.0,1.0)));
    ARIADNE_TEST_FAIL(set1.subset(Box(0)));
     // testing w.rt. to a box with different dimension should raise an error
    ARIADNE_TEST_FAIL(set1.subset(Box(3)));
    // check test with an empty set or empty box
    set1 = BDDTreeSet(Grid(2), false);
    ARIADNE_TEST_ASSERT(definitely(set1.subset(Box(2, 0.0,1.0, 0.0,1.0))));
    Box ebx = Box::empty_box(2);
    ARIADNE_TEST_ASSERT(!possibly(set2.subset(ebx)));
    // test with a general set
    enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) & bdd_ithvar(4));    
    set2 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), enabled_cells);
    Box bx1(2, 2.25,5.0, 4.25,6.5);
    Box bx2(2, 2.25,5.0, 4.75,6.5);
    ARIADNE_TEST_ASSERT(definitely(set2.subset(bx1)));
    ARIADNE_TEST_ASSERT(!possibly(set2.subset(bx2)));
       
}

void test_operations() {
    ARIADNE_PRINT_TEST_COMMENT("Testing operations.");
    // Test minimize_height
    BDDTreeSet set0;
    // raise an error if the set is zero dimensional
    ARIADNE_TEST_FAIL(set0.minimize_height());
    BDDTreeSet set1(Grid(3), true);
    BDDTreeSet set2 = set1;
    // No changes if the height is zero.
    ARIADNE_TEST_CHECK(set1.minimize_height(), 0);
    ARIADNE_TEST_EQUAL(set1, set2);

    bdd enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) | bdd_ithvar(4));    
    set1 = BDDTreeSet(Grid(3), 4, array<int>(3, 1,1,1), enabled_cells);
    set2 = set1;
    ARIADNE_TEST_CHECK(set1.minimize_height(), 2);
    ARIADNE_TEST_ASSERT(set1 != set2);
    ARIADNE_TEST_EQUAL(set1.root_cell_height(), 2);
    ARIADNE_TEST_EQUAL(set1.root_cell_coordinates(), array<int>(3, 2,1,2));
    enabled_cells = bdd_ithvar(1) | bdd_ithvar(2);
    ARIADNE_TEST_EQUAL(set1.enabled_cells(), enabled_cells);
    ARIADNE_TEST_EQUAL(set1.root_cell(), Box(3, 2.0,3.0, 2.0,4.0, 4.0,6.0));
    
    // Test increase height
    // raise an error if the set is zero dimensional
    ARIADNE_TEST_FAIL(set0.increase_height(5));
    BDDTreeSet set3 = BDDTreeSet(2, true);
    ARIADNE_TEST_CHECK(set3.increase_height(5), 5);
    ARIADNE_TEST_EQUAL(set3.grid(),Grid(2));
    ARIADNE_TEST_EQUAL(set3.root_cell_height(), 5);
    ARIADNE_TEST_EQUAL(set3.root_cell_coordinates(), array<int>(2, 0,0));
    ARIADNE_TEST_EQUAL(set3.root_cell(), Box(2, -2.0,2.0, -2.0,6.0));
    enabled_cells = bdd_nithvar(0) & bdd_ithvar(1) & bdd_ithvar(2) & bdd_nithvar(3) & bdd_nithvar(4);
    ARIADNE_TEST_EQUAL(set3.enabled_cells(), enabled_cells);
    
    ARIADNE_TEST_CHECK(set1.increase_height(4), 4);
    ARIADNE_TEST_EQUAL(set1, set2);
    // if the new height is smaller than the current one the set is not changed
    ARIADNE_TEST_EQUAL(set1.increase_height(1), 4);
    ARIADNE_TEST_EQUAL(set1, set2);
        
}

int main() {

    test_constructors();
    test_properties();
    test_predicates();
    test_operations();
    
    return ARIADNE_TEST_FAILURES;
}

