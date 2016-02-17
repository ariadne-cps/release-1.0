/***************************************************************************
 *            test_reachability_restriction.cc
 *
 *  Copyright  2011  Luca Geretti
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

#include "test.h"
#include "ariadne.h"

using namespace std;
using namespace Ariadne;


class TestReachabilityRestriction {
  public:
    void test();
  private:
    void test_accessors();
    void test_refine();
    void test_apply_to();
    void test_restricts();
    void test_update();
    void test_copy();
    void test_forward_jump_set();
    void test_backward_jump_set();
    void test_projection();

    HybridAutomaton _get_system();
    ReachabilityRestriction _get_reference_restriction();
};


HybridAutomaton
TestReachabilityRestriction::_get_system()
{
	HybridAutomaton sys;

	DiscreteLocation q1(1);
	DiscreteLocation q2(2);

	DiscreteEvent q1q2("q1q2");

	RealVariable x("x");
	RealVariable y("y");
	List<RealVariable> varlist;
	varlist.append(x);
	varlist.append(y);

	RealExpression idx = x;
	RealExpression zero = 0.0;
	RealExpression one = 1.0;

	RealExpression x_q1 = 1.0;
	RealExpression y_q1 = 1.0;

	RealExpression x_q2 = 0.0;
	RealExpression y_q2 = 0.0;

	List<RealExpression> exprlist;
	exprlist.append(x_q1);
	exprlist.append(y_q1);
	VectorFunction q1_d(exprlist, varlist);
	exprlist[0] = x_q2;
	exprlist[1] = y_q2;
	VectorFunction q2_d(exprlist, varlist);

	exprlist[0] = zero;
	exprlist[1] = zero;
	VectorFunction reset_zero(exprlist, varlist);
	exprlist[0] = x+2;
	exprlist[1] = y+2;
	VectorFunction reset_plus_two(exprlist, varlist);

	RealExpression guard_q1q2_expr = x+y-6.0;
	ScalarFunction guard_q1q2(guard_q1q2_expr, varlist);

	sys.new_mode(q1,q1_d);
	sys.new_mode(q2,q2_d);

	sys.new_forced_transition(q1q2,q1,q2,reset_plus_two,guard_q1q2);

	return sys;
}


ReachabilityRestriction
TestReachabilityRestriction::_get_reference_restriction()
{
	DiscreteLocation q1(1);
	DiscreteLocation q2(2);

	HybridBoxes domain;
	domain.insert(make_pair(q1,Box(2,-1.0,1.0,-2.0,3.0)));
	domain.insert(make_pair(q2,Box(3,0.5,2.0,-1.0,0.0,1.0,1.5)));

	HybridGrid grid;
	grid[q1] = Grid(Vector<Float>(2,2.0,3.0));
	grid[q2] = Grid(Vector<Float>(3,1.0,0.5,1.0));

	int accuracy = 1;

	ReachabilityRestriction rr(domain,grid,accuracy);

	return rr;
}


void 
TestReachabilityRestriction::test() 
{
    ARIADNE_TEST_CALL(test_accessors());
    ARIADNE_TEST_CALL(test_refine());
    ARIADNE_TEST_CALL(test_apply_to());
    ARIADNE_TEST_CALL(test_restricts());
    ARIADNE_TEST_CALL(test_update());
    ARIADNE_TEST_CALL(test_copy());
    ARIADNE_TEST_CALL(test_forward_jump_set());
    ARIADNE_TEST_CALL(test_backward_jump_set());
    ARIADNE_TEST_CALL(test_projection());
}


void
TestReachabilityRestriction::test_accessors() {

	ReachabilityRestriction rr1 = _get_reference_restriction();

	DiscreteLocation q1(1);
	DiscreteLocation q2(2);
	DiscreteLocation q3(3);

	HybridBoxes domain;
	domain.insert(make_pair(q1,Box(2,-1.0,1.0,-2.0,3.0)));
	domain.insert(make_pair(q2,Box(3,0.5,2.0,-1.0,0.0,1.0,1.5)));

	ARIADNE_PRINT_TEST_CASE_TITLE("Location presence");

	ARIADNE_TEST_ASSERT(rr1.has_location(q1));
	ARIADNE_TEST_ASSERT(rr1.has_location(q2));
	ARIADNE_TEST_ASSERT(!rr1.has_location(q3));

	ARIADNE_PRINT_TEST_CASE_TITLE("Discretisation absence at startup");

	ARIADNE_TEST_ASSERT(!rr1.has_discretised(q1));
	ARIADNE_TEST_ASSERT(!rr1.has_discretised(q2));
    ARIADNE_TEST_FAIL(rr1.has_discretised(q3));

    ARIADNE_PRINT_TEST_CASE_TITLE("Bounding box call effect at startup");

    ARIADNE_TEST_ASSERT(rr1.bounding_box().superset(domain));
	ARIADNE_TEST_ASSERT(!rr1.has_discretised(q1));
	ARIADNE_TEST_ASSERT(!rr1.has_discretised(q2));
}


void
TestReachabilityRestriction::test_refine() {

	ReachabilityRestriction rr1 = _get_reference_restriction();

	DiscreteLocation q1(1);
	DiscreteLocation q2(2);

    HybridBoxes prerefine_boxes = rr1.bounding_box();
    rr1.refine_at(2);

	ARIADNE_PRINT_TEST_CASE_TITLE("Accuracy change");

    ARIADNE_TEST_EQUAL(rr1.accuracy(),2);

    ARIADNE_PRINT_TEST_CASE_TITLE("Discretisation invariance");

    ARIADNE_TEST_ASSERT(!rr1.has_discretised(q1));
    ARIADNE_TEST_ASSERT(!rr1.has_discretised(q2));

    ARIADNE_PRINT_TEST_CASE_TITLE("Bounding boxes change");

    ARIADNE_TEST_ASSERT(prerefine_boxes != rr1.bounding_box());
    ARIADNE_TEST_ASSERT(prerefine_boxes.superset(rr1.bounding_box()));
}


void
TestReachabilityRestriction::test_apply_to() {

	ReachabilityRestriction rr1 = _get_reference_restriction();

	DiscreteLocation q1(1);
	DiscreteLocation q2(2);

	HybridGrid grid = rr1.grid();
	int accuracy = rr1.accuracy();

	ARIADNE_PRINT_TEST_CASE_TITLE("Apply to a HybridDenotableSet that would not be restricted thanks to inclusion test");

	HybridDenotableSet set_to_restrict1(grid);
	HybridBoxes set_to_restrict1_boxToAdjoin;
	set_to_restrict1_boxToAdjoin.insert(make_pair(q1,Box(2,-0.4,0.5,-1.0,0.0)));
	set_to_restrict1_boxToAdjoin.insert(make_pair(q2,Box(3,0.9,1.0,-0.4,-0.2,1.2,1.3)));
	set_to_restrict1.adjoin_outer_approximation(set_to_restrict1_boxToAdjoin,accuracy);
	HybridDenotableSet restricted_set1 = set_to_restrict1;
	rr1.apply_to(restricted_set1);
	ARIADNE_TEST_ASSERT(restricted_set1 == set_to_restrict1);
	ARIADNE_TEST_ASSERT(!rr1.has_discretised(q1));
	ARIADNE_TEST_ASSERT(!rr1.has_discretised(q2));

	ARIADNE_PRINT_TEST_CASE_TITLE("Apply to a HybridDenotableSet that would not be restricted, but on q1 the inclusion test is too coarse to avoid discretisation");

	ReachabilityRestriction rr2 = _get_reference_restriction();

	HybridDenotableSet set_to_restrict2(grid);
	HybridBoxes set_to_restrict2_boxToAdjoin;
	set_to_restrict2_boxToAdjoin.insert(make_pair(q1,Box(2,-0.4,0.5,-1.0,4.0)));
	set_to_restrict2_boxToAdjoin.insert(make_pair(q2,Box(3,0.9,1.0,-0.4,-0.2,1.2,1.3)));
	set_to_restrict2.adjoin_outer_approximation(set_to_restrict2_boxToAdjoin,accuracy);
	HybridDenotableSet restricted_set2 = set_to_restrict2;
	rr2.apply_to(restricted_set2);
	ARIADNE_TEST_ASSERT(set_to_restrict2 == restricted_set2);
	ARIADNE_TEST_ASSERT(rr2.has_discretised(q1));
	ARIADNE_TEST_ASSERT(!rr2.has_discretised(q2));

	ARIADNE_PRINT_TEST_CASE_TITLE("Apply to a HybridDenotableSet that WOULD be restricted thanks to a finer restriction accuracy");

	rr2.refine_at(2);
	HybridDenotableSet restricted_set3 = set_to_restrict2;
	rr2.apply_to(restricted_set3);
	ARIADNE_TEST_ASSERT(set_to_restrict2 != restricted_set3);

	ARIADNE_PRINT_TEST_CASE_TITLE("Apply to a list of enclosures that would not be restricted thanks to inclusion test");

	ReachabilityRestriction rr3 = _get_reference_restriction();

	std::list<LocalisedEnclosureType> enclosures1;
	enclosures1.push_back(LocalisedEnclosureType(q1,EnclosureType(Box(2,-0.4,0.5,-1.0,0.0))));
	enclosures1.push_back(LocalisedEnclosureType(q2,EnclosureType(Box(3,0.9,1.0,-0.4,-0.2,1.2,1.3))));

	std::list<LocalisedEnclosureType> restricted_enclosures1 = rr3.filter(enclosures1);
	ARIADNE_TEST_ASSERT(restricted_enclosures1.size() == enclosures1.size());
	ARIADNE_TEST_ASSERT(!rr3.has_discretised(q1));
	ARIADNE_TEST_ASSERT(!rr3.has_discretised(q2));

	ARIADNE_PRINT_TEST_CASE_TITLE("Apply to a list of enclosures that would be restricted");

	std::list<LocalisedEnclosureType> enclosures2;
	enclosures1.push_back(LocalisedEnclosureType(q1,EnclosureType(Box(2,-0.4,0.5,-1.0,0.0))));
	enclosures2.push_back(LocalisedEnclosureType(q2,EnclosureType(Box(3,-5.0,-4.0,-10.0,-9.0,8.0,8.1))));

	std::list<LocalisedEnclosureType> restricted_enclosures2 = rr3.filter(enclosures2);
	ARIADNE_TEST_ASSERT(restricted_enclosures2.size() != enclosures2.size());
	ARIADNE_TEST_ASSERT(!rr3.has_discretised(q1));
	ARIADNE_TEST_ASSERT(rr3.has_discretised(q2));

}


void
TestReachabilityRestriction::test_restricts() {

	ReachabilityRestriction rr1 = _get_reference_restriction();

	DiscreteLocation q1(1);
	DiscreteLocation q2(2);

	HybridGrid grid = rr1.grid();
	int accuracy = rr1.accuracy();

	HybridDenotableSet set1(grid);
	HybridDenotableSet set2(grid);

	ARIADNE_PRINT_TEST_CASE_TITLE("Check a HybridDenotableSet that clearly would not be restricted");

	HybridDenotableSet set_to_restrict1(grid);
	HybridBoxes set_to_restrict1_boxToAdjoin;
	set_to_restrict1_boxToAdjoin.insert(make_pair(q1,Box(2,-0.4,0.5,-1.0,0.0)));
	set_to_restrict1_boxToAdjoin.insert(make_pair(q2,Box(3,0.9,1.0,-0.4,-0.2,1.2,1.3)));
	set_to_restrict1.adjoin_outer_approximation(set_to_restrict1_boxToAdjoin,accuracy);
	ARIADNE_TEST_ASSERT(!rr1.restricts(set_to_restrict1));
	ARIADNE_TEST_ASSERT(!rr1.has_discretised(q1));
	ARIADNE_TEST_ASSERT(!rr1.has_discretised(q2));

	ARIADNE_PRINT_TEST_CASE_TITLE("Check a HybridDenotableSet that would be restricted on q2");

	HybridDenotableSet set_to_restrict2(grid);
	HybridBoxes set_to_restrict2_boxToAdjoin;
	set_to_restrict2_boxToAdjoin.insert(make_pair(q1,Box(2,-0.4,0.5,-1.0,0.0)));
	set_to_restrict2_boxToAdjoin.insert(make_pair(q2,Box(3,0.9,1.0,-0.4,-0.2,1.2,3.3)));
	set_to_restrict2.adjoin_outer_approximation(set_to_restrict2_boxToAdjoin,accuracy);
	ARIADNE_TEST_ASSERT(rr1.restricts(set_to_restrict2));
	ARIADNE_TEST_ASSERT(!rr1.has_discretised(q1));
	ARIADNE_TEST_ASSERT(rr1.has_discretised(q2));
}


void
TestReachabilityRestriction::test_update() {

	ReachabilityRestriction rr1 = _get_reference_restriction();

	DiscreteLocation q1(1);
	DiscreteLocation q2(2);
	DiscreteLocation q3(3);

	HybridGrid grid = rr1.grid();
	int accuracy = rr1.accuracy();

	ARIADNE_PRINT_TEST_CASE_TITLE("Update with a matching grid");

	HybridDenotableSet set1(grid);
	HybridBoxes full_boxes;
	full_boxes.insert(make_pair(q1,Box(2,0.4,0.6,-1.0,2.0)));
	full_boxes.insert(make_pair(q2,Box(3,0.8,1.0,-0.4,-0.3,1.1,1.2)));
	set1.adjoin_outer_approximation(full_boxes,accuracy);
	HybridBoxes preupdate_bb = rr1.bounding_box();
	rr1.update_with(set1);
	ARIADNE_TEST_ASSERT(preupdate_bb.superset(rr1.bounding_box()));
	ARIADNE_TEST_ASSERT(rr1.has_discretised(q1));
	ARIADNE_TEST_ASSERT(rr1.has_discretised(q2));

	ARIADNE_PRINT_TEST_CASE_TITLE("Update with a mismatching grid");

	ReachabilityRestriction rr2 = _get_reference_restriction();

	HybridGrid q1_grid;
	q1_grid[q1] = Grid(1);
	HybridDenotableSet set2(q1_grid);
	HybridBoxes q3_boxes;
	q3_boxes.insert(make_pair(q3,Box(1,0.5,0.6)));
	set2.adjoin_outer_approximation(q3_boxes,accuracy);
	ARIADNE_TEST_FAIL(rr2.update_with(set2));
}


void
TestReachabilityRestriction::test_copy() {

	ReachabilityRestriction rr1 = _get_reference_restriction();

	DiscreteLocation q1(1);
	DiscreteLocation q2(2);

	HybridGrid grid = rr1.grid();

    ReachabilityRestriction rr2 = rr1;
    ARIADNE_TEST_ASSERT(rr2.accuracy() == rr1.accuracy());
    ARIADNE_TEST_ASSERT(rr2.grid() == rr1.grid());

	ARIADNE_PRINT_TEST_CASE_TITLE("Invariance of the origin after a refinement of the copy");

    rr2.refine_at(2);
    ARIADNE_TEST_EQUAL(rr1.accuracy(),1);

    ARIADNE_PRINT_TEST_CASE_TITLE("Invariance of the origin after a discretisation on the copy");

	HybridDenotableSet set1(grid);
	HybridBoxes boxes;
	boxes.insert(make_pair(q1,Box(2,0.0,0.1,-1.0,1.0)));
	boxes.insert(make_pair(q2,Box(3,-1.0,1.1,0.0,0.1,-1.0,1.0)));

	set1.adjoin_outer_approximation(boxes,rr1.accuracy());
	rr2.update_with(set1);
	ARIADNE_TEST_ASSERT(!rr1.has_discretised(q1));
	ARIADNE_TEST_ASSERT(!rr1.has_discretised(q2));
}


void
TestReachabilityRestriction::test_forward_jump_set() {

	HybridAutomaton sys = _get_system();

	DiscreteLocation q1(1);
	DiscreteLocation q2(2);

	HybridBoxes domain(sys.state_space());
	domain[q1] = Box(2,0.0,10.0,0.0,10.0);
	domain[q2] = Box(2,0.0,10.0,0.0,10.0);

	HybridGrid grid;
	grid[q1] = Grid(Vector<Float>(2,1.0,1.0));
	grid[q2] = Grid(Vector<Float>(2,1.0,1.0));

	int accuracy = 3;

	ReachabilityRestriction rr1(domain,grid,accuracy);

	HybridDenotableSet update_set(grid);
	Box update_set_bx1(2,0.0,6.0,0.0,6.0);
	Box update_set_bx2 = domain[q2];
	update_set[q1].adjoin_over_approximation(update_set_bx1,accuracy);
	update_set[q2].adjoin_over_approximation(update_set_bx2,accuracy);
	rr1.update_with(update_set);

    ARIADNE_PRINT_TEST_CASE_TITLE("The set has no forward jump set due to no activations.");

	HybridDenotableSet set1(grid);
	Box bx1(2,0.0,1.0,0.0,1.0);
	set1[q1].adjoin_over_approximation(bx1,accuracy);

	HybridDenotableSet forward_jump_set1 = rr1.forward_jump_set(set1,sys);

	ARIADNE_TEST_ASSERT(forward_jump_set1.empty());

    ARIADNE_PRINT_TEST_CASE_TITLE("The set has no forward jump set due to invariants associated with urgent transitions.");

	HybridDenotableSet set2(grid);
	Box bx2(2,5.0,8.0,6.0,7.0);
	set2[q1].adjoin_over_approximation(bx2,accuracy);
	HybridDenotableSet forward_jump_set2 = rr1.forward_jump_set(set2,sys);

	ARIADNE_TEST_ASSERT(forward_jump_set2.empty());

	ARIADNE_PRINT_TEST_CASE_TITLE("The set has a forward jump set");

	HybridDenotableSet set4(grid);
	Box bx4(2,3.0,5.0,2.0,4.0);
	set4[q1].adjoin_over_approximation(bx4,accuracy);
	HybridDenotableSet forward_jump_set4 = rr1.forward_jump_set(set4,sys);

	ARIADNE_TEST_ASSERT(!forward_jump_set4.empty());
}


void
TestReachabilityRestriction::test_backward_jump_set() {

	HybridAutomaton sys = _get_system();

	DiscreteLocation q1(1);
	DiscreteLocation q2(2);

	HybridBoxes domain(sys.state_space());
	domain[q1] = Box(2,0.0,10.0,0.0,10.0);
	domain[q2] = Box(2,0.0,10.0,0.0,10.0);

	HybridGrid grid;
	grid[q1] = Grid(Vector<Float>(2,1.0,1.0));
	grid[q2] = Grid(Vector<Float>(2,1.0,1.0));

	int accuracy = 3;

	ReachabilityRestriction rr1(domain,grid,accuracy);

	HybridDenotableSet update_set(grid);
	Box update_set_bx1(2,1.0,6.0,1.0,6.0);
	Box update_set_bx2 = domain[q2];
	update_set[q1].adjoin_over_approximation(update_set_bx1,accuracy);
	update_set[q2].adjoin_over_approximation(update_set_bx2,accuracy);
	rr1.update_with(update_set);

    ARIADNE_PRINT_TEST_CASE_TITLE("The set has no backward jump set.");

	HybridDenotableSet set1(grid);
	Box bx1(2,5.0,9.0,6.0,10.0);
	set1[q2].adjoin_over_approximation(bx1,accuracy);
	HybridDenotableSet backward_jump_set1 = rr1.backward_jump_set(set1,sys);

	ARIADNE_TEST_ASSERT(backward_jump_set1.empty());

	ARIADNE_PRINT_TEST_CASE_TITLE("The set has a backward jump set");

	HybridDenotableSet set2(grid);
	Box bx2(2,5.0,7.0,2.0,5.0);
	set2[q2].adjoin_over_approximation(bx2,accuracy);
	HybridDenotableSet backward_jump_set2 = rr1.backward_jump_set(set2,sys);

	ARIADNE_TEST_ASSERT(!backward_jump_set2.empty());
}


void
TestReachabilityRestriction::test_projection() {

	ReachabilityRestriction rr1 = _get_reference_restriction();

	DiscreteLocation q1(1);
	DiscreteLocation q2(2);

	RealVariable x("x");
	RealVariable y("y");
	RealVariable z("z");
	List<RealVariable> varlist_q1, varlist_q2;
	varlist_q1.append(x);
	varlist_q1.append(y);
	varlist_q2.append(x);
	varlist_q2.append(y);
	varlist_q2.append(z);

	HybridSpace space;
	space[q1] = 2;

	ARIADNE_PRINT_TEST_CASE_TITLE("The set has no outer intersection with the constraint set.");

	RealExpression expr1 = x;
	List<RealExpression> consexpr1;
	consexpr1.append(expr1);
	VectorFunction cons_f1(consexpr1,varlist_q1);
	Box codomain1(1,2.5,4.0);
	HybridConstraintSet constraint1(space,ConstraintSet(cons_f1,codomain1));

	HybridDenotableSet outer_intersection1 = rr1.outer_intersection_with(constraint1);
	HybridDenotableSet inner_intersection1 = rr1.inner_intersection_with(constraint1);
	HybridDenotableSet outer_difference1 = rr1.outer_difference_from(constraint1);
	HybridDenotableSet inner_difference1 = rr1.inner_difference_from(constraint1);

	ARIADNE_TEST_ASSERT(outer_intersection1.empty());
	ARIADNE_TEST_ASSERT(inner_intersection1.empty());
	ARIADNE_TEST_ASSERT(!outer_difference1.empty());
	ARIADNE_TEST_ASSERT(!inner_difference1.empty());

    ARIADNE_PRINT_TEST_CASE_TITLE("The denotable set has only an outer intersection with the constraint set.");

	RealExpression expr2 = x;
	List<RealExpression> consexpr2;
	consexpr2.append(expr2);
	VectorFunction cons_f2(consexpr2,varlist_q1);
	Box codomain2(1,0.5,1.0);

	HybridConstraintSet constraint2(space,ConstraintSet(cons_f2,codomain2));

	HybridDenotableSet outer_intersection2 = rr1.outer_intersection_with(constraint2);
	HybridDenotableSet inner_intersection2 = rr1.inner_intersection_with(constraint2);
	HybridDenotableSet outer_difference2 = rr1.outer_difference_from(constraint2);
	HybridDenotableSet inner_difference2 = rr1.inner_difference_from(constraint2);

	ARIADNE_TEST_ASSERT(!outer_intersection2.empty());
	ARIADNE_TEST_ASSERT(inner_intersection2.empty());
	ARIADNE_TEST_ASSERT(!outer_difference2.empty());
	ARIADNE_TEST_ASSERT(!inner_difference2.empty());

    ARIADNE_PRINT_TEST_CASE_TITLE("The denotable set has inner intersection with the constraint set.");

	ReachabilityRestriction rr3 = _get_reference_restriction();
    rr3.refine_at(4);

	HybridDenotableSet outer_intersection3 = rr3.outer_intersection_with(constraint2);
	HybridDenotableSet inner_intersection3 = rr3.inner_intersection_with(constraint2);
	HybridDenotableSet outer_difference3 = rr3.outer_difference_from(constraint2);
	HybridDenotableSet inner_difference3 = rr3.inner_difference_from(constraint2);

	ARIADNE_TEST_ASSERT(!outer_intersection3.empty());
	ARIADNE_TEST_ASSERT(!inner_intersection3.empty());
	ARIADNE_TEST_ASSERT(!outer_difference3.empty());
	ARIADNE_TEST_ASSERT(!outer_difference3.empty());

    ARIADNE_TEST_ASSERT(!inner_intersection3.bounding_box().superset(outer_intersection3.bounding_box()));

    ARIADNE_PRINT_TEST_CASE_TITLE("The denotable set has no inner difference from the constraint set.");

	ReachabilityRestriction rr4 = _get_reference_restriction();
    rr4.refine_at(4);

    HybridSpace space4 = rr4.space();

	RealExpression expr4_x = x;
	RealExpression expr4_y = y;
	RealExpression expr4_z = z;
	List<RealExpression> consexpr_q1_4, consexpr_q2_4;
	consexpr_q1_4.append(expr4_x);
	consexpr_q1_4.append(expr4_y);
	consexpr_q2_4.append(expr4_x);
	consexpr_q2_4.append(expr4_y);
	consexpr_q2_4.append(expr4_z);
	VectorFunction cons_q1_f4(consexpr_q1_4,varlist_q1);
	VectorFunction cons_q2_f4(consexpr_q2_4,varlist_q2);
	Box codomain_q1_4(2,-1.1,1.1,-2.1,3.1);
	Box codomain_q2_4(3,0.4,2.1,-1.1,0.1,0.9,1.6);
	HybridConstraintSet constraint4;
	constraint4[q1] = ConstraintSet(cons_q1_f4,codomain_q1_4);
	constraint4[q2] = ConstraintSet(cons_q2_f4,codomain_q2_4);

	HybridDenotableSet outer_intersection4 = rr4.outer_intersection_with(constraint4);
	HybridDenotableSet inner_intersection4 = rr4.inner_intersection_with(constraint4);
	HybridDenotableSet outer_difference4 = rr4.outer_difference_from(constraint4);
	HybridDenotableSet inner_difference4 = rr4.inner_difference_from(constraint4);

	ARIADNE_TEST_ASSERT(!outer_intersection4.empty());
	ARIADNE_TEST_ASSERT(!inner_intersection4.empty());
	ARIADNE_TEST_ASSERT(!outer_difference4.empty());
	ARIADNE_TEST_ASSERT(inner_difference4.empty());

    ARIADNE_PRINT_TEST_CASE_TITLE("The denotable set has no outer difference from the constraint set.");

	ReachabilityRestriction rr5 = _get_reference_restriction();
    rr5.refine_at(5);

	HybridDenotableSet outer_intersection5 = rr5.outer_intersection_with(constraint4);
	HybridDenotableSet inner_intersection5 = rr5.inner_intersection_with(constraint4);
	HybridDenotableSet outer_difference5 = rr5.outer_difference_from(constraint4);
	HybridDenotableSet inner_difference5 = rr5.inner_difference_from(constraint4);

	ARIADNE_TEST_ASSERT(!outer_intersection5.empty());
	ARIADNE_TEST_ASSERT(!inner_intersection5.empty());
	ARIADNE_TEST_ASSERT(outer_difference5.empty());
	ARIADNE_TEST_ASSERT(inner_difference5.empty());
}

int main() {
    TestReachabilityRestriction().test();
    return ARIADNE_TEST_FAILURES;
}

