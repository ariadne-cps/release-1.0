/***************************************************************************
 *            test_function_set.cc
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

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <fenv.h>

#include "config.h"
#include "taylor_model.h"
#include "formula.h"
#include "function.h"
#include "predicate.h"
#include "real.h"
#include "expression.h"
#include "function_set.h"

#include "test.h"


using namespace std;
using namespace Ariadne;

class TestFunctionSet
{
  public:
    void test();
  private:
    void test_constraint_set_1D();
};

void TestFunctionSet::test()
{
    ARIADNE_TEST_CALL(test_constraint_set_1D());
}

void TestFunctionSet::test_constraint_set_1D()
{
	RealVariable x("x");
	List<RealVariable> varlist;
	varlist.append(x);

	RealExpression id_x = x;
	List<RealExpression> consexpr;
	consexpr.append(id_x);

	VectorFunction cons_f(consexpr,varlist);
	Box codomain(1,0.0,1.0);

	ConstraintSet cons(cons_f,codomain);

	Box check_box1(1,-0.1,1.1);
    ARIADNE_TEST_ASSERT(definitely(cons.overlaps(check_box1)));
    ARIADNE_TEST_ASSERT(!possibly(cons.disjoint(check_box1)));
    ARIADNE_TEST_ASSERT(!possibly(cons.covers(check_box1)));
    Box check_box2(1,-1.0,-0.1);
    ARIADNE_TEST_ASSERT(!possibly(cons.overlaps(check_box2)));
    ARIADNE_TEST_ASSERT(definitely(cons.disjoint(check_box2)));
    ARIADNE_TEST_ASSERT(!possibly(cons.covers(check_box2)));
    Box check_box3(1,0.1,0.9);
    ARIADNE_TEST_ASSERT(definitely(cons.overlaps(check_box3)));
    ARIADNE_TEST_ASSERT(!possibly(cons.disjoint(check_box3)));
    ARIADNE_TEST_ASSERT(definitely(cons.covers(check_box3)));
    Box check_box4(1,1.0,2.0);
    ARIADNE_TEST_ASSERT(indeterminate(cons.overlaps(check_box4)));
    ARIADNE_TEST_ASSERT(indeterminate(cons.disjoint(check_box4)));
    ARIADNE_TEST_ASSERT(!possibly(cons.covers(check_box4)));
    Box check_box5(1,0.0,1.0);
    ARIADNE_TEST_ASSERT(definitely(cons.overlaps(check_box5)));
    ARIADNE_TEST_ASSERT(!possibly(cons.disjoint(check_box5)));
    ARIADNE_TEST_ASSERT(indeterminate(cons.covers(check_box5)));
}


int main() {
    TestFunctionSet().test();

    return ARIADNE_TEST_FAILURES;
}

