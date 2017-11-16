/***************************************************************************
 *            test_integrator.cc
 *
 *  Copyright 2017  Luca Geretti
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
#include <iomanip>
#include "multi_index.h"
#include "taylor_model.h"
#include "zonotope.h"
#include "taylor_set.h"
#include "taylor_function.h"
#include "graphics.h"
#include "integrator.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

class TestIntegrator {
  public:
    TestIntegrator(const IntegratorInterface& i)
            : integrator(i)
    {
        o=ScalarFunction::constant(2,1);
        x=ScalarFunction::coordinate(2,0);
        y=ScalarFunction::coordinate(2,1);
        x0=ScalarFunction::coordinate(3,0);
        y0=ScalarFunction::coordinate(3,1);
        t=ScalarFunction::coordinate(3,2);
    }

    void test();
  private:
    const IntegratorInterface& integrator;
    ScalarFunction o,x,y,x0,y0,t;
private:
    void test_bounds();
    void test_constant_derivative();
    void test_flow_at_step();
    void test_continuous_step();
};

void
TestIntegrator::test()
{
    ARIADNE_TEST_CALL(test_bounds());
    ARIADNE_TEST_CALL(test_constant_derivative());
    ARIADNE_TEST_CALL(test_flow_at_step());
    ARIADNE_TEST_CALL(test_continuous_step());
}

void TestIntegrator::test_bounds() {
    VectorFunction f=join(o*2,o*3);
    ARIADNE_TEST_PRINT(f);

    Box initial(2, 0.0,1.0, -0.5,1.5);
    float h=1.0;

    Pair<float,Box> step_and_bounds = integrator.flow_bounds(f, initial, h);
    float actual_step = step_and_bounds.first;
    Box computed_bounds = step_and_bounds.second;

    Box expected_bounds(2, 0.0,3.0, -0.5,4.5);

    ARIADNE_TEST_PRINT(step_and_bounds);

    ARIADNE_TEST_EQUAL(actual_step,h);
    ARIADNE_TEST_ASSERT(computed_bounds.superset(expected_bounds));
}

void TestIntegrator::test_constant_derivative() {
    VectorFunction f=join(o*2,o*3);
    ARIADNE_TEST_PRINT(f);

    Box d(2, 0.0,1.0, -0.5,1.5);
    float h=1.0;
    VectorTaylorFunction flow=integrator.flow(f,d,h);
    VectorTaylorFunction expected_flow(flow.domain(),join(x0+2*t,y0+3*t));

    ARIADNE_TEST_PRINT(flow);
    ARIADNE_TEST_PRINT(expected_flow);
    ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow - expected_flow),1e-10);
}


void TestIntegrator::test_flow_at_step() {
    VectorFunction f=join(o*2,o*3);
    ARIADNE_TEST_PRINT(f);

    Box d(2, 0.0,1.0, -0.5,1.5);
    float h=1.0;
    VectorTaylorFunction flow_at_step = integrator.flow_at_step(f,d,h);
    VectorTaylorFunction expected_flow_at_step(d,join(x+2,y+3));

    ARIADNE_TEST_PRINT(flow_at_step);
    ARIADNE_TEST_PRINT(expected_flow_at_step);
    ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow_at_step - expected_flow_at_step),1e-10);
}

void TestIntegrator::test_continuous_step() {
    VectorFunction f=join(o*2,o*3);
    ARIADNE_TEST_PRINT(f);

    Vector<TaylorModel> models(2);
    models[0] = TaylorModel(Expansion<Real>(3,3, 0,0,0,0.5, 1,0,0,0.375, 0,1,0,0.125),0.0);
    models[1] = TaylorModel(Expansion<Real>(3,4, 0,0,0,0.5, 1,0,0,0.25, 0,1,0,0.625, 0,0,1,0.125),0.0);
    TaylorSet starting_set(models);

    ARIADNE_TEST_PRINT(starting_set);

    float h=1.0;
    VectorTaylorFunction flow = integrator.flow(f,starting_set.bounding_box(),h);

    ARIADNE_TEST_PRINT(flow);

    VectorTaylorFunction flow_at_step = partial_evaluate(flow,f.result_size(),h);
    VectorTaylorFunction flow_at_deltat = partial_evaluate(flow,f.result_size(),Interval(0,h));

    TaylorSet finishing_set = apply(flow_at_step, starting_set);

    Box expected_finishing_set_bounds(2, 2.0,3.0, 2.5,4.5);

    ARIADNE_TEST_PRINT(finishing_set);
    ARIADNE_TEST_PRINT(finishing_set.bounding_box());
    ARIADNE_TEST_ASSERT(finishing_set.bounding_box().superset(expected_finishing_set_bounds));

    TaylorSet reached_set = apply(flow_at_deltat, starting_set);

    Box expected_flow_set_bounds(2, 0.0,3.0, -0.5,4.5);

    ARIADNE_TEST_PRINT(reached_set);
    ARIADNE_TEST_PRINT(reached_set.bounding_box());
    ARIADNE_TEST_ASSERT(reached_set.bounding_box().superset(expected_flow_set_bounds));
}

int main() {

    TaylorIntegrator integrator(1);
    TestIntegrator(integrator).test();
    return ARIADNE_TEST_FAILURES;
}
