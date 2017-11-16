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
    void test_constant_derivative();
};

void
TestIntegrator::test()
{
    ARIADNE_TEST_CALL(test_constant_derivative());
}

void TestIntegrator::test_constant_derivative() {
    /*VectorFunction f={o*2,o*3};
    ARIADNE_TEST_PRINT(f);
    ExactBoxType d={ExactIntervalType(0.0,1.0),ExactIntervalType(-0.5,1.5)};
    FloatDP h=0.25;
    ValidatedVectorFunctionModelDP flow=integrator_ptr->flow_step(f,d,h);
    EffectiveVectorFunction expected_flow={x0+2*t,y0+3*t};
    ARIADNE_TEST_PRINT(flow);
    ARIADNE_TEST_PRINT(expected_flow);
    ARIADNE_TEST_PRINT(flow.errors());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,norm(flow-expected_flow),1e-8);
     */
}


int main() {

    TaylorIntegrator integrator(1);
    TestIntegrator(integrator).test();
    return ARIADNE_TEST_FAILURES;
}
