/***************************************************************************
 *            test_taylor_set.cc
 *
 *  Copyright 2008  Pieter Collins
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

#include "test.h"
using namespace std;
using namespace Ariadne;

class TestTaylorSet {
  public:
    void test();
  private:
    void test_linearise();
    void test_split();
    void test_subsume();
};

void
TestTaylorSet::test()
{
    ARIADNE_TEST_CALL(test_linearise());
    ARIADNE_TEST_CALL(test_split());
    ARIADNE_TEST_CALL(test_subsume());
}


void
TestTaylorSet::test_subsume()
{
    TaylorSet ts1=TaylorSet(2,2,2, 0.0,1.0,0.5,0.0,0.0,0.0, 0.25, 0.0,0.5,1.0,1.0,0.0,0.0, 0.375);
    TaylorSet cts1=TaylorSet(2,4,2, 0.0, 1.0,0.5,0.25,0.0,  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,  0.0,
                                    0.0, 0.5,1.0,0.0,0.375, 1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,  0.0);
    ARIADNE_TEST_EQUAL(ts1.subsume(),cts1);

    TaylorSet ts2=TaylorSet(2,2,2, 0.0,1.0,0.5,0.0,0.0,0.0, 0.0, 0.0,0.5,1.0,1.0,0.0,0.0, 0.375);
    TaylorSet cts2=TaylorSet(2,3,2, 0.0, 1.0,0.5,0.0,  0.0,0.0,0.0,0.0,0.0,0.0,  0.0,
                                    0.0, 0.5,1.0,0.375, 1.0,0.0,0.0,0.0,0.0,0.0,  0.0);
    ARIADNE_TEST_EQUAL(ts2.subsume(),cts2);
}


void
TestTaylorSet::test_linearise()
{
    TaylorSet ts(2,2,2, 0.0,1.0,0.25,0.0,0.0,0.0, 0.0, 0.0,0.5,1.0,1.0,0.0,0.0, 0.0);

    Zonotope z=zonotope(ts);
    Box b=ts.bounding_box();

    ARIADNE_TEST_EQUAL(ts.linearise(),TaylorSet(2,2,1, 0.0,1.0,0.25, 0.0, 0.0,0.5,1.0, 1.0));

    Box bounding_box=ts.bounding_box()+Vector<Interval>(2,Interval(-1,1));
    plot("test_taylor_set-linearise",PlanarProjectionMap(2,0,1),bounding_box,
         Colour(1,0,0),b,Colour(1,0,1),z,Colour(0,0,1),ts);
}

void plot(const char* filename, const TaylorSet& set) {
    Figure fig;
    fig.set_bounding_box(set.bounding_box());
    fig.set_line_width(0.0);
    draw(fig,set);
    fig.write(filename);
}


void
TestTaylorSet::test_split()
{
    TaylorSet ts(2,2,2, 0.0,1.0,0.25,0.0,0.0,0.0, 0.0, 0.0,0.5,1.0,1.0,0.0,0.0, 0.0);

    TaylorSet ts1,ts2,ts3,ts4,ts5,ts6;
    make_lpair(ts1,ts2)=ts.split();
    make_lpair(ts3,ts4)=ts2.split();
    make_lpair(ts5,ts6)=ts4.split();

    ARIADNE_TEST_EQUAL(ts.split().first,
        TaylorSet(2,2,2, -0.5,0.5,0.25,0.0,0.0,0.0, 0.0, +0.0,-0.25,1.0,0.25,0.0,0.0, 0.0));
    ARIADNE_TEST_EQUAL(ts.split().second,
        TaylorSet(2,2,2, +0.5,0.5,0.25,0.0,0.0,0.0, 0.0, +0.5,+0.75,1.0,0.25,0.0,0.0, 0.0));

    Box bounding_box=ts.bounding_box()+Vector<Interval>(2,Interval(-1,1));
    plot("test_taylor_set-split",bounding_box,
         Colour(0,0.0,1),ts1,
         Colour(0,0.4,1),ts3,
         Colour(0,0.8,1),ts5,
         Colour(0,0.9,1),ts6);

    // Test split with an error term
    ts=TaylorSet(2,2,2, 0.5,1.0,0.25,0.0,0.0,0.0, 2.5, 0.0,0.5,1.0,1.0,0.0,0.0, 1.0);
    ARIADNE_TEST_EQUAL(ts.split().first,
        TaylorSet(2,2,2, -0.75,1.0,0.25,0.0,0.0,0.0, 1.25, 0.0,0.5,1.0,1.0,0.0,0.0, 1.0));
    ARIADNE_TEST_EQUAL(ts.split().second,
        TaylorSet(2,2,2, 1.75,1.0,0.25,0.0,0.0,0.0, 1.25, 0.0,0.5,1.0,1.0,0.0,0.0, 1.0));
    

}


int main() {
    TestTaylorSet().test();
    return ARIADNE_TEST_FAILURES;
}
