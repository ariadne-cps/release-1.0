/***************************************************************************
 *            test_discretised_evolution.cc
 *
 *  Copyright  2006-11  Pieter Collins, Alberto Casagrande
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

#include <fstream>

#include "tuple.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "taylor_set.h"
#include "zonotope.h"
#include "list_set.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "denotable_set.h"
#include "hybrid_automaton.h"
#include "hybrid_evolver.h"
#include "orbit.h"
#include "graphics.h"

#include "models.h"

#include "logging.h"

#include "test.h"

using namespace Ariadne;
using namespace std;
using Models::VanDerPol;

class TestDiscretisedEvolution
{
  public:
    void test() const;
  private:
    void test_hybrid_time() const;
};

int main()
{
    TestDiscretisedEvolution().test();
    return ARIADNE_TEST_FAILURES;
}

void TestDiscretisedEvolution::test() const
{
    ARIADNE_TEST_CALL(test_hybrid_time());
}


void TestDiscretisedEvolution::test_hybrid_time() const
{
    typedef TaylorSet EnclosureType;
    typedef HybridBasicSet<EnclosureType> LocalisedEnclosureType;

    cout << __PRETTY_FUNCTION__ << endl;

    // Set up the evolution parameters and grid
    Float time(1.0);
    uint steps(6);
    Float maximum_step_size(0.125);
    int depth=8;
    DiscreteLocation location(1);
    DiscreteEvent event(1);
    Grid grid(2);

    // Set up the vector field
    Float a=1.5; Float b=0.375;
    Vector<Float> p(2); p[0]=a; p[1]=b;

    VectorUserFunction<Henon> henon(p);
    cout << "henon=" << henon << endl;
    HybridAutomaton ha("Henon");
    ha.new_mode(location,IdentityFunction(2));
    ha.new_transition(event,location,location,henon,VectorConstantFunction(Vector<Float>(1,1.0),2),PERMISSIVE);

    // Set up the evaluators
    HybridEvolver evolver(ha);

    evolver.settings().minimum_discretised_enclosure_widths[location]=Vector<Float>(2,0.5);
    evolver.settings().set_hybrid_maximum_step_size(maximum_step_size);

    // Define a bounding box for the evolution
    std::cout<<"making bounding_box"<<std::endl;
    Box bounding_box=make_box("[-4,4]x[-4,4]") ;
    std::cout<<"bounding_box="<<bounding_box<<"\n"<<std::endl;
    //Box eps_bounding_box=bounding_box.neighbourhood(0.1);

    // Define the initial cell
    Box box=make_box("[1.0001,1.0002]x[0.5001,0.5002]");
    LocalisedBox hybrid_initial_set(location,box);
    cout << "hybrid_initial_set=" << hybrid_initial_set << endl << endl;
    //[1.00098:1.00122],
    HybridTime htime(time,steps);
    cout << "hybrid_time=" << htime << endl << endl;

    // Compute the reachable sets
    cout << "Computing evolution... " << flush;
    // evolver.verbosity=1;
    LocalisedEnclosureType localised_initial_enclosure(location,EnclosureType(box));
    Orbit<LocalisedEnclosureType> evolve_orbit = evolver.orbit(localised_initial_enclosure,htime,UPPER_SEMANTICS);
    cout << "done." << endl;

    cout << "enclosure_orbit="<<evolve_orbit<<endl;
    cout << evolve_orbit.reach()<<endl;

    EnclosureType const& initial_set=evolve_orbit.initial().second.range();
    ListSet<EnclosureType> const& reach_set=evolve_orbit.reach()[location];
    ListSet<EnclosureType> const& intermediate_set=evolve_orbit.intermediate()[location];
    ListSet<EnclosureType> const& final_set=evolve_orbit.final()[location];
    HybridGrid hagrid(ha.state_space(),Float(1.0));

    // Compute the reachable sets
    cout << "Computing discretised evolution... " << flush;


    HybridDenotableSet reach(outer_approximation(evolve_orbit.reach(),hagrid,depth));
    HybridDenotableSet final(outer_approximation(evolve_orbit.final(),hagrid,depth));
    cout << "done." << endl;

    DenotableSetType const& reach_cells=reach[location];
    DenotableSetType const& final_cells=final[location];

    cout << "initial_set=" << initial_set.range() << endl << endl;
    cout << "reach_set=" << reach_set << endl << endl;
    cout << "reach_cells=" << reach_cells << endl << endl;
    cout << "final_set=" << final_set << endl << endl;
    cout << "final_cells=" << final_cells << endl << endl;

    // Print the initial, evolve and reach sets
    {
        Figure fig;
        fig.set_bounding_box(Box(2,Interval(-3,3)));
        fig << line_style(true);
        fig << line_style(true);
        fig << fill_colour(cyan) << reach_cells;
        fig << fill_colour(green) << final_cells;
        fig.write("test_discrete_evolver-hybrid-cells");
    }

    // Print the initial, evolve and reach sets
    {
        Figure fig;
        fig.set_bounding_box(Box(2,Interval(-3,3)));
        fig << line_style(true);
        fig << fill_colour(cyan) << reach_set;
        fig << fill_colour(magenta) << intermediate_set;
        fig << fill_colour(yellow) << initial_set;
        fig << fill_colour(green) << final_set;
        fig.write("test_discrete_evolver-hybrid-sets");
    }
}

