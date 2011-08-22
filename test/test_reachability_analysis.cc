/***************************************************************************
 *            test_reachability_analysis.cc
 *
 *  Copyright  2006-8  Pieter Collins
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
#include <string>

#include "taylor_set.h"
#include "grid_set.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"
#include "settings.h"
#include "hybrid_evolver.h"
#include "reachability_analyser.h"
#include "graphics.h"
#include "logging.h"

#include "models.h"

#include "test.h"

using namespace std;
using namespace Ariadne;
using Ariadne::Models::Henon;



class TestReachabilityAnalysis 
{  
    Grid grid;
    Interval bound;
    HybridImageSet initial_set;
    HybridTime reach_time;
 
    typedef TaylorSet EnclosureType;
    typedef HybridBasicSet<TaylorSet> HybridEnclosureType;

  public:
    static HybridReachabilityAnalyser build_analyser()
    {
        HybridReachabilityAnalyser analyser(build_system());
        cout << "Done building analyser\n";
        return analyser;
    }

    static HybridAutomaton build_system() {

    	HybridAutomaton system;

    	std::cout<<std::setprecision(20);
		std::cerr<<std::setprecision(20);
		std::clog<<std::setprecision(20);
		DiscreteLocation location(1);

		Matrix<Float> A=Matrix<Float>("[-0.5,-1.0;1.0,-0.5]");
		Vector<Float> b=Vector<Float>("[0.0,0.0]");
		VectorAffineFunction aff(A,b);
		system.new_mode(location,aff);
		cout << "Done building system\n";

		cout << "system=" << system << endl;

		return system;
    }

    TestReachabilityAnalysis()
        : grid(2),
          bound(-4,4),
          reach_time(4.0,3)
    {
        cout << "Done initialising variables\n";

        ImageSet initial_box(make_box("[2.0,2.0]x[0.0,0.0]"));
        DiscreteLocation location(1);
        initial_set[location]=initial_box;
        cout << "Done creating initial set\n" << endl;

        cout << "initial_set=" << initial_set << endl;
    }

    template<class S> void plot(const char* name, const Box& bounding_box, const S& set) {
        Figure g;
        g << fill_colour(white) << bounding_box << line_style(true);
        g << fill_colour(blue) << set;
        g.write(name);
    }

    template<class S, class IS> void plot(const char* name, const Box& bounding_box, const S& set, const IS& initial_set) {
        Figure g;
        g << fill_colour(white) << bounding_box;
        g << line_style(true);
        g << fill_colour(red) << set;
        g << fill_colour(blue);
        g << initial_set;
        g.write(name);
    }

    template<class ES, class RS, class IS> void plot(const char* name, const Box& bounding_box, 
                                                     const ES& evolve_set, const RS& reach_set, const IS& initial_set) {
        Figure g;
        g << fill_colour(white) << bounding_box;
        g << line_style(true);
        g << fill_colour(green) << reach_set;
        g << fill_colour(red) << evolve_set;
        g << fill_colour(blue) << initial_set;
        g.write(name);
    }

    void test_lower_reach_evolve() {  

    	HybridReachabilityAnalyser analyser(build_analyser());

        DiscreteLocation loc(1);
        Box bounding_box(2,bound);
        analyser.verbosity=0;
        analyser.settings().domain_bounds[loc] = bounding_box;
        cout << "Computing timed evolve set" << endl;
        HybridGridTreeSet hybrid_lower_evolve=analyser.lower_evolve(initial_set,reach_time);
        cout << "Computing timed reachable set" << endl;
        HybridGridTreeSet hybrid_lower_reach=analyser.lower_reach(initial_set,reach_time);
        GridTreeSet& lower_evolve=hybrid_lower_evolve[loc];
        GridTreeSet& lower_reach=hybrid_lower_reach[loc];
        cout << "Evolved to " << lower_evolve.size() << " cells " << endl << endl;
        cout << "Reached " << lower_reach.size() << " cells " << endl << endl;
        plot("test_reachability_analyser-map_lower_reach_evolve.png",bounding_box,lower_evolve,lower_reach,initial_set);
    }
  
    void test_upper_reach_evolve() {  

    	HybridReachabilityAnalyser analyser(build_analyser());

        cout << "Computing timed reachable set" << endl;
        DiscreteLocation loc(1);
        Box bounding_box(2,bound);
        analyser.settings().domain_bounds[loc] = Box(2,-10.0,10.0,-10.0,10.0);
        HybridGridTreeSet upper_evolve_set=analyser.upper_evolve(initial_set,reach_time);
        cout << "upper_evolve_set="<<upper_evolve_set<<std::endl;
        HybridGridTreeSet upper_reach_set=analyser.upper_reach(initial_set,reach_time);
        cout << "upper_reach_set="<<upper_reach_set<<std::endl;
 
        const GridTreeSet& upper_evolve=upper_evolve_set[loc];
        const GridTreeSet& upper_reach=upper_reach_set[loc];
        ImageSet& initial=initial_set[loc];
        //cout << "Reached " << upper_reach.size() << " cells out of " << upper_reach.capacity() << endl << endl;
        plot("test_reachability_analyser-map_upper_reach_evolve.png",bounding_box,upper_evolve,upper_reach,initial);
    }
  
    void test_chain_reach() {  

    	HybridReachabilityAnalyser analyser(build_analyser());

        cout << "Computing chain reachable set" << endl;
        DiscreteLocation loc(1);
        HybridBoxes bounding_boxes
            =Ariadne::bounding_boxes(build_system().state_space(),bound);
        Box bounding_box=bounding_boxes[loc];

        analyser.verbosity=0;
        analyser.settings().transient_time=4.0;
        analyser.settings().lock_to_grid_time=1.0;
        analyser.settings().domain_bounds = bounding_boxes;
        cout << analyser.settings();
        HybridGridTreeSet chain_reach_set=analyser.outer_chain_reach(initial_set);
        plot("test_reachability_analyser-map_chain_reach.png",bounding_box,chain_reach_set[loc],initial_set[loc]);
    }
  
    void test() {
        ARIADNE_TEST_CALL(test_lower_reach_evolve());
        ARIADNE_TEST_CALL(test_upper_reach_evolve());
        ARIADNE_TEST_CALL(test_chain_reach());
    }

};


int main(int nargs, const char* args[]) 
{
    TestReachabilityAnalysis().test();
    if(ARIADNE_TEST_SKIPPED) { cerr << "INCOMPLETE "; }
    return ARIADNE_TEST_FAILURES;
}

