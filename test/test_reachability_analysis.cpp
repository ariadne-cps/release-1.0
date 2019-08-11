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
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"
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
    HybridBoundedConstraintSet initial_set;
    HybridTime reach_time;
    HybridAutomaton system;
 
    typedef TaylorSet EnclosureType;
    typedef HybridBasicSet<TaylorSet> LocalisedEnclosureType;

  public:

    HybridReachabilityAnalyser build_analyser()
    {
    	HybridBoxes domain;
    	DiscreteLocation loc(1);
    	domain.insert(make_pair(loc,Box(2,-1.0,2.1,-0.5,1.1)));
    	int accuracy = 6;

    	cout << "Started building analyser\n";
        HybridReachabilityAnalyser analyser(build_system(),domain,accuracy);
        cout << "Done building analyser\n";
        return analyser;
    }

    static HybridAutomaton build_system() {

    	HybridAutomaton sys;

    	std::cout<<std::setprecision(20);
		std::cerr<<std::setprecision(20);
		std::clog<<std::setprecision(20);
		DiscreteLocation loc(1);

		Matrix<Float> A=Matrix<Float>("[-0.5,-1.0;1.0,-0.5]");
		Vector<Float> b=Vector<Float>("[0.0,0.0]");
		VectorAffineFunction aff(A,b);
		sys.new_mode(loc,aff);
		cout << "Done building system\n";

		cout << "system=" << sys << endl;

		return sys;
    }

    TestReachabilityAnalysis()
        : reach_time(4.0,3)
    {
        cout << "Done initialising variables\n";

        DiscreteLocation loc(1);
        HybridSpace hspace;
        hspace.insert(make_pair(loc,2));
        initial_set = HybridBoundedConstraintSet(hspace);
        initial_set[loc]=Box(2,2.0,2.0,0.0,0.0);
        cout << "Done creating initial set\n" << endl;

        cout << "initial_set=" << initial_set << endl;
    }

    template<class S> void plot(const char* name, const Box& bounding_box, const S& set) {
        Figure g;
        g.set_bounding_box(bounding_box);
        g << fill_colour(white) << bounding_box << line_style(true);
        g << fill_colour(blue) << set;
        g.write(name);
    }

    template<class ES, class RS> void plot(const char* name, const Box& bounding_box,
                                                     const ES& evolve_set, const RS& reach_set) {
        Figure g;
        g.set_bounding_box(bounding_box);
        g << fill_colour(white) << bounding_box;
        g << line_style(true);
        g << fill_colour(green) << reach_set;
        g << fill_colour(red) << evolve_set;
        g.write(name);
    }

    void test_lower_reach_evolve() {  

    	HybridReachabilityAnalyser analyser(build_analyser());

        DiscreteLocation loc(1);
        Box graphics_box(2,-1.0,2.1,-0.5,1.1);
        cout << "Computing timed evolve set" << endl;
        cout << "System: " << build_system() << endl;
        cout << "Initial set: " << initial_set << endl;
        HybridDenotableSet hybrid_lower_evolve=analyser.lower_evolve(initial_set,reach_time);
        cout << "Computing timed reachable set" << endl;
        HybridDenotableSet hybrid_lower_reach=analyser.lower_reach(initial_set,reach_time);
        DenotableSetType& lower_evolve=hybrid_lower_evolve[loc];
        DenotableSetType& lower_reach=hybrid_lower_reach[loc];
        cout << "Evolved to " << lower_evolve.size() << " cells " << endl;
        cout << "Reached " << lower_reach.size() << " cells " << endl << endl;
        //ARIADNE_TEST_EQUAL(lower_evolve.size(),1);
        //ARIADNE_TEST_EQUAL(lower_reach.size(),173);
        plot("test_reachability_analyser-map_lower_reach_evolve.png",graphics_box,lower_evolve,lower_reach);
    }
  
    void test_upper_reach_evolve() {  

    	HybridReachabilityAnalyser analyser(build_analyser());

        cout << "Computing timed reachable set" << endl;
        DiscreteLocation loc(1);
        Box graphics_box(2,-1.0,2.1,-0.5,1.1);
        HybridDenotableSet upper_evolve_set=analyser.upper_evolve(initial_set,reach_time);
        // cout << "upper_evolve_set="<<upper_evolve_set<<std::endl;
        HybridDenotableSet upper_reach_set=analyser.upper_reach(initial_set,reach_time);
        // cout << "upper_reach_set="<<upper_reach_set<<std::endl;
 
        const DenotableSetType& upper_evolve=upper_evolve_set[loc];
        const DenotableSetType& upper_reach=upper_reach_set[loc];;
        cout << "Evolved to " << upper_evolve.size() << " cells"  << endl;
        cout << "Reached " << upper_reach.size() << " cells" << endl << endl;
        //ARIADNE_TEST_EQUAL(upper_evolve.size(),6);
        //ARIADNE_TEST_EQUAL(upper_reach.size(),240);
        plot("test_reachability_analyser-map_upper_reach_evolve.png",graphics_box,upper_evolve,upper_reach);
    }
  
    void test_chain_reach() {  

    	HybridReachabilityAnalyser analyser(build_analyser());

        cout << "Computing chain reachable set" << endl;
        DiscreteLocation loc(1);
        Box graphics_box(2,-1.0,2.1,-0.5,1.1);
        // analyser.verbosity=4;
        HybridDenotableSet chain_reach_set=analyser.outer_chain_reach(initial_set);
        cout << "Reached " << chain_reach_set.size() << " cells" << endl << endl;
        //ARIADNE_TEST_EQUAL(chain_reach_set.size(),277);
        plot("test_reachability_analyser-map_chain_reach.png",graphics_box,chain_reach_set[loc]);
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

