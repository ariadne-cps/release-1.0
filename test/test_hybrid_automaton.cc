/***************************************************************************
 *            test_hybrid_automaton.cc
 *
 *  Copyright  2006-7  Alberto Casagrande,  Pieter Collins
 *  Email  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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

#include "test_float.h"
#include "geometry/set_interface.h"
#include "geometry/set_reference.h"
#include "geometry/hybrid_set.h"
#include "geometry/polyhedron.h"
#include "geometry/polyhedral_set.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/hybrid_automaton.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

template<class R> int test_hybrid_automaton();
  
int main() {
  return test_hybrid_automaton<Flt>();
}

template<class R>
int test_hybrid_automaton() 
{
  
  Box<R> r("[-1,1]x[-1,1]");
  cout << "r=" << r << endl;

  AffineFunction<R> dynamic(Matrix<R>("[0.25,-1.00;1.00,0.25]"),Vector<R>("[0.00,0.00]"));
  cout << "dynamic=" << dynamic << endl;
  AffineFunction<R> reset(Matrix<R>("[-0.125,0;0,-0.125]"),Vector<R>("[0,0]"));
  cout << "reset=" << reset << endl;
  
  AffineFunction<R> invariant1(Matrix<R>("[-1,0]"),Vector<R>("[1]"));
  AffineFunction<R> invariant2(Matrix<R>("[-1,0]"),Vector<R>("[4]"));
  cout << "invariant1=" << invariant1 << endl;
  cout << "invariant2=" << invariant2 << endl;
  AffineFunction<R> activation12(Matrix<R>("[-1,0]"),Vector<R>("[1.125]"));
  AffineFunction<R> guard21(Matrix<R>("[1,0]"),Vector<R>("[-1]"));
  cout << "activation12=" << activation12 << endl;
  cout << "activation21=" << guard21 << endl;
  cout << endl;
  
  HybridAutomaton<R> automaton("Constraint-based affine test automaton");
  DiscreteState dstate1(0);
  DiscreteState dstate2(1);
  const DiscreteMode<R>& mode1=automaton.new_mode(dstate1,dynamic,invariant1);
  const DiscreteMode<R>& mode2=automaton.new_mode(dstate2,dynamic,invariant2);
  DiscreteEvent event(5);
  const DiscreteTransition<R>& transition12=automaton.new_transition(event,dstate1,dstate2,reset,activation12);
  const DiscreteTransition<R>& transition21=automaton.new_transition(event,dstate2,dstate1,reset,guard21);
  
  cout << mode1  <<  "\n" << mode2 << "\n" << transition12 << "\n" << transition21 << endl;

  return 0;
}
