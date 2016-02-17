/***************************************************************************
 *            test_vector.cc
 *
 *  Copyright  2006  Pieter Collins, Alberto Casagrande
 *  Email Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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
#include <cassert>


#include "config.h"
#include "numeric.h"
#include "vector.h"

#include "test.h"

using namespace std;
using namespace Ariadne;

class TestVector {
  public:
    void test();
  private:
    void test_concept();
    void test_misc();
};

void 
TestVector::test() 
{
    ARIADNE_TEST_CALL(test_concept());
    ARIADNE_TEST_CALL(test_misc());
}

void
TestVector::test_concept()
{
    Float fx;
    Interval ix;
    Vector<Float> fv, fvv;
    Vector<Interval> iv, ivv;
  
    fv=fv+fv; iv=fv+fv; iv=fv+iv; iv=iv+fv; iv=iv+iv; 
    fv=fv-fv; iv=fv-fv; iv=fv-iv; iv=iv-fv; iv=iv-iv; 

    fv=fx*fv; iv=fx*fv; iv=fx*iv; iv=ix*fv; iv=ix*iv; 
    fv=fv*fx; iv=fv*fx; iv=fv*ix; iv=iv*fx; iv=iv*ix; 
    fv=fv/fx; iv=fv/fx; iv=fv/ix; iv=iv/fx; iv=iv/ix; 
    
    // Test variadic constructor and comma operator
    fv = Vector<Float>(3);
    fv[0] = 1.0; fv[1] = 2.0; fv[2] = 3.3;
    fvv = Vector<Float>(3, 1.0, 2.0, 3.3);
    ARIADNE_TEST_EQUAL(fv,fvv);
    fvv = Vector<Float>(3);
    ARIADNE_TEST_COMPARE(fv,!=,fvv);
    fvv = 1.0, 2.0, 3.3;
    ARIADNE_TEST_EQUAL(fv,fvv);
       
    iv = Vector<Interval>(3);
    iv[0] = 1.0; iv[1] = 2.0; iv[2] = 3.3;
    ivv = Vector<Interval>(3, 1.0,1.0, 2.0,2.0, 3.3,3.3);
    ARIADNE_TEST_EQUAL(iv,ivv);
    ivv = Vector<Interval>(3);
    ARIADNE_TEST_COMPARE(iv,!=,ivv);
    ivv = Interval(1.0,1.0), Interval(2.0,2.0), Interval(3.3,3.3);
    ARIADNE_TEST_EQUAL(iv,ivv);
       
#ifdef HAVE_RATIONAL
    Vector<Rational> rv(3);
    rv[0] = 1.0; rv[1] = 2.0; rv[2] = 3.3;
    Vector<Rational> rvv (3, 1.0, 2.0, 3.3);
    ARIADNE_TEST_EQUAL(rv,rvv);
    rvv = Vector<Rational>(3);
    ARIADNE_TEST_COMPARE(rv,!=,rvv);
    rvv = 1.0, 2.0, 3.3;
    ARIADNE_TEST_EQUAL(rv,rvv);
#endif // HAVE_RATIONAL
     
}


void
TestVector::test_misc()
{
   
    uint n=3;
    Float vptr[3]={-4.0,3.0,1.0};
    Float x=1.5;
    
    // Test constructors for Vector<Float>
    Vector<Float> v0;
    ARIADNE_TEST_EQUAL(v0.size(),0);
    Vector<Float> v1(n,vptr);
    ARIADNE_TEST_EQUAL(v1.size(),n);
    Vector<Float> v2=Vector<Float>("[2.375,4.25,-1.25]");
    ARIADNE_TEST_EQUAL(norm(v1),4);
    ARIADNE_TEST_EQUAL(norm(v2),4.25);

    Vector<Float> v3(1);
    ARIADNE_TEST_EQUAL(v3.size(),1);
    ARIADNE_TEST_EQUAL(v3,Vector<Float>("[0.0]"));
    Vector<Float> v4=v2;
    ARIADNE_TEST_EQUAL(v4,v2);
    
    // Test arithmetic operators for Vector<Float>
    Vector<Float> vf0;
    v1=Vector<Float>("[0.25,-1.5]");
    v2=Vector<Float>("[-0.5,2.25]");
    vf0=-v1;
    ARIADNE_TEST_EQUAL(vf0,Vector<Float>("[-0.25,1.5]"));    
    vf0=v1+v2;
    ARIADNE_TEST_EQUAL(vf0,Vector<Float>("[-0.25,0.75]"));    
    vf0=v1-v2;
    ARIADNE_TEST_EQUAL(vf0,Vector<Float>("[0.75,-3.75]"));    
    vf0=x*v2;
    ARIADNE_TEST_EQUAL(vf0,Vector<Float>("[-0.75,3.375]"));    
    vf0=Vector<Float>(v1)*x;
    ARIADNE_TEST_EQUAL(vf0,Vector<Float>("[0.375,-2.25]"));
    x = 2.0;
    vf0=v1/x;
    ARIADNE_TEST_EQUAL(vf0,Vector<Float>("[0.125,-0.75]"));    
    // Test addition and multiplication with Integers
    int i = 2;
    vf0=v1*i;
    ARIADNE_TEST_EQUAL(vf0,Vector<Float>("[0.5,-3.0]"));    
    vf0=i*v2;
    ARIADNE_TEST_EQUAL(vf0,Vector<Float>("[-1.0,4.5]"));    
    vf0=v1/i;
    ARIADNE_TEST_EQUAL(vf0,Vector<Float>("[0.125,-0.75]"));    
  
    // Test Vector<Interval>
    Vector< Interval > iv1=Vector<Interval>("[[0.984375,1.015625],[2.25,2.375],[4.0,4.375],[-0.03125,0.015625]]");
    ARIADNE_TEST_EQUAL(iv1,Vector<Interval>("[[0.984375,1.015625],[2.25,2.375],[4.0,4.375],[-0.03125,0.015625]]"));
    ARIADNE_TEST_EQUAL(norm(iv1).lower(),4.0);
    ARIADNE_TEST_EQUAL(norm(iv1).upper(),4.375);

    Vector< Interval > iv2=Vector<Interval>("[[-1,1],[-1,1]]");
    ARIADNE_TEST_EQUAL(iv2,Vector<Interval>("[[-1,1],[-1,1]]"));
    Vector< Interval > iv3(3);
    ARIADNE_TEST_EQUAL(iv3,Vector<Interval>("[[0,0],[0,0],[0,0]]"));
    iv3=Vector<Interval>("[[4.25,4.25],[2.375,2.375]]");
    ARIADNE_TEST_EQUAL(iv3,Vector<Interval>("[[4.25,4.25],[2.375,2.375]]"));
    Interval ix=Interval(-2,1);
 
    Vector< Interval > iv0;
    iv1=iv0;
    ARIADNE_TEST_EQUAL(iv1,iv0);
    iv1=iv2;
    ARIADNE_TEST_EQUAL(iv1,iv2);
  

    iv1=iv2+iv3;
    ARIADNE_TEST_EQUAL(iv1,Vector<Interval>("[[3.25,5.25],[1.375,3.375]]"));  
    iv1=iv2-iv3;
    ARIADNE_TEST_EQUAL(iv1,Vector<Interval>("[[-5.25,-3.25],[-3.375,-1.375]]"));  
    iv1=ix*iv3;
    ARIADNE_TEST_EQUAL(iv1,Vector<Interval>("[[-8.5,4.25],[-4.75,2.375]]"));  
    iv1=iv2*ix;
    ARIADNE_TEST_EQUAL(iv1,Vector<Interval>("[[-2,2],[-2,2]]"));  
    ix=Interval(1,2);
    iv1=iv2/ix;
    ARIADNE_TEST_EQUAL(iv1,iv2);  

    ARIADNE_TEST_ASSERT(v1==Vector<Float>("[0.25,-1.5]"));
    iv0=iv1+v1;
    ARIADNE_TEST_EQUAL(iv0,Vector<Interval>("[[-0.75,1.25],[-2.5,-0.5]]"));
    iv0=v1+iv1;
    ARIADNE_TEST_EQUAL(iv0,Vector<Interval>("[[-0.75,1.25],[-2.5,-0.5]]"));
    iv0=iv1-v1;
    ARIADNE_TEST_EQUAL(iv0,Vector<Interval>("[[-1.25,0.75],[0.5,2.5]]"));
    iv0=v1-iv1;
    ARIADNE_TEST_EQUAL(iv0,Vector<Interval>("[[-0.75,1.25],[-2.5,-0.5]]"));
    iv0=x*iv1;
    ARIADNE_TEST_EQUAL(iv0,Vector<Interval>("[[-2,2],[-2,2]]"));
    iv0=ix*v1;
    ARIADNE_TEST_EQUAL(iv0,Vector<Interval>("[[0.25,0.5],[-3,-1.5]]"));
    iv0=iv1*x;
    ARIADNE_TEST_EQUAL(iv0,Vector<Interval>("[[-2,2],[-2,2]]"));
    iv0=v1*ix;
    ARIADNE_TEST_EQUAL(iv0,Vector<Interval>("[[0.25,0.5],[-3,-1.5]]"));
    iv0=iv1/x;
    ARIADNE_TEST_EQUAL(iv0,Vector<Interval>("[[-0.5,0.5],[-0.5,0.5]]"));
    iv0=v1/ix;
    ARIADNE_TEST_EQUAL(iv0,Vector<Interval>("[[0.125,0.25],[-1.5,-0.75]]"));
    // Test addition and multiplication with Integers
    iv1=Vector<Interval>("[[-0.25,0.25],[-3.0,1.5]]");
    iv0=iv1*i;
    ARIADNE_TEST_EQUAL(iv0,Vector<Interval>("[[-0.5,0.5],[-6.0,3.0]]"));
    iv1=Vector<Interval>("[[-2,2],[-0.75,1.25]]");
    iv0=i*iv1;
    ARIADNE_TEST_EQUAL(iv0,Vector<Interval>("[[-4.0,4.0],[-1.5,2.5]]"));
    iv0=iv0/i;
    ARIADNE_TEST_EQUAL(iv0,iv1);


/*
    iv0=v1;
    iv0/=ix;
    iv0=Vector<Interval>("[2,1]");
    iv1=Vector<Interval>("[0,1]");
    ARIADNE_TEST_EQUAL( (iv0+=iv1), Vector<Interval>("[2,2]") );
//    ARIADNE_TEST_ASSERT( (iv0-=Vector<Interval>("[0,1]")) == Vector<Interval>("[2,1]") );
//    ARIADNE_TEST_ASSERT( (iv0*=2) == Vector<Interval>("[4,2]") );
//    ARIADNE_TEST_ASSERT( (iv0/=4) == Vector<Interval>("[1,0.5]") );

*/

}


int main() {
    TestVector().test();

    return ARIADNE_TEST_FAILURES;
}  

