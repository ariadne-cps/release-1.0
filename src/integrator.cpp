/***************************************************************************
 *            integrator.cc
 *
 *  Copyright  2006-9  Pieter Collins
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

#include <iomanip>

#include "integrator.h"

#include "logging.h"
#include "container.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "taylor_function.h"
#include "box.h"

#include "polynomial.h"

namespace Ariadne {

typedef Vector<Interval> IVector;

Pair<Float,IVector>
IntegratorBase::flow_bounds(const VectorFunction& vf, const IVector& dx, const Float& hmax) const
{

    ARIADNE_ASSERT(vf.result_size()==dx.size());
    ARIADNE_ASSERT(hmax>0);

    // Set up constants of the method.
    const double INITIAL_MULTIPLIER=2;
    const double MULTIPLIER=1.125;
    const double BOX_RADIUS_MULTIPLIER=1.03125;
    const uint EXPANSION_STEPS=8;
    const uint REDUCTION_STEPS=8;
    const uint REFINEMENT_STEPS=4;

    IVector delta=(dx-midpoint(dx))*(BOX_RADIUS_MULTIPLIER-1);

    Float h=hmax;
    Float hmin=hmax/(1<<REDUCTION_STEPS);
    bool success=false;
    IVector bx,nbx,df;
    Interval ih(0,h);

    while(!success) {
        ARIADNE_ASSERT_MSG(h>hmin," h="<<h<<", hmin="<<hmin);
        bx=dx+INITIAL_MULTIPLIER*ih*vf.evaluate(dx)+delta;
        for(uint i=0; i!=EXPANSION_STEPS; ++i) {
            df=vf.evaluate(bx);
            nbx=dx+delta+ih*df;
            if (!Box(nbx).bounded()) {
                success = false;
                break;
            } else if(subset(nbx,bx)) {
                success=true;
                break;
            } else {
                bx=dx+delta+MULTIPLIER*ih*df;
            }
        }
        if(!success) {
            h/=2;
            ih=Interval(0,h);
        }
    }

    ARIADNE_ASSERT(subset(nbx,bx));

    IVector vfbx;
    vfbx=vf.evaluate(bx);

    for(uint i=0; i!=REFINEMENT_STEPS; ++i) {
        bx=nbx;
        vfbx=vf.evaluate(bx);
        nbx=dx+delta+ih*vfbx;
        ARIADNE_ASSERT_MSG(subset(nbx,bx),std::setprecision(20)<<"refinement "<<i<<": "<<nbx<<" is not a inside of "<<bx);
    }


    // Check result of operation
    // We use subset rather than inner subset here since the bound may touch
    ARIADNE_ASSERT(subset(nbx,bx));

    bx=nbx;

    ARIADNE_ASSERT(subset(dx,bx));

    ARIADNE_ASSERT_MSG(subset(dx+h*vf.evaluate(bx),bx),
        "d="<<dx<<"\nh="<<h<<"\nf(b)="<<vf.evaluate(bx)<<"\nd+hf(b)="<<Vector<Interval>(dx+h*vf.evaluate(bx))<<"\nb="<<bx<<"\n");

    return std::make_pair(h,bx);
}


VectorTaylorFunction
TaylorIntegrator::flow(const VectorFunction& f, const IVector& dx, const Float& hmax) const
{
    ARIADNE_LOG(1,"TaylorIntegrator::flow, f="<<f<<" dx="<<dx<<" hmax="<<hmax<<"\n");
    const uint nx=dx.size();

    IVector bx(nx);
    Float h;
    make_lpair(h,bx)=this->flow_bounds(f,dx,hmax);
    ARIADNE_LOG(2,"h="<<h<<" bx="<<bx<<"\n");

    IVector dom=join(dx,Interval(-h,h));
    ARIADNE_LOG(2,"dom="<<dom<<"\n");

    VectorTaylorFunction phi0(nx,ScalarTaylorFunction(dom));
    for(uint i=0; i!=nx; ++i) { phi0[i]=ScalarTaylorFunction::variable(dom,i); }
    ARIADNE_LOG(2,"phi0="<<phi0<<"\n");

    VectorTaylorFunction phi(nx,ScalarTaylorFunction(dom));
    for(uint i=0; i!=nx; ++i) { phi[i]=ScalarTaylorFunction::constant(dom,bx[i]); }

    ARIADNE_LOG(2,"phi="<<phi<<"\n");
    for(uint k=0; k!=this->_temporal_order; ++k) {
        VectorTaylorFunction fphi=compose(f,phi);
        ARIADNE_LOG(3,"fphi="<<fphi<<"\n");
        for(uint i=0; i!=nx; ++i) {
            phi[i]=antiderivative(fphi[i],nx)+phi0[i];
        }
        ARIADNE_LOG(3,"phi="<<phi<<"\n");
    }

    VectorTaylorFunction res(nx,ScalarTaylorFunction(dom));
    for(uint i=0; i!=nx; ++i) { res[i]=phi[i]; }
    ARIADNE_LOG(2,"res="<<res<<"\n");
    return res;

}


} // namespace Ariadne
