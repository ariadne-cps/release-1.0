/***************************************************************************
 *            set_splitting.cc
 *
 *  Copyright  2017  Luca Geretti
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

#include <cstdarg>
#include "ariadne.h"
#include "taylor_calculus.h"

using namespace Ariadne;

TaylorSet get_finishing_set(TaylorSet initial_set_model, Float step, VectorFunction dynamic) {

    TaylorCalculus tc(10u,10u,1e-12,1e-16);
    Box initial_set_bounding_box = initial_set_model.bounding_box();
    Box flow_bounds;
    make_lpair(step,flow_bounds)=tc.flow_bounds(dynamic,initial_set_bounding_box,step,1e3);
    TaylorCalculus::FlowModelType flow_model= tc.flow_model(dynamic,initial_set_bounding_box,step,flow_bounds);
    ScalarTaylorFunction identity_time_expression=ScalarTaylorFunction::variable(Box(1u,Interval(-step,+step)),0u);
    TaylorSet flow_set_model = unchecked_apply(flow_model,combine(initial_set_model.models(),identity_time_expression.model()));
    return partial_evaluate(flow_set_model.models(),initial_set_model.argument_size(),1.0);
}

int main(int argc, char* argv[])
{
    HybridIOAutomaton vanderpol("vanderpol");

    RealParameter mu("mu",1.0);

    RealVariable x("x"), y("y");
    vanderpol.add_internal_var(x);
    vanderpol.add_internal_var(y);

    RealExpression x_d = y;
    RealExpression y_d = mu * (1.0 - Ariadne::sqr(x))*y - x;

    DiscreteLocation loc("loc");
    vanderpol.new_mode(loc);
    vanderpol.set_dynamics(loc,x,x_d);
    vanderpol.set_dynamics(loc,y,y_d);

    Float eps = 1e-3;
    Box initial_box(2, 2.0-eps,2.0+eps, 0.0-eps,0.0+eps);
    LocalisedEnclosureType initial_enclosure(loc,initial_box);

    RealVectorFunction dynamic = vanderpol.dynamic_function(loc);
    TaylorSet initial_set_model(initial_box);
    Float step = 1e0;

    TaylorCalculus tc(6u,6u,1e-10,1e-16);
    Box initial_set_bounding_box = initial_set_model.bounding_box();
    Box flow_bounds;
    make_lpair(step,flow_bounds)=tc.flow_bounds(dynamic,initial_set_bounding_box,step,1e3);
    TaylorCalculus::FlowModelType flow_model= tc.flow_model(dynamic,initial_set_bounding_box,step,flow_bounds);
    ScalarTaylorFunction identity_time_expression=ScalarTaylorFunction::variable(Box(1u,Interval(-step,+step)),0u);

    TaylorSet flow_set_model_unsplit = unchecked_apply(flow_model,combine(initial_set_model.models(),identity_time_expression.model()));
    TaylorSet finishing_set_unsplit = partial_evaluate(flow_set_model_unsplit.models(),initial_set_model.argument_size(),1.0);

    std::cout << "Not splitted: " << finishing_set_unsplit.widths() - initial_set_model.widths() << std::endl;

    std::pair<TaylorSet,TaylorSet> split_set = initial_set_model.split();

    TaylorSet finishing_set_split1 = get_finishing_set(split_set.first,step,dynamic);
    TaylorSet finishing_set_split2 = get_finishing_set(split_set.second,step,dynamic);

    std::cout << "Splitted, dedicated flow model: " << finishing_set_split1.widths()-split_set.first.widths() << " and " << finishing_set_split2.widths()-split_set.second.widths() << std::endl;

    Box union_set = hull(finishing_set_split1.bounding_box(),finishing_set_split2.bounding_box());

    std::cout << "Joined, dedicated flow model: " << union_set.widths() - initial_set_model.widths() << std::endl;

    TaylorSet flow_set_model_shared_split1 = unchecked_apply(flow_model,combine(split_set.first.models(),identity_time_expression.model()));
    TaylorSet finishing_set_shared_split1 = partial_evaluate(flow_set_model_shared_split1.models(),initial_set_model.argument_size(),1.0);

    TaylorSet flow_set_model_shared_split2 = unchecked_apply(flow_model,combine(split_set.second.models(),identity_time_expression.model()));
    TaylorSet finishing_set_shared_split2 = partial_evaluate(flow_set_model_shared_split2.models(),initial_set_model.argument_size(),1.0);

    std::cout << "Splitted, shared flow model: " << finishing_set_shared_split1.widths()-split_set.first.widths() << " and " << finishing_set_shared_split2.widths()-split_set.second.widths() << std::endl;
    Box union_shared_set = hull(finishing_set_shared_split1.bounding_box(),finishing_set_shared_split2.bounding_box());

    std::cout << "Joined, shared flow model: " << union_shared_set.widths() - initial_set_model.widths() << std::endl;
}
