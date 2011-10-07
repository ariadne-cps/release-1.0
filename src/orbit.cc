/***************************************************************************
 *            orbit.cc
 *
 *  Copyright 2007-8  Pieter Collins
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
 
#include <utility>

#include "orbit.h"

#include "box.h"
#include "point.h"
#include "curve.h"
#include "taylor_set.h"
#include "list_set.h"
#include "hybrid_set.h"

#include "hybrid_time.h"

namespace Ariadne {


Orbit<Point>::Orbit(const Point& pt)
    : _curve(new InterpolatedCurve(0.0,pt))
{ }

void 
Orbit<Point>::insert(Time t, const Point& pt)
{
    this->_curve->insert(t,pt);
}


Orbit<LocalisedPoint>::Orbit(const LocalisedPoint& pt)
    : _curves(new std::vector<LocalisedInterpolatedCurve>(1u,make_pair(pt.first,InterpolatedCurve(pt.second))))
{ }

uint
Orbit<LocalisedPoint>::size() const
{
    return this->_curves->size();
}

const InterpolatedCurve& 
Orbit<LocalisedPoint>::curve(uint m) const
{
    return (*this->_curves)[m].second; 
}

void 
Orbit<LocalisedPoint>::insert(HybridTime ht, LocalisedPoint& hpt)
{
    ARIADNE_ASSERT((uint)ht.discrete_time()<=this->size());
    if(ht.discrete_time()==(int)this->size()) {
        this->_curves->push_back(make_pair(hpt.location(),InterpolatedCurve(hpt.continuous_state_set())));
    } else {
        (*this->_curves)[ht.discrete_time()].second.insert(ht.continuous_time(),hpt.continuous_state_set());
    }
}


template<> 
std::ostream& 
operator<<(std::ostream& os, const Orbit< LocalisedPoint >& orb)
{
    return os << orb.curves();
}


struct Orbit<TaylorSet>::Data {
    Data(const TaylorSet& initial_set) 
        : initial(initial_set) { }
    TaylorSet initial;
    TaylorSetList reach;
    TaylorSetList intermediate;
    TaylorSetList final;
};

Orbit<TaylorSet>::
Orbit(const TaylorSet& initial_set)
    : _data(new Data(initial_set))
{
}

void
Orbit<TaylorSet>::
adjoin_reach(const TaylorSet& set)
{
    this->_data->reach.adjoin(set);
}

void
Orbit<TaylorSet>::
adjoin_intermediate(const TaylorSet& set)
{
    this->_data->intermediate.adjoin(set);
}

void
Orbit<TaylorSet>::
adjoin_final(const TaylorSet& set)
{
    this->_data->final.adjoin(set);
}


void
Orbit<TaylorSet>::
adjoin_reach(const TaylorSetList& list_set)
{
    this->_data->reach.adjoin(list_set);
}

void
Orbit<TaylorSet>::
adjoin_intermediate(const TaylorSetList& list_set)
{
    this->_data->intermediate.adjoin(list_set);
}

void
Orbit<TaylorSet>::
adjoin_final(const TaylorSetList& list_set)
{
    this->_data->final.adjoin(list_set);
}


TaylorSet const&
Orbit<TaylorSet>::
initial() const
{
    return this->_data->initial;
}

TaylorSetList const&
Orbit<TaylorSet>::
reach() const
{
    return this->_data->reach;
}

TaylorSetList const&
Orbit<TaylorSet>::
intermediate() const
{
    return this->_data->intermediate;
}

TaylorSetList const&
Orbit<TaylorSet>::
final() const
{
    return this->_data->final;
}




struct Orbit<LocalisedTaylorSet>::Data {
    Data(const LocalisedTaylorSet& initial_set) 
        : initial(initial_set) { }
    LocalisedTaylorSet initial;
    HybridTaylorSetList reach;
    HybridTaylorSetList intermediate;
    HybridTaylorSetList final;
};

Orbit<LocalisedTaylorSet>::
Orbit(const LocalisedTaylorSet& initial_set)
    : _data(new Data(initial_set))
{
}

void
Orbit<LocalisedTaylorSet>::
adjoin_reach(const LocalisedTaylorSet& set)
{
    this->_data->reach.adjoin(set);
}

void
Orbit<LocalisedTaylorSet>::
adjoin_intermediate(const LocalisedTaylorSet& set)
{
    this->_data->intermediate.adjoin(set);
}

void
Orbit<LocalisedTaylorSet>::
adjoin_final(const LocalisedTaylorSet& set)
{
    this->_data->final.adjoin(set);
}


void
Orbit<LocalisedTaylorSet>::
adjoin_reach(const HybridTaylorSetList& list_set)
{
    this->_data->reach.adjoin(list_set);
}

void
Orbit<LocalisedTaylorSet>::
adjoin_intermediate(const HybridTaylorSetList& list_set)
{
    this->_data->intermediate.adjoin(list_set);
}

void
Orbit<LocalisedTaylorSet>::
adjoin_final(const HybridTaylorSetList& list_set)
{
    this->_data->final.adjoin(list_set);
}


LocalisedTaylorSet const&
Orbit<LocalisedTaylorSet>::
initial() const
{
    return this->_data->initial;
}

HybridTaylorSetList const&
Orbit<LocalisedTaylorSet>::
reach() const
{
    return this->_data->reach;
}

HybridTaylorSetList const&
Orbit<LocalisedTaylorSet>::
intermediate() const
{
    return this->_data->intermediate;
}

HybridTaylorSetList const&
Orbit<LocalisedTaylorSet>::
final() const
{
    return this->_data->final;
}


template<> 
std::ostream& 
operator<<(std::ostream& os, const Orbit<TaylorSet>& orb)
{
    os << "Orbit(\n  initial=" << Orbit<TaylorSet>::EnclosureListType(orb.initial()).bounding_boxes()
       << "\n  intermediate=" << orb.intermediate().bounding_boxes()
       << "\n  reach=" << orb.reach().bounding_boxes()
       << "\n  final=" << orb.final().bounding_boxes()
       << ")\n";
    return os;
}


template<> 
std::ostream& 
operator<<(std::ostream& os, const Orbit<LocalisedTaylorSet>& orb)
{
    os << "Orbit(\n  initial=" << Orbit<LocalisedTaylorSet>::EnclosureListType(orb.initial()).bounding_boxes()
       << "\n  intermediate=" << orb.intermediate().bounding_boxes()
       << "\n  reach=" << orb.reach().bounding_boxes()
       << "\n  final=" << orb.final().bounding_boxes()
       << ")\n";
    return os;
}

void draw(CanvasInterface& graphic, const Orbit<TaylorSet>& orbit)
{
    orbit.reach().draw(graphic);
    orbit.initial().draw(graphic);
    orbit.final().draw(graphic);
}



void draw(CanvasInterface& graphic, const Orbit<LocalisedPoint>& orbit)
{
    for(uint i=0; i<=orbit.size(); ++i) {
        orbit.curve(i).draw(graphic);
    }
}


} // namespace Ariadne
