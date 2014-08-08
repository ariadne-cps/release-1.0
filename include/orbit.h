/***************************************************************************
 *            orbit.h
 *
 *  Copyright 2007  Pieter Collins
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
 
/*! \file orbit.h
 *  \brief Orbits of dynamic systems
 */

#ifndef ARIADNE_ORBIT_H
#define ARIADNE_ORBIT_H

#include <utility>
#include <iostream>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include "numeric.h"
#include "graphics_interface.h"
#include "taylor_set.h"

namespace Ariadne {

#ifdef DOXYGEN
//! \brief Class for storing evolution data.
template<class E> class Orbit {
  public:
    //! The type used to store a single enclosure set.
    typedef E EnclosureType;
    //! The type used to store a list of enclosure sets.
    typedef EL EnclosureListType;
    //! The initial set of the orbit.
    EnclosureType const& initial() const;
    //! The set of points reached by evolution from the given initial set over the evolution time.
    EnclosureListType const& reach() const;
    //! The set of points reached by evolution from the given initial set at the final evolution time.
    EnclosureListType const& final() const;
};
#endif

typedef double Time;

template<class ES> class Orbit;

template<class BS> class ListSet;
class Point;
class InterpolatedCurve;
class HybridDenotableSet;
class HybridTime;

class DiscreteLocation;
template<class BS> class HybridBasicSet;

typedef HybridBasicSet<Point> LocalisedPoint;
typedef HybridBasicSet<Box> LocalisedBox;
typedef HybridBasicSet<TaylorSet> LocalisedTaylorSet;
typedef HybridBasicSet<InterpolatedCurve> LocalisedInterpolatedCurve;
typedef ListSet<TaylorSet> TaylorSetList;
typedef ListSet<LocalisedTaylorSet> HybridTaylorSetList;

template<class ES> std::ostream& operator<<(std::ostream&, const Orbit<ES>&);

template<>
class Orbit<Point>
{
  public:
    Orbit(const Point& pt);
    void insert(Time t, const Point& hpt);
    const InterpolatedCurve& curve() const { return *this->_curve; }
  private:
    boost::shared_ptr< InterpolatedCurve > _curve;
};

template<>
class Orbit<LocalisedPoint>
{
  public:
    Orbit(const LocalisedPoint& hpt);
    void insert(HybridTime ht, LocalisedPoint& hpt);
    uint size() const;
    const InterpolatedCurve& curve(uint m) const;
    const std::vector<LocalisedInterpolatedCurve>& curves() const { return *this->_curves; }
  private:
    boost::shared_ptr<std::vector<LocalisedInterpolatedCurve> > _curves;
};

template<class ES>
class Orbit
{
    typedef ListSet<ES> ESL;
  public:
    typedef ES EnclosureType;
    typedef ListSet<ES> EnclosureListType;

    Orbit(const ES& set) : _initial(set) { }
    void adjoin_reach(const EnclosureType& set) { this->_reach.adjoin(set); }
    void adjoin_intermediate(const EnclosureType& set) { this->_intermediate.adjoin(set); }
    void adjoin_final(const EnclosureType& set) { this->_final.adjoin(set); }

    void adjoin_reach(const EnclosureListType& set) { this->_reach.adjoin(set); }
    void adjoin_intermediate(const EnclosureListType& set) { this->_intermediate.adjoin(set); }
    void adjoin_final(const EnclosureListType& set) { this->_final.adjoin(set); }

    EnclosureType const& initial() const { return this->_initial; }
    EnclosureListType const& reach() const { return this->_reach; }
    EnclosureListType const& intermediate() const { return this->_intermediate; }
    EnclosureListType const& final() const { return this->_final; }
  private:
    ES _initial;
    ESL _reach;
    ESL _intermediate;
    ESL _final;
};

template<>
class Orbit<TaylorSet>
{
    class Data;
    typedef TaylorSetList list_set_const_iterator;
  public:
    typedef TaylorSet EnclosureType;
    typedef TaylorSetList EnclosureListType;

    Orbit(const TaylorSet&);
    void adjoin_reach(const TaylorSet& set);
    void adjoin_intermediate(const TaylorSet& set);
    void adjoin_final(const TaylorSet& set);

    void adjoin_reach(const TaylorSetList& set);
    void adjoin_intermediate(const TaylorSetList& set);
    void adjoin_final(const TaylorSetList& set);

    TaylorSet const& initial() const;
    TaylorSetList const& reach() const;
    TaylorSetList const& intermediate() const;
    TaylorSetList const& final() const;
  private:
    boost::shared_ptr<Data> _data;
};

template<>
class Orbit<LocalisedTaylorSet>
{
    class Data;
    typedef HybridTaylorSetList list_set_const_iterator;
  public:
    typedef LocalisedTaylorSet EnclosureType;
    typedef HybridTaylorSetList EnclosureListType;

    Orbit(const LocalisedTaylorSet&);
    void adjoin_reach(const LocalisedTaylorSet& set);
    void adjoin_intermediate(const LocalisedTaylorSet& set);
    void adjoin_final(const LocalisedTaylorSet& set);

    void adjoin_reach(const HybridTaylorSetList& set);
    void adjoin_intermediate(const HybridTaylorSetList& set);
    void adjoin_final(const HybridTaylorSetList& set);

    LocalisedTaylorSet const& initial() const;
    HybridTaylorSetList const& reach() const;
    HybridTaylorSetList const& intermediate() const;
    HybridTaylorSetList const& final() const;

  private:
    boost::shared_ptr<Data> _data;
};

template<class ES> std::ostream& operator<<(std::ostream& os, const Orbit< ES >& orb);
template<> std::ostream& operator<<(std::ostream& os, const Orbit<TaylorSet>& orb);
template<> std::ostream& operator<<(std::ostream& os, const Orbit<LocalisedTaylorSet>& orb);

template<class ES>
std::ostream& 
operator<<(std::ostream& os, const Orbit< ES >& orb)
{
    os << "Orbit(\n  initial=" << orb.initial()
       << "\n  intermediate=" << orb.intermediate()
       << "\n  reach=" << orb.reach()
       << "\n  final=" << orb.final()
       << ")\n";
    return os;
}

template<> 
std::ostream& 
operator<<(std::ostream& os, const Orbit< LocalisedPoint >& orb);



} // namespace Ariadne

#endif // ARIADNE_ORBIT_H
