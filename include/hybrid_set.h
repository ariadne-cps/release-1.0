/***************************************************************************
 *            hybrid_set.h
 *
 *  Copyright 2008-11  Pieter Collins, Alberto Casagrande
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

/*! \file hybrid_set.h
 *  \brief Sets in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_SET_H
#define ARIADNE_HYBRID_SET_H

#include <map>

#include <boost/iterator.hpp>
#include <boost/iterator_adaptors.hpp>

#include <boost/shared_ptr.hpp>

#include "macros.h"
#include "stlio.h"
#include "function_set.h"
#include "list_set.h"
#include "denotable_set.h"
#include "curve.h"

#include "hybrid_set_interface.h"
#include "point.h"
#include "box.h"
#include "orbit.h"

#include "serialization.h"

#include "graphics_interface.h"

namespace Ariadne {

class DiscreteLocation;

class HybridDenotableSet;
class HybridImageSet;
class HybridConstraintSet;

typedef std::map<DiscreteLocation,Vector<Float> > HybridFloatVector;
typedef std::map<DiscreteLocation,VectorFunction > HybridVectorFunction;

template<class HBS> class HybridBasicSetExpression { };
template<class HDS> class HybridDenotableSetExpression { };

//! \brief A hybrid space \f$\bigsqcup_{q\in Q} \R^{d_q}\f$ with discrete states \f$Q\f$.
class HybridSpace
    : public std::map<DiscreteLocation,uint>
{
  public:
    //! \brief The interface satisfied by bounded sets in the space.
    typedef HybridBoundedSetInterface BoundedSetInterfaceType;
    //! \brief The interface satisfied by overt sets in the space.
    typedef HybridOvertSetInterface OvertSetInterfaceType;
    //! \brief The interface satisfied by over sets in the space.
    typedef HybridOpenSetInterface OpenSetInterfaceType;
    //! \brief The interface satisfied by closed sets in the space.
    typedef HybridClosedSetInterface ClosedSetInterfaceType;
    //! \brief The interface satisfied by compact sets in the space.
    typedef HybridCompactSetInterface CompactSetInterfaceType;
    //! \brief The interface satisfied by regular sets in the space.
    typedef HybridRegularSetInterface RegularSetInterfaceType;
    //! \brief The interface satisfied by located sets in the space.
    typedef HybridLocatedSetInterface LocatedSetInterfaceType;
    //! \brief The type of approximations to sets in the space.
    typedef HybridDenotableSet SetApproximationType;

    typedef std::map<DiscreteLocation,uint>::const_iterator
    locations_const_iterator;

    HybridSpace() : std::map<DiscreteLocation,uint>() { }
    template<class SET> HybridSpace(const std::map<DiscreteLocation,SET>& qsmap) {
        for(typename std::map<DiscreteLocation,SET>::const_iterator loc_iter
                =qsmap.begin(); loc_iter!=qsmap.end(); ++loc_iter) {
            this->insert(make_pair(loc_iter->first,loc_iter->second.dimension())); }
    }
    template<class HSET> HybridSpace(const HSET& set) {
        for(typename HSET::locations_const_iterator loc_iter
                =set.locations_begin(); loc_iter!=set.locations_end(); ++loc_iter) {
            this->insert(make_pair(loc_iter->first,loc_iter->second.dimension())); }
    }

    locations_const_iterator locations_begin() const {
        return this->std::map<DiscreteLocation,uint>::begin(); }
    locations_const_iterator locations_end() const {
        return this->std::map<DiscreteLocation,uint>::end(); }
};



HybridBoxes
hull(const HybridBoxes& box1, const HybridBoxes& box2);

bool
superset(const HybridBoxes& box1, const HybridBoxes& box2);

bool
covers(const HybridBoxes& hboxes, const LocalisedBox& hbox);

HybridBoxes
shrink_in(const HybridBoxes& box, const HybridFloatVector& epsilon);

HybridBoxes
shrink_out(const HybridBoxes& box, const HybridFloatVector& epsilon);

HybridBoxes
widen(const HybridBoxes& box);

HybridBoxes
unbounded_hybrid_boxes(const HybridSpace& hspace);

Box
project(const HybridBoxes& box, const std::vector<uint>& dimensions);

HybridBoxes
project(const Box& box, const std::vector<uint>& dimensions, const HybridSpace& target_space);

inline
HybridBoxes
bounding_boxes(const std::map<DiscreteLocation,uint> space, Interval bound)
{
    HybridBoxes result;
    for(std::map<DiscreteLocation,uint>::const_iterator loc_iter=space.begin();
        loc_iter!=space.end(); ++loc_iter)
        {
            result.insert(make_pair(loc_iter->first,Box(loc_iter->second, bound)));
        }
    return result;
}

inline
HybridBoxes
bounding_boxes(const std::map<DiscreteLocation,uint> space, const Box& bbox)
{
    HybridBoxes result;
    for(std::map<DiscreteLocation,uint>::const_iterator loc_iter=space.begin();
        loc_iter!=space.end(); ++loc_iter)
        {
            result.insert(make_pair(loc_iter->first,bbox));
        }
    return result;
}

inline
HybridSpace
state_space(const HybridBoxes& hboxes)
{
	HybridSpace result;

	for (HybridBoxes::const_iterator hb_it = hboxes.begin(); hb_it != hboxes.end(); ++hb_it) {
		result.insert(make_pair(hb_it->first,hb_it->second.dimension()));
	}

	return result;
}


template<class BS>
class HybridBasicSet
    : public std::pair<DiscreteLocation,BS>
{
  public:
    typedef BS ContinuousStateSetType;
    HybridBasicSet() : std::pair<DiscreteLocation,BS>() { }
    HybridBasicSet(const DiscreteLocation& q, const BS& s) : std::pair<DiscreteLocation,BS>(q,s) { }
    HybridBasicSet(const std::pair<DiscreteLocation,BS>& p) : std::pair<DiscreteLocation,BS>(p) { }
    const DiscreteLocation& location() const { return this->first; }
    const ContinuousStateSetType& continuous_state_set() const { return this->second; }
};

template< class DS, class HBS >
class HybridSetConstIterator
    : public boost::iterator_facade<HybridSetConstIterator<DS,HBS>,
                                    HBS,
                                    boost::forward_traversal_tag,
                                    HBS const&
                                    >

{
  public:
    typedef HBS const& reference;
  public:
    HybridSetConstIterator(const std::map<DiscreteLocation,DS>&, bool);
    bool equal(const HybridSetConstIterator<DS,HBS>&) const;
    const HBS& dereference() const;
    void increment();
  private:
    void increment_loc();
  private:
    typename std::map< DiscreteLocation,DS>::const_iterator loc_begin;
    typename std::map< DiscreteLocation,DS>::const_iterator loc_end;
    typename std::map< DiscreteLocation,DS>::const_iterator loc_iter;
    typename DS::const_iterator bs_iter;
    mutable HBS hybrid_set;
};


//! A set comprising an ImageSet in each location.
class HybridImageSet
    : public std::map<DiscreteLocation,ImageSet>
    , public HybridLocatedSetInterface
{
  public:
    typedef std::map<DiscreteLocation,ImageSet>::iterator locations_iterator;
    typedef std::map<DiscreteLocation,ImageSet>::const_iterator locations_const_iterator;
    locations_iterator locations_begin() {
        return this->std::map<DiscreteLocation,ImageSet>::begin(); }
    locations_iterator locations_end() {
        return this->std::map<DiscreteLocation,ImageSet>::end(); }
    locations_const_iterator locations_begin() const {
        return this->std::map<DiscreteLocation,ImageSet>::begin(); }
    locations_const_iterator locations_end() const {
        return this->std::map<DiscreteLocation,ImageSet>::end(); }

    using std::map<DiscreteLocation,ImageSet>::insert;
    
    virtual HybridImageSet* clone() const { return new HybridImageSet(*this); }
    virtual HybridSpace space() const { return HybridSpace(*this); }
    virtual ImageSet& operator[](DiscreteLocation q) {
        return this->std::map<DiscreteLocation,ImageSet>::operator[](q); }
    virtual ImageSet const& operator[](DiscreteLocation q) const {
        ARIADNE_ASSERT(this->find(q)!=this->locations_end());
        return this->find(q)->second; }
    virtual tribool overlaps(const LocalisedBox& lbx) const {
        locations_const_iterator loc_iter=this->find(lbx.first);
        if (loc_iter!=this->locations_end()) return loc_iter->second.overlaps(lbx.second);
		else return false; }
    virtual tribool disjoint(const LocalisedBox& lbx) const {
        locations_const_iterator loc_iter=this->find(lbx.first);
        if (loc_iter!=this->locations_end()) return loc_iter->second.disjoint(lbx.second);
		else return true; } 
    virtual tribool inside(const HybridBoxes& hbx) const  {
		tribool result = true; // Initially assumed as true
        for(locations_const_iterator loc_iter=this->begin(); loc_iter!=this->locations_end(); ++loc_iter) {
            if(!loc_iter->second.empty()) {
                HybridBoxes::const_iterator hbx_loc_iter=hbx.find(loc_iter->first);
                if(hbx_loc_iter!=hbx.end()) // If the locations match
				{
					tribool temp_result = loc_iter->second.inside(hbx_loc_iter->second); 
					if (!possibly(temp_result)) return false; // If the ImageSet is not inside the box, returns false
					else if (indeterminate(temp_result)) result = indeterminate; // If it is not decided, the result is reduced to indeterminate
				}
				else return false; // Otherwise the result is certainly false
            } 
		} 
		return result; }
    virtual HybridBoxes bounding_box() const {
        HybridBoxes result;
        for(locations_const_iterator loc_iter=this->begin(); loc_iter!=this->locations_end(); ++loc_iter) {
            if(!loc_iter->second.empty()) { result.insert(std::make_pair(loc_iter->first,loc_iter->second.bounding_box())); } }
        return result; }
    virtual std::ostream& write(std::ostream& os) const { 
        os << "HybridImageSet(";
        for(locations_const_iterator loc_iter=this->begin(); loc_iter!=this->locations_end(); ++loc_iter) {
          if (loc_iter!=this->begin()) os << ", ";
          os << loc_iter->first << ":" << loc_iter->second;
        }
        os << ")";
        return os;
    }
};


//! A set comprising a ConstraintSet in each location.
//! The semantics for an absent location is that no constraint is present (i.e., satisfaction returns true for each point of a set)
class HybridConstraintSet
    : public std::map<DiscreteLocation,ConstraintSet>
    , public HybridRegularSetInterface
{
  public:
	typedef std::map<DiscreteLocation,ConstraintSet>::iterator locations_iterator;
	typedef std::map<DiscreteLocation,ConstraintSet>::const_iterator locations_const_iterator;
	locations_iterator locations_begin() {
		return this->std::map<DiscreteLocation,ConstraintSet>::begin(); }
	locations_iterator locations_end() {
		return this->std::map<DiscreteLocation,ConstraintSet>::end(); }
	locations_const_iterator locations_begin() const {
		return this->std::map<DiscreteLocation,ConstraintSet>::begin(); }
	locations_const_iterator locations_end() const {
		return this->std::map<DiscreteLocation,ConstraintSet>::end(); }

	using std::map<DiscreteLocation,ConstraintSet>::insert;

	HybridConstraintSet() { }
	//! \brief Constructs from a hybrid \a codomain and some functions \a func.
	HybridConstraintSet(const HybridVectorFunction& func, const HybridBoxes& codomain);

	/** \brief Constructs from a \a constraint, copied in all locations from \a hspace */
	HybridConstraintSet(const HybridSpace& hspace, ConstraintSet constraint);

	//! \brief Constructs a full-space constraint set in the given \a hspace.
	HybridConstraintSet(const HybridSpace& hspace);

    virtual HybridConstraintSet* clone() const { return new HybridConstraintSet(*this); }
    virtual HybridSpace space() const { return HybridSpace(*this); }
    virtual ConstraintSet& operator[](DiscreteLocation q) {
    	return this->std::map<DiscreteLocation,ConstraintSet>::operator[](q); }
    virtual ConstraintSet const& operator[](DiscreteLocation q) const {
        ARIADNE_ASSERT(this->find(q)!=this->locations_end());
        return this->find(q)->second; }

    virtual HybridVectorFunction functions() const {
    	HybridVectorFunction result;
    	for (HybridConstraintSet::const_iterator loc_it = this->begin(); loc_it != this->end(); ++loc_it)
    		result.insert(make_pair(loc_it->first,loc_it->second.function()));
    	return result;
    }

    virtual HybridBoxes codomain() const {
    	HybridBoxes result;
    	for (HybridConstraintSet::const_iterator loc_it = this->begin(); loc_it != this->end(); ++loc_it)
    		result.insert(make_pair(loc_it->first,loc_it->second.codomain()));
    	return result;
    }

    virtual tribool overlaps(const LocalisedBox& hbx) const {
        locations_const_iterator loc_iter=this->find(hbx.first);
        if (loc_iter!=this->locations_end()) return loc_iter->second.overlaps(hbx.second);
		else return false; }
    virtual tribool disjoint(const LocalisedBox& hbx) const {
        locations_const_iterator loc_iter=this->find(hbx.first);
        if (loc_iter!=this->locations_end()) return loc_iter->second.disjoint(hbx.second);
		else return true; }
    virtual tribool covers(const LocalisedBox& hbx) const {
        locations_const_iterator loc_iter=this->find(hbx.first);
        if (loc_iter!=this->locations_end()) return loc_iter->second.covers(hbx.second);
		else return false; }

    virtual std::ostream& write(std::ostream& os) const {
        os << "HybridConstraintSet(";
        for(locations_const_iterator loc_iter=this->begin(); loc_iter!=this->locations_end(); ++loc_iter) {
          if (loc_iter!=this->begin()) os << ", ";
          os << loc_iter->first << ":" << loc_iter->second;
        }
        os << ")";
        return os;
    }

};


//! A set comprising a BoundedConstraintSet in each location.
//! The semantics for an absent location is that the domain is empty, or equivalently that no set is present
//! (satisfaction of the constraint would yield false for all points of another set)
class HybridBoundedConstraintSet
	: public std::map<DiscreteLocation,BoundedConstraintSet>
	, public HybridLocatedSetInterface
{
  public:
	typedef std::map<DiscreteLocation,BoundedConstraintSet>::iterator locations_iterator;
	typedef std::map<DiscreteLocation,BoundedConstraintSet>::const_iterator locations_const_iterator;
	locations_iterator locations_begin() {
		return this->std::map<DiscreteLocation,BoundedConstraintSet>::begin(); }
	locations_iterator locations_end() {
		return this->std::map<DiscreteLocation,BoundedConstraintSet>::end(); }
	locations_const_iterator locations_begin() const {
		return this->std::map<DiscreteLocation,BoundedConstraintSet>::begin(); }
	locations_const_iterator locations_end() const {
		return this->std::map<DiscreteLocation,BoundedConstraintSet>::end(); }

	using std::map<DiscreteLocation,BoundedConstraintSet>::insert;

	HybridBoundedConstraintSet() { }
	//! \brief Constructs from a \a domain, \a codomain and some functions \a func.
	HybridBoundedConstraintSet(
			const HybridBoxes& domain,
			const HybridVectorFunction& func,
			const HybridBoxes& codomain);

	/** \brief Constructs from a \a constraint, copied in all locations from \a hspace */
	HybridBoundedConstraintSet(
			const HybridSpace& hspace,
			BoundedConstraintSet constraint);

	//! \brief Constructs from domain boxes.
	HybridBoundedConstraintSet(const HybridBoxes& domain_boxes);

	//! \brief Constructs an empty set, where all domains are empty.
	HybridBoundedConstraintSet(const HybridSpace& hspace);

	//! \brief Constructs from a single box extended to the \a hspace.
	HybridBoundedConstraintSet(
			const HybridSpace& hspace,
			const Box& domain_box);

	virtual HybridBoundedConstraintSet* clone() const { return new HybridBoundedConstraintSet(*this); }
	virtual HybridSpace space() const { return HybridSpace(*this); }
	virtual BoundedConstraintSet& operator[](DiscreteLocation q) {
		return this->std::map<DiscreteLocation,BoundedConstraintSet>::operator[](q); }
	virtual BoundedConstraintSet const& operator[](DiscreteLocation q) const {
		ARIADNE_ASSERT(this->find(q)!=this->locations_end());
		return this->find(q)->second; }

    virtual HybridBoxes domain() const {
    	HybridBoxes result;
    	for (HybridBoundedConstraintSet::const_iterator loc_it = this->begin(); loc_it != this->end(); ++loc_it)
    		result.insert(make_pair(loc_it->first,loc_it->second.domain()));
    	return result;
    }

	virtual HybridVectorFunction functions() const {
		HybridVectorFunction result;
		for (HybridBoundedConstraintSet::const_iterator loc_it = this->begin(); loc_it != this->end(); ++loc_it)
			result.insert(std::pair<DiscreteLocation,VectorFunction>(loc_it->first,loc_it->second.function()));
		return result;
	}

    virtual HybridBoxes codomain() const {
    	HybridBoxes result;
    	for (HybridBoundedConstraintSet::const_iterator loc_it = this->begin(); loc_it != this->end(); ++loc_it)
    		result.insert(make_pair(loc_it->first,loc_it->second.codomain()));
    	return result;
    }

	virtual tribool overlaps(const LocalisedBox& hbx) const {
		locations_const_iterator loc_iter=this->find(hbx.first);
		if (loc_iter!=this->locations_end()) return loc_iter->second.overlaps(hbx.second);
		else return false; }
	virtual tribool disjoint(const LocalisedBox& hbx) const {
		locations_const_iterator loc_iter=this->find(hbx.first);
		if (loc_iter!=this->locations_end()) return loc_iter->second.disjoint(hbx.second);
		else return true; }
	virtual tribool covers(const LocalisedBox& hbx) const {
		locations_const_iterator loc_iter=this->find(hbx.first);
		if (loc_iter!=this->locations_end()) return loc_iter->second.covers(hbx.second);
		else return false; }
    virtual tribool inside(const HybridBoxes& hbx) const  {
		tribool result = true;
        for(HybridBoxes::const_iterator loc_iter=hbx.begin(); loc_iter!=hbx.end(); ++loc_iter) {
        	locations_const_iterator this_iter = this->find(loc_iter->first);
        	if (this_iter == this->locations_end())
        		return false;
        	else {
        		tribool current_result = possibly(this_iter->second.inside(loc_iter->second));
        		if (!possibly(current_result))
        			return false;
        		else if (indeterminate(current_result))
        			result = indeterminate;
        	}
		}
		return result; }

	virtual HybridBoxes bounding_box() const {
		HybridBoxes result;
		for (HybridBoundedConstraintSet::const_iterator loc_it = this->begin(); loc_it != this->end(); ++loc_it)
			result.insert(std::pair<DiscreteLocation,Box>(loc_it->first,loc_it->second.bounding_box()));
		return result;
	}

	virtual std::ostream& write(std::ostream& os) const {
		os << "HybridBoundedConstraintSet(";
		for(locations_const_iterator loc_iter=this->begin(); loc_iter!=this->locations_end(); ++loc_iter) {
		  if (loc_iter!=this->begin()) os << ", ";
		  os << loc_iter->first << ":" << loc_iter->second;
		}
		os << ")";
		return os;
	}

};


//! A set comprising a %ListSet in each location.
template<class ES>
class HybridListSet
    : public std::map<DiscreteLocation,ListSet<ES> >
{
  public:
    typedef typename std::map< DiscreteLocation,ListSet<ES> >::iterator locations_iterator;
    typedef typename std::map< DiscreteLocation,ListSet<ES> >::const_iterator locations_const_iterator;
    typedef HybridSetConstIterator< ListSet<ES>, std::pair<DiscreteLocation,ES> > const_iterator;

    HybridListSet() { }
    HybridListSet(const std::pair<DiscreteLocation,ES>& hes) { this->adjoin(hes); }

    locations_iterator locations_begin() {
        return this->std::map<DiscreteLocation,ListSet<ES> >::begin(); }
    locations_iterator locations_end() {
        return this->std::map<DiscreteLocation,ListSet<ES> >::end(); }
    locations_const_iterator locations_begin() const {
        return this->std::map<DiscreteLocation,ListSet<ES> >::begin(); }
    locations_const_iterator locations_end() const {
        return this->std::map<DiscreteLocation,ListSet<ES> >::end(); }
    const_iterator begin() const {
        return const_iterator(*this,false); }
    const_iterator end() const {
        return const_iterator(*this,true); }

    /*! \brief Returns the number of basic hybrid sets forming this object. */
    size_t size() const {
        size_t s = 0;
        for(locations_const_iterator loc_iter=this->locations_begin();
            loc_iter!=this->locations_end(); ++loc_iter) {
            s += loc_iter->second.size();
        }
        return s;
    }

    //using std::map< DiscreteLocation,ListSet<ES> >::insert;

    //using std::map<DiscreteLocation,ListSet<ES> >::operator[];
    ListSet<ES>& operator[](const DiscreteLocation& q) {
        return this->std::map<DiscreteLocation,ListSet<ES> >::operator[](q); }
    const ListSet<ES>& operator[](const DiscreteLocation& q) const {
        ARIADNE_ASSERT_MSG(this->find(q)!=this->locations_end(),(*this)<<" has no location "<<q);
        return this->find(q)->second; }

    void adjoin(const DiscreteLocation& q, const ES& es) {
        (*this)[q].adjoin(es); }
    void adjoin(const std::pair<DiscreteLocation,ES>& hes) {
        (*this)[hes.first].adjoin(hes.second); }
    void adjoin(const HybridListSet<ES>& hls) {
        for(locations_const_iterator loc_iter=hls.locations_begin();
            loc_iter!=hls.locations_end(); ++loc_iter) {
            (*this)[loc_iter->first].adjoin(loc_iter->second); } }
    void adjoin(const ListSet<std::pair<DiscreteLocation,ES> >& hls) {
        for(locations_const_iterator loc_iter=hls.locations_begin();
            loc_iter!=hls.locations_end(); ++loc_iter) {
            (*this)[loc_iter->first].adjoin(loc_iter->second); } }
    void adjoin(const HybridBasicSet<ES>& hbs) {
        (*this)[hbs.location()].adjoin(hbs.continuous_state_set()); }

    HybridListSet<Box> bounding_boxes() const {
        HybridListSet<Box> result;
        for(locations_const_iterator loc_iter=this->locations_begin();
            loc_iter!=this->locations_end(); ++loc_iter) {
            result[loc_iter->first]=loc_iter->second.bounding_boxes(); }
        return result; }

    // A box for the whole set
    Box bounding_box() const {
        if(this->size()==0) { return Box(0); }
        locations_const_iterator loc_iter=this->locations_begin();
        Box result(loc_iter->second.bounding_box());
        for(loc_iter++; loc_iter!=this->locations_end(); ++loc_iter) {
            Box bx=loc_iter->second.bounding_box();
            result=hull(result,bx);             
        }
        return result;
    }


    HybridSpace space() const { return HybridSpace(*this); }
};


template<class ES>
class ListSet< HybridBasicSet<ES> >
    : public HybridListSet<ES>
{
  public:
    ListSet() { }
    ListSet(const HybridBasicSet<ES>& hes) { this->adjoin(hes); }
};


template<class ES>
std::ostream&
operator<<(std::ostream& os,
           const ListSet< HybridBasicSet<ES> >& ls) {
    const std::map< DiscreteLocation, ListSet<ES> >& hls=ls;
    return os << "HybridListSet" << hls;
}


class HybridGrid
    : public std::map<DiscreteLocation,Grid>
{
  public:
    typedef std::map<DiscreteLocation,Grid>::const_iterator
    locations_const_iterator;

    HybridGrid() { }

    HybridGrid(const HybridSpace& hspc, const Float l=1.0) {
        for(HybridSpace::locations_const_iterator loc_iter=hspc.begin();
            loc_iter!=hspc.end(); ++loc_iter) {
            this->insert(make_pair(loc_iter->first,Grid(loc_iter->second,l)));
        }
    }

    HybridGrid(const HybridSpace& hspc, const Grid& grid) {
        for(HybridSpace::locations_const_iterator loc_iter=hspc.begin();
            loc_iter!=hspc.end(); ++loc_iter) {
            this->insert(make_pair(loc_iter->first,Grid(grid)));
        }
    }

    template<class HGSET> HybridGrid(const HGSET& set) {
        for(typename HGSET::locations_const_iterator loc_iter=set.locations_begin();
        		loc_iter!=set.locations_end(); ++loc_iter) {
            this->insert(make_pair(loc_iter->first,loc_iter->second.grid())); 
        }
    }

    HybridSpace state_space() const {
        return HybridSpace(*this);
    }

    locations_const_iterator locations_begin() const {
        return this->std::map<DiscreteLocation,Grid>::begin(); }
    locations_const_iterator locations_end() const {
        return this->std::map<DiscreteLocation,Grid>::end(); }
    const Grid& operator[](DiscreteLocation q) const {
        ARIADNE_ASSERT(this->find(q)!=this->locations_end());
        return this->find(q)->second; }
    Grid& operator[](DiscreteLocation q) {
        return this->std::map<DiscreteLocation,Grid>::operator[](q); }
    std::map<DiscreteLocation,Vector<Float> > lengths() const;
};


class HybridDenotableSet;


template<class HDS1, class HDS2>
void adjoin_denotable_set(HDS1& hds1, const HDS2 hds2) {
    for(typename HDS2::locations_const_iterator loc2_iter=
            hds2.locations_begin(); loc2_iter!=hds2.locations_end(); ++loc2_iter)
        {
            typename HDS1::locations_iterator loc1_iter=hds1.find(loc2_iter->first);
            if(loc1_iter==hds1.locations_end()) {
                hds1.insert(make_pair(loc2_iter->first,
                                      typename HDS2::value_type::second_type(loc2_iter->second)));
            } else {
                loc1_iter->second.adjoin(loc2_iter->second);
            }
        }
}


//! A set comprising a denotable set in each location.
class HybridDenotableSet
    : public std::map<DiscreteLocation,DenotableSetType>
{
  public:
    typedef std::map<DiscreteLocation,DenotableSetType>::iterator locations_iterator;
    typedef std::map<DiscreteLocation,DenotableSetType>::const_iterator locations_const_iterator;
    typedef HybridSetConstIterator<DenotableSetType,LocalisedBox> const_iterator;
  public:
    locations_iterator locations_begin() {
        return this->std::map<DiscreteLocation,DenotableSetType>::begin(); }
    locations_iterator locations_end() {
        return this->std::map<DiscreteLocation,DenotableSetType>::end(); }
    locations_const_iterator locations_begin() const {
        return this->std::map<DiscreteLocation,DenotableSetType>::begin(); }
    locations_const_iterator locations_end() const {
        return this->std::map<DiscreteLocation,DenotableSetType>::end(); }
    const_iterator begin() const {
        return const_iterator(*this,false); }
    const_iterator  end() const {
        return const_iterator(*this,true); }
  public:
    HybridDenotableSet() { }

    virtual ~HybridDenotableSet() {}

    HybridDenotableSet(const HybridSpace& hspace) {
        for(HybridSpace::locations_const_iterator loc_iter = hspace.
                locations_begin(); loc_iter!=hspace.locations_end(); ++loc_iter) {
            this->insert(make_pair(loc_iter->first,Grid(loc_iter->second))); } }
    HybridDenotableSet(const HybridSpace& hspace, const Vector<Float>& lengths) {
        for(HybridSpace::locations_const_iterator loc_iter = hspace.
                locations_begin(); loc_iter!=hspace.locations_end(); ++loc_iter) {
            this->insert(make_pair(loc_iter->first,Grid(lengths))); } }
    HybridDenotableSet(const HybridGrid& hgrid) {
        for(HybridGrid::locations_const_iterator loc_iter = hgrid.
                locations_begin(); loc_iter!=hgrid.locations_end(); ++loc_iter) {
            this->insert(make_pair(loc_iter->first,Grid(loc_iter->second))); } }

    void adjoin(const HybridDenotableSet& hgts) {
        for(HybridDenotableSet::locations_const_iterator loc_iter=hgts.locations_begin(); loc_iter!=hgts.locations_end(); ++loc_iter) {
            if(!this->has_location(loc_iter->first)) {
                DenotableSetType new_gts(loc_iter->second.grid());
                this->insert(make_pair(loc_iter->first,new_gts)); }
            this->find(loc_iter->first)->second.adjoin(loc_iter->second); } }

    void remove(const HybridDenotableSet& hgts) {
        for(HybridDenotableSet::locations_const_iterator loc_iter=hgts.locations_begin(); loc_iter!=hgts.locations_end(); ++loc_iter) {
            if(this->has_location(loc_iter->first)) {
                this->find(loc_iter->first)->second.remove(loc_iter->second); } } }

    void restrict(const HybridDenotableSet& hgts) {
        for(HybridDenotableSet::locations_const_iterator loc_iter=hgts.locations_begin(); loc_iter!=hgts.locations_end(); ++loc_iter) {
            if(this->has_location(loc_iter->first)) {
                this->find(loc_iter->first)->second.restrict(loc_iter->second); } } }

    void adjoin_inner_approximation(const HybridBoxes& hbxs, const int depth) {
        for(HybridBoxes::const_iterator loc_iter=hbxs.begin();
            loc_iter!=hbxs.end(); ++loc_iter)
            {
                DiscreteLocation loc=loc_iter->first;
                if(!this->has_location(loc)) {
                    this->insert(make_pair(loc,DenotableSetType(loc_iter->second.dimension()))); }
                (*this)[loc].adjoin_inner_approximation(loc_iter->second,loc_iter->second,depth); } }

    void adjoin_lower_approximation(const HybridOvertSetInterface& hs, const int height, const int depth) {
        HybridSpace hspc=hs.space();
        for(HybridSpace::const_iterator loc_iter=hspc.begin();
            loc_iter!=hspc.end(); ++loc_iter)
            {
                DiscreteLocation loc=loc_iter->first;
                if(!this->has_location(loc)) {
                    this->insert(make_pair(loc,DenotableSetType(loc_iter->second))); }
                (*this)[loc].adjoin_lower_approximation(hs[loc],height,depth); } }

    void adjoin_outer_approximation(const HybridCompactSetInterface& hs, const int depth) {
        HybridSpace hspc=hs.space();
        for(HybridSpace::const_iterator loc_iter=hspc.begin();
            loc_iter!=hspc.end(); ++loc_iter)
            {
                DiscreteLocation loc=loc_iter->first;
                if(!this->has_location(loc)) {
                    this->insert(make_pair(loc,DenotableSetType(loc_iter->second))); }
                (*this)[loc].adjoin_outer_approximation(hs[loc],depth); } }

    void adjoin_outer_approximation(const HybridBoxes& hbxs, const int depth) {
        for(HybridBoxes::const_iterator loc_iter=hbxs.begin();
            loc_iter!=hbxs.end(); ++loc_iter)
            {
                DiscreteLocation loc=loc_iter->first;
                if(!this->has_location(loc)) {
                    this->insert(make_pair(loc,DenotableSetType(loc_iter->second.dimension()))); }
                (*this)[loc].adjoin_outer_approximation(loc_iter->second,depth); } }

    virtual HybridGrid grid() const { return HybridGrid(*this); }

    virtual bool has_location(DiscreteLocation q) const {
        return this->find(q)!=this->std::map<DiscreteLocation,DenotableSetType>::end(); }

    virtual DenotableSetType& operator[](DiscreteLocation q) {
        ARIADNE_ASSERT(this->has_location(q));
        return this->find(q)->second;
    }

    virtual const DenotableSetType& operator[](DiscreteLocation q) const {
        ARIADNE_ASSERT(this->has_location(q));
        return this->find(q)->second;
    }

    virtual size_t size() const {
        size_t result=0;
        for(locations_const_iterator loc_iter=this->locations_begin();
            loc_iter!=this->locations_end(); ++loc_iter) {
            result+=loc_iter->second.size(); }
        return result; }

    virtual void clear() {
    	for(locations_iterator loc_iter=this->locations_begin(); loc_iter!=this->locations_end(); ++loc_iter) {
    		loc_iter->second.clear();
    	}
    }

    bool empty() const {
        for(locations_const_iterator loc_iter=this->locations_begin();
            loc_iter!=this->locations_end(); ++loc_iter) {
            if(!loc_iter->second.empty()) { return false; } }
        return true; }


    HybridListSet<Box> boxes() const {
        HybridListSet<Box> result;
        for(const_iterator iter=this->begin();
            iter!=this->end(); ++iter) {
            result[iter->first].adjoin(iter->second); }
        return result; }
    void mince(int depth) {
        for(locations_iterator loc_iter=this->locations_begin();
            loc_iter!=this->locations_end(); ++loc_iter) {
            loc_iter->second.mince(depth); } }
    void recombine() {
        for(locations_iterator loc_iter=this->locations_begin();
            loc_iter!=this->locations_end(); ++loc_iter) {
            loc_iter->second.recombine(); } }
  public:
    // HybridSetInterface methods
    HybridDenotableSet* clone() const { return new HybridDenotableSet(*this); }
    HybridSpace space() const { return HybridSpace(*this); }

    tribool disjoint(const LocalisedBox& hbx) const {
        locations_const_iterator loc_iter = this->find( hbx.first );
        if (loc_iter != this->locations_end())
			return loc_iter->second.disjoint( hbx.second );
		else
			return true;
    }

    tribool overlaps(const LocalisedBox& hbx) const {
        locations_const_iterator loc_iter = this->find( hbx.first );
        if (loc_iter != this->locations_end())
			return loc_iter->second.overlaps( hbx.second );
		else
			return false;
    }

    tribool superset(const LocalisedBox& hbx) const {
        locations_const_iterator loc_iter=this->find(hbx.first);
        if (loc_iter != this->locations_end())
			return loc_iter->second.superset( hbx.second );
		else
			return false;
    }

    tribool subset(const HybridBoxes& hbx) const {
		bool result = true;
        for( locations_const_iterator loc_iter = this->locations_begin(); loc_iter != this->locations_end(); ++loc_iter ) {
		    if( !loc_iter->second.empty() ) {
                HybridBoxes::const_iterator hbx_loc_iter = hbx.find( loc_iter->first );
                if( hbx_loc_iter != hbx.end()) {
					const tribool temp_result = loc_iter->second.subset( hbx_loc_iter->second ); // Temporarily store the result
					if (!possibly(temp_result)) return false; // If the cells are not included in the box of the location, directly return false
					else if (indeterminate(temp_result)) result = indeterminate; // Otherwise if indeterminate, set the result as indeterminate
                }
				else // Otherwise it cannot be a subset of the hybrid box
					return false;
            }
        }
        return result;
    }

    HybridBoxes bounding_box() const {
        HybridBoxes result;
        for( locations_const_iterator loc_iter = this->locations_begin(); loc_iter != this->locations_end(); ++loc_iter ) {
            if( !loc_iter->second.empty() ) {
                result.insert( std::make_pair( loc_iter->first, loc_iter->second.bounding_box() ) );
            }
        }
        return result;
    }

    std::ostream& write(std::ostream& os) const {
        return os << static_cast<const std::map<DiscreteLocation,DenotableSetType>&>(*this);
    }

  private:
    /*
      friend class boost::serialization::access;
      template<class Archive> void serialize(Archive & ar, const unsigned int version) {
      ar & static_cast<std::map<int,DenotableSetType>&>(*this); }
    */
};



template<class DS, class HBS> inline
HybridSetConstIterator<DS,HBS>::
HybridSetConstIterator(const std::map<DiscreteLocation,DS>& map, bool end)
    : loc_begin(map.begin()),
      loc_end(map.end()),
      loc_iter(end?loc_end:loc_begin)
{
    if(loc_iter!=loc_end) {
        bs_iter=loc_iter->second.begin();
        this->increment_loc();
    }
}


template<class DS, class HBS> inline
bool
HybridSetConstIterator<DS,HBS>::equal(const HybridSetConstIterator<DS,HBS>& other) const
{
    return this->loc_iter==other.loc_iter && (this->loc_iter==this->loc_end || this->bs_iter==other.bs_iter);
}


template<class DS, class HBS> inline
HBS const&
HybridSetConstIterator<DS,HBS>::dereference() const
{
    this->hybrid_set=HBS(loc_iter->first,this->bs_iter->box());
    return this->hybrid_set;
}


template<class DS, class HBS> inline
void
HybridSetConstIterator<DS,HBS>::increment()
{
    ++this->bs_iter;
    this->increment_loc();
}

template<class DS, class HBS> inline
void
HybridSetConstIterator<DS,HBS>::increment_loc()
{
    while(bs_iter==loc_iter->second.end()) {
        ++loc_iter;
        if(loc_iter==loc_end) { return; }
        bs_iter=loc_iter->second.begin();
    }
}


template<class A> void serialize(A& archive, HybridDenotableSet& set, const unsigned int version) {
    archive & static_cast<std::map<DiscreteLocation,DenotableSetType>&>(set); }

template<class BS> inline
void
draw(FigureInterface& figure, const Orbit< BS >& orbit) {
    draw(figure,orbit.reach());
    draw(figure,orbit.initial());
    draw(figure,orbit.final());
}

template<class BS> inline
void
draw(FigureInterface& figure, const HybridBasicSet<BS>& hs) {
    draw(figure,hs.continuous_state_set());
}



template<class DS> inline
void
draw(FigureInterface& figure, const std::map<DiscreteLocation,DS>& hds) {
    for(typename std::map<DiscreteLocation,DS>::const_iterator loc_iter=hds.begin();
        loc_iter!=hds.end(); ++loc_iter) {
        draw(figure,loc_iter->second);
        //figure.draw(loc_iter->second);
        //figure << loc_iter->second;
    }
}

template<class BS> inline FigureInterface& operator<<(FigureInterface& figure, const Orbit< HybridBasicSet<BS> >& horb) {
    draw(figure,horb); return figure;
}

template<class BS> inline FigureInterface& operator<<(FigureInterface& figure, const HybridBasicSet<BS>& hs) {
    draw(figure,hs); return figure;
}

template<class DS> inline FigureInterface& operator<<(FigureInterface& figure, const std::map<DiscreteLocation,DS>& hs) {
    draw(figure,hs); return figure;
}

template<class X> std::map<DiscreteLocation,Vector<X> > max_elementwise(
		const std::map<DiscreteLocation,Vector<X> >& hv1, const std::map<DiscreteLocation,Vector<X> >& hv2) {
	std::map<DiscreteLocation,Vector<X> > result;
	for (typename std::map<DiscreteLocation,Vector<X> >::const_iterator hv1_it = hv1.begin(); hv1_it != hv1.end(); ++hv1_it) {
		typename std::map<DiscreteLocation,Vector<X> >::const_iterator hv2_it = hv2.find(hv1_it->first);
		ARIADNE_ASSERT_MSG(hv2_it != hv2.end(), "The location " << hv1_it->first << " was not found in the second hybrid vector.");

		result.insert(std::pair<DiscreteLocation,Vector<X> >(hv1_it->first,max_elementwise(hv1_it->second,hv2_it->second)));
	}
	return result;
}

template<class X> std::map<DiscreteLocation,Vector<X> > min_elementwise(
		const std::map<DiscreteLocation,Vector<X> >& hv1, const std::map<DiscreteLocation,Vector<X> >& hv2) {
	std::map<DiscreteLocation,Vector<X> > result;
	for (typename std::map<DiscreteLocation,Vector<X> >::const_iterator hv1_it = hv1.begin(); hv1_it != hv1.end(); ++hv1_it) {
		typename std::map<DiscreteLocation,Vector<X> >::const_iterator hv2_it = hv2.find(hv1_it->first);
		ARIADNE_ASSERT_MSG(hv2_it != hv2.end(), "The location " << hv1_it->first << " was not found in the second hybrid vector.");

		result.insert(std::pair<DiscreteLocation,Vector<X> >(hv1_it->first,min_elementwise(hv1_it->second,hv2_it->second)));
	}
	return result;
}

template<class ES>
HybridDenotableSet
outer_approximation(const ListSet<HybridBasicSet<ES> >& hls,
                    const HybridGrid& hgr,
                    const int accuracy)
{
    HybridDenotableSet result(hgr);
    for(typename HybridListSet<ES>::const_iterator
            iter=hls.begin(); iter!=hls.end(); ++iter)
        {
            DiscreteLocation loc=iter->first;
            const ES& es=iter->second;
            if(result.find(loc)==result.locations_end()) {
                result.insert(make_pair(loc,DenotableSetType(hgr[loc])));
            }
            DenotableSetType& gts=result[loc];
            Box bbx = es.bounding_box();
            gts.adjoin_outer_approximation(es.bounding_box(),accuracy);
        }
    return result;
}

template<class ES>
HybridDenotableSet
outer_approximation(const HybridBasicSet<ES>& hs,
                    const HybridGrid& hgr,
                    const int accuracy)
{
    HybridDenotableSet result(hgr);
    DiscreteLocation loc=hs.location();
    const ES& es=hs.continuous_state_set();
    if(result.find(loc)==result.locations_end()) {
        result.insert(make_pair(loc,DenotableSetType(hgr[loc])));
    }
    DenotableSetType& gts=result[loc];
    gts.adjoin_outer_approximation(es.bounding_box(),accuracy);

    return result;
}


//! \brief Checks if \a theSet1 is a subset of \a theSet2.
//! \details If the sets do not have the same grid, an error is raised.
bool subset(const HybridDenotableSet& theSet1, const HybridDenotableSet& theSet2);

//! \brief Whether \a cons_set is disjoint from \a grid_set.
tribool disjoint(const HybridConstraintSet& cons_set, const HybridDenotableSet& grid_set);
//! \brief Whether \a cons_set overlaps with \a grid_set.
tribool overlaps(const HybridConstraintSet& cons_set, const HybridDenotableSet& grid_set);
//! \brief Whether \a cons_set covers \a grid_set.
tribool covers(const HybridConstraintSet& cons_set, const HybridDenotableSet& grid_set);

//! \brief An outer approximation of the intersection of \a den_set with \a cons_set.
HybridDenotableSet outer_intersection(const HybridDenotableSet& den_set, const HybridConstraintSet& cons_set);
//! \brief An inner approximation of the intersection of \a den_set with \a cons_set.
HybridDenotableSet inner_intersection(const HybridDenotableSet& den_set, const HybridConstraintSet& cons_set);
//! \brief An outer approximation of the difference of \a den_set with \a cons_set.
HybridDenotableSet outer_difference(const HybridDenotableSet& den_set, const HybridConstraintSet& cons_set);
//! \brief An inner approximation of the difference of \a den_set with \a cons_set.
HybridDenotableSet inner_difference(const HybridDenotableSet& den_set, const HybridConstraintSet& cons_set);

//! \brief Evaluates the codomain of the function \a func applied on the cells of \a grid_set, each widened by \a eps.
HybridBoxes eps_codomain(const HybridDenotableSet& grid_set, const HybridFloatVector& eps, const HybridVectorFunction& func);

//! \brief Projects \a grid_set using the given \a indices, also flattening the set in respect to the hybrid space.
//! \details The \a grid_set must have the same grid for all locations.
DenotableSetType flatten_and_project_down(const HybridDenotableSet& den_set, const Vector<uint>& indices);

} // namespace Ariadne

namespace boost { namespace serialization {
template<class A> void serialize(A& archive, const Ariadne::HybridDenotableSet& set, const uint version);
template<class A> void serialize(A& archive, const Ariadne::DiscreteLocation& state, const uint version);
}}

#endif // ARIADNE_HYBRID_SET_H
