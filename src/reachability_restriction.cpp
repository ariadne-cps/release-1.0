/***************************************************************************
 *            reachability_restriction.cc
 *
 *  Copyright 2011  Luca Geretti
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
 
#include "reachability_restriction.h"
#include "hybrid_automaton_interface.h"
#include "set_checker.h"
#include "denotable_set.h"

namespace Ariadne {


ReachabilityRestriction::
ReachabilityRestriction(
		const HybridBoxes& domain,
		const HybridGrid& grid,
		int accuracy) :
		_domain(domain),
		_grid(grid),
		_accuracy(accuracy)
{
	for (HybridBoxes::const_iterator loc_it = domain.begin(); loc_it != domain.end(); ++loc_it) {
		HybridGrid::const_iterator grid_it = grid.find(loc_it->first);
		ARIADNE_ASSERT_MSG(grid_it != grid.end(),
				"The location " << loc_it->first.name() << " is not present in the grid space.");
		ARIADNE_ASSERT_MSG(grid_it->second.dimension() == loc_it->second.dimension(),
				"The dimensions in location " << loc_it->first.name() << " between domain and grid do not match.");
	}
}

ReachabilityRestriction::
ReachabilityRestriction(const ReachabilityRestriction& other) :
	_domain(other._domain),
	_grid(other._grid),
	_accuracy(other._accuracy),
	_set(other._set)
{

}


HybridBoxes
ReachabilityRestriction::
bounding_box() const
{
	HybridBoxes result;

	for (HybridBoxes::const_iterator domain_it = _domain.begin(); domain_it != _domain.end(); ++domain_it) {
		const DiscreteLocation& loc = domain_it->first;
		result.insert(make_pair(loc,this->_bounding_box(loc)));
	}
	return result;
}


HybridBoxes
ReachabilityRestriction::
outer_domain_box() const
{
	HybridBoxes result;

	for (HybridBoxes::const_iterator domain_it = _domain.begin(); domain_it != _domain.end(); ++domain_it) {
		const DiscreteLocation& loc = domain_it->first;
		result.insert(make_pair(loc,this->_outer_domain_box(loc)));
	}
	return result;
}


bool
ReachabilityRestriction::
has_location(DiscreteLocation q) const
{
	return (_domain.find(q) != _domain.end());
}


bool
ReachabilityRestriction::
has_discretised(DiscreteLocation q) const
{
    ARIADNE_ASSERT_MSG(this->has_location(q), "The location " << q << " was not found in the HybridRestrictionSet.");
    return _set.find(q) != _set.locations_end();
}


bool
ReachabilityRestriction::
has_empty(DiscreteLocation q) const
{
	ARIADNE_ASSERT_MSG(this->has_location(q), "The location " << q << " was not found in the HybridRestrictionSet.");

	if (_domain.find(q)->second.empty())
		return true;

	if (!this->has_discretised(q))
		return false;

	return _set[q].empty();
}


void
ReachabilityRestriction::
update_with(const HybridDenotableSet& set)
{
	ARIADNE_ASSERT_MSG(_grid == set.grid(), "You cannot update a restriction with a set having a mismatched grid.");

	for (HybridDenotableSet::locations_const_iterator set_it = set.locations_begin(); set_it != set.locations_end(); ++set_it) {
		const DiscreteLocation& loc = set_it->first;

		if (this->has_discretised(loc))
			_set[loc] = set_it->second;
		else
			_set.insert(make_pair(loc,set_it->second));
	}
}


void
ReachabilityRestriction::
refine_at(int accuracy)
{
	_accuracy = accuracy;

	for (HybridDenotableSet::locations_iterator set_it = _set.locations_begin(); set_it != _set.locations_end(); ++set_it) {
		const DiscreteLocation& loc = set_it->first;
		const DenotableSetType& local_set = set_it->second;
		const Box& domain_box = _domain.find(loc)->second;
		const Grid& local_grid = _grid.find(loc)->second;

		if (!local_set.bounding_box().inside(domain_box)) {

			DenotableSetType refined_local_restriction(local_grid);
			refined_local_restriction.adjoin_outer_approximation(domain_box,accuracy);
			_set[loc] = refined_local_restriction;
		}
	}
}


void
ReachabilityRestriction::
apply_to(HybridDenotableSet& set) const
{
	ARIADNE_ASSERT_MSG(_grid == set.grid(), "You cannot apply a restriction to a set with a mismatching grid.");

	for (HybridDenotableSet::locations_iterator result_it = set.locations_begin(); result_it != set.locations_end(); ++result_it) {
		const DiscreteLocation& loc = result_it->first;

		DenotableSetType& local_set = result_it->second;

		if (!local_set.empty()) {

			// If the domain covers the result, no restriction is required. Otherwise we discretise, if necessary, and restrict.
			if (!this->has_discretised(loc) && this->_bounding_box(loc).covers(local_set.bounding_box())) {

			} else {
				if (!this->has_discretised(loc))
					_insert_domain_discretisation(loc);

				local_set.restrict(_set[loc]);
			}
		}
	}
}


std::list<LocalisedEnclosureType>
ReachabilityRestriction::
filter(const std::list<LocalisedEnclosureType>& enclosures) const
{
	std::list<LocalisedEnclosureType> result;

	for (std::list<LocalisedEnclosureType>::const_iterator encl_it = enclosures.begin(); encl_it != enclosures.end(); ++encl_it) {
		const DiscreteLocation& loc = encl_it->first;
		ARIADNE_ASSERT_MSG(this->has_location(loc),
				"The location " << loc.name() << " is not present in the ReachabilityRestriction: cannot apply the restriction.");

		Box encl_bb = encl_it->second.bounding_box();

		// If the domain covers the local enclosure, no restriction is required. Otherwise we discretise, if necessary, and restrict.
		if (!this->has_discretised(loc) && this->_bounding_box(loc).covers(encl_bb)) {
			result.push_back(LocalisedEnclosureType(loc,encl_it->second));
		} else {
			if (!this->has_discretised(loc))
				_insert_domain_discretisation(loc);
/*
			cout << "Set bb: " << _set[loc].bounding_box() << ", Encl bb: " << encl_bb;

			cout << "Set cells bbs: ";

			for (DenotableSetType::const_iterator cell_it = _set[loc].begin(); cell_it != _set[loc].end(); ++cell_it) {
				cout << cell_it->box() << " ";
				if (possibly(!cell_it->box().disjoint(encl_bb)))
					cout << "Not disjoint!";
			}
			cout << endl;
*/
			if (possibly(!_set[loc].disjoint(encl_bb))) {
	//			cout << " ADDED\n";
				result.push_back(LocalisedEnclosureType(loc,encl_it->second));
			} /*else
				cout << " DISCARDED\n";*/
		}
	}

    return result;
}


bool
ReachabilityRestriction::
restricts(const HybridDenotableSet& set) const
{
	// This is a sanity assertion: in principle you could ignore differing locations, but it would silently
	// let pass all improper uses of this method
	ARIADNE_ASSERT_MSG(_grid == set.grid(),
			"You are checking if a set is restricted by a ReachabilityRestriction with a different grid.")

    for (HybridDenotableSet::locations_const_iterator set_it = set.locations_begin(); set_it != set.locations_end(); ++set_it) {
		const DiscreteLocation& loc = set_it->first;
		const DenotableSetType& local_set = set_it->second;

		// If the domain covers the result, no restriction is required
		if (!this->has_discretised(loc) && this->_bounding_box(loc).covers(local_set.bounding_box())) {

		} else {
			if (!this->has_discretised(loc))
				_insert_domain_discretisation(loc);

			DenotableSetType local_set_copy = local_set;
			local_set_copy.restrict(_set[loc]);

			// Equivalent to inequality checking, if done between one set and its subset
			if (!subset(local_set,local_set_copy))
				return true;
		}
    }
	return false;
}


tribool
ReachabilityRestriction::
disjoint(const LocalisedBox& hbx) const
{
	tribool result;

	const DiscreteLocation& loc = hbx.first;
	const Box& bx = hbx.second;

    if (!this->has_location(loc))
    	return true;
	else {
		if (this->_bounding_box(loc).disjoint(bx))
			return true;
		else {
			if (!this->has_discretised(loc))
				_insert_domain_discretisation(loc);

			result = _set.find(loc)->second.disjoint(bx);
		}
	}

    return result;
}


tribool
ReachabilityRestriction::
overlaps(const LocalisedBox& hbx) const
{
	tribool result;

	const DiscreteLocation& loc = hbx.first;
	const Box& bx = hbx.second;

    if (!this->has_location(loc))
    	return false;
	else {
		if (!this->_bounding_box(loc).overlaps(bx))
			return false;
		else {
			if (!this->has_discretised(loc))
				_insert_domain_discretisation(loc);

			result = _set.find(loc)->second.overlaps(bx);
		}
	}

    return result;
}


tribool
ReachabilityRestriction::
superset(const LocalisedBox& hbx) const
{
	tribool result;

	const DiscreteLocation& loc = hbx.first;
	const Box& bx = hbx.second;

    if (!this->has_location(loc))
    	return false;
	else {
		if (!this->_bounding_box(loc).superset(bx))
			return false;
		else {
			if (!this->has_discretised(loc))
				_insert_domain_discretisation(loc);

			result = _set.find(loc)->second.superset(bx);
		}
	}

    return result;
}


HybridDenotableSet
ReachabilityRestriction::
forward_jump_set(
		const HybridDenotableSet& set,
		const HybridAutomatonInterface& sys) const
{
	ARIADNE_ASSERT_MSG(set.grid() == _grid, "To create a forward jump set, the set grid must match the restriction grid.");

	HybridDenotableSet result(_grid);

	for (HybridDenotableSet::locations_const_iterator src_loc_it = set.locations_begin();
			src_loc_it != set.locations_end(); ++src_loc_it) {

		const DiscreteLocation& src_location = src_loc_it->first;

		Set<DiscreteEvent> events = sys.events(src_location);
		for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {
			const DiscreteEvent event = *event_it;
			EventKind kind = sys.event_kind(src_location,event);
			if (kind == URGENT || kind == PERMISSIVE) {
				const DiscreteLocation trg_location = sys.target(src_location,event);
				const Grid& trg_grid = _grid.find(trg_location)->second;

				DenotableSetType event_feasible_src_set = src_loc_it->second;
				ForwardDiscreteJumpSetChecker checker(src_location,sys,event);
				event_feasible_src_set.outer_restrict(checker,_accuracy);

				DenotableSetType event_feasible_trg_set = checker.get_reset(
						event_feasible_src_set,event,trg_grid,_accuracy);

				result[trg_location].adjoin(event_feasible_trg_set);
			}
		}
	}

	// It can still be the case that the discretisation in the target location introduces cells outside the restriction
	this->apply_to(result);

	return result;
}

HybridDenotableSet
ReachabilityRestriction::
backward_jump_set(
		const HybridDenotableSet& set,
		const HybridAutomatonInterface& sys) const
{
	ARIADNE_ASSERT_MSG(set.grid() == _grid, "To create a backward jump set, the target set grid must match the restriction grid.");

	HybridDenotableSet result(_grid);

	HybridSpace space = _grid.state_space();
	for (HybridSpace::const_iterator loc_it = space.begin(); loc_it != space.end(); ++loc_it) {

		const DiscreteLocation& src_location = loc_it->first;

		//! TODO: should be improved, since we need all locations discretised
		if (!has_discretised(src_location))
			_insert_domain_discretisation(src_location);

		DenotableSetType src_set = _set[src_location];

		BackwardDiscreteJumpSetChecker checker(src_location,sys,set);
		src_set.outer_restrict(checker,_accuracy);
		result[src_location] = src_set;
	}

	return result;
}


HybridDenotableSet
ReachabilityRestriction::
outer_intersection_with(const HybridConstraintSet& constraint_set) const
{
	HybridDenotableSet result;

	for (HybridConstraintSet::const_iterator cons_it = constraint_set.begin(); cons_it != constraint_set.end(); ++cons_it) {
		const DiscreteLocation& loc = cons_it->first;
		if (this->has_location(loc)) {
			if (this->has_empty(loc)) {
				result.insert(make_pair(loc,DenotableSetType(_grid.find(loc)->second)));
			} else {
				if (!this->has_discretised(loc))
					_insert_domain_discretisation(loc);
				result.insert(make_pair(loc,Ariadne::outer_intersection(_set[loc],cons_it->second)));
			}
		}
	}

	return result;
}


HybridDenotableSet
ReachabilityRestriction::
inner_intersection_with(const HybridConstraintSet& constraint_set) const
{
	HybridDenotableSet result;

	for (HybridConstraintSet::const_iterator cons_it = constraint_set.begin(); cons_it != constraint_set.end(); ++cons_it) {
		const DiscreteLocation& loc = cons_it->first;
		if (this->has_location(loc)) {
			if (this->has_empty(loc)) {
				result.insert(make_pair(loc,DenotableSetType(_grid.find(loc)->second)));
			} else {
				if (!this->has_discretised(loc))
					_insert_domain_discretisation(loc);

				result.insert(make_pair(loc,Ariadne::inner_intersection(_set[loc],cons_it->second)));
			}
		}
	}

	return result;
}


HybridDenotableSet
ReachabilityRestriction::
outer_difference_from(const HybridConstraintSet& constraint_set) const
{
	HybridDenotableSet result;

	HybridSpace space = _grid.state_space();

	for (HybridGrid::const_iterator grid_it = _grid.begin(); grid_it != _grid.end(); ++grid_it) {
		const DiscreteLocation& loc = grid_it->first;

		if (this->has_empty(loc)) {
			result.insert(make_pair(loc,DenotableSetType(grid_it->second)));
		} else {
			if (!this->has_discretised(loc))
				_insert_domain_discretisation(loc);

			DenotableSetType local_result = _set[loc];
			HybridConstraintSet::const_iterator cons_it = constraint_set.find(loc);
		    if (cons_it != constraint_set.end()) {
		    	DenotableSetType inner_intersection = Ariadne::inner_intersection(_set[loc],cons_it->second);
		    	local_result.remove(inner_intersection);
		    }
			result.insert(make_pair(loc,local_result));
		}
	}

	return result;
}


HybridDenotableSet
ReachabilityRestriction::
inner_difference_from(const HybridConstraintSet& constraint_set) const
{
	HybridDenotableSet result;

	HybridSpace space = _grid.state_space();

	for (HybridGrid::const_iterator grid_it = _grid.begin(); grid_it != _grid.end(); ++grid_it) {
		const DiscreteLocation& loc = grid_it->first;

		if (this->has_empty(loc)) {
			result.insert(make_pair(loc,DenotableSetType(grid_it->second)));
		} else {
			if (!this->has_discretised(loc))
				_insert_domain_discretisation(loc);

			DenotableSetType local_result = _set[loc];
			HybridConstraintSet::const_iterator cons_it = constraint_set.find(loc);
		    if (cons_it != constraint_set.end()) {
		    	DenotableSetType outer_intersection = Ariadne::outer_intersection(_set[loc],cons_it->second);
		    	local_result.remove(outer_intersection);
		    }
			result.insert(make_pair(loc,local_result));
		}
	}

	return result;
}


HybridDenotableSet
ReachabilityRestriction::
outer_intersection_with(const HybridBoundedConstraintSet& constraint_set) const
{
    ReachabilityRestriction domain_restriction(constraint_set.domain(),_grid,_accuracy);

    HybridConstraintSet unbounded_constraint_set(constraint_set.functions(),constraint_set.codomain());
    HybridDenotableSet result = domain_restriction.outer_intersection_with(unbounded_constraint_set);

    this->apply_to(result);

	return result;
}


HybridDenotableSet
ReachabilityRestriction::
inner_intersection_with(const HybridBoundedConstraintSet& constraint_set) const
{
    ReachabilityRestriction domain_restriction(constraint_set.domain(),_grid,_accuracy);

    HybridConstraintSet unbounded_constraint_set(constraint_set.functions(),constraint_set.codomain());
    HybridDenotableSet result = domain_restriction.inner_intersection_with(unbounded_constraint_set);
    this->apply_to(result);

	return result;
}


HybridDenotableSet
ReachabilityRestriction::
outer_difference_from(const HybridBoundedConstraintSet& constraint_set) const
{
    ReachabilityRestriction domain_restriction(constraint_set.domain(),_grid,_accuracy);

    HybridConstraintSet unbounded_constraint_set(constraint_set.functions(),constraint_set.codomain());
    HybridDenotableSet result = domain_restriction.outer_difference_from(unbounded_constraint_set);
    this->apply_to(result);

	return result;
}


HybridDenotableSet
ReachabilityRestriction::
inner_difference_from(const HybridBoundedConstraintSet& constraint_set) const
{
    ReachabilityRestriction domain_restriction(constraint_set.domain(),_grid,_accuracy);

    HybridConstraintSet unbounded_constraint_set(constraint_set.functions(),constraint_set.codomain());
    HybridDenotableSet result = domain_restriction.inner_difference_from(unbounded_constraint_set);
    this->apply_to(result);

	return result;
}


std::ostream&
ReachabilityRestriction::
write(std::ostream& os) const
{
    return os << "ReachabilityRestriction( domain=" << _domain <<
    		", grid=" << _grid << ", accuracy=" << _accuracy << ")";
}


void
ReachabilityRestriction::
_insert_domain_discretisation(DiscreteLocation q) const
{
	DenotableSetType discretisation(_grid.find(q)->second);
	discretisation.adjoin_outer_approximation(_domain.find(q)->second,_accuracy);
	const_cast<HybridDenotableSet&>(_set).insert(make_pair(q,discretisation));
}


Box
ReachabilityRestriction::
_bounding_box(DiscreteLocation q) const
{
    ARIADNE_ASSERT_MSG(this->has_location(q), "The location " << q << " was not found in the HybridRestrictionSet.");

    return (this->has_discretised(q) ? _set[q].bounding_box() : _outer_domain_box(q));
}

Box
ReachabilityRestriction::
_outer_domain_box(DiscreteLocation q) const
{
    const Float accuracy_divider = (1<<_accuracy);

	Box result = _domain.find(q)->second;

	const Vector<Float>& grid_lengths = _grid.find(q)->second.lengths();
	// Add the minimum cell size on each extreme of the domain box
	for (uint i=0; i<result.size();++i) {
		Float inaccuracy = grid_lengths[i]/accuracy_divider;
		result[i] = Interval(result[i].lower()-inaccuracy,result[i].upper()+inaccuracy);
	}
	return result;
}


} // namespace Ariadne
