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
#include "taylor_calculus.h"

namespace Ariadne {


ReachabilityRestriction::
ReachabilityRestriction(
		const HybridBoxes& domain,
		const HybridGrid& grid,
		int accuracy) :
		_domain(domain),
		_grid(grid),
		_accuracy(accuracy),
		_calculus(new TaylorCalculus())
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
	_set(other._set),
	_calculus(other._calculus)
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


bool
ReachabilityRestriction::
has_discretised(DiscreteLocation q) const
{
    ARIADNE_ASSERT_MSG(this->has_location(q), "The location " << q << " was not found in the HybridRestrictionSet.");
    return _set.find(q) != _set.locations_end();
}


bool
ReachabilityRestriction::
has_location(DiscreteLocation q) const
{
	return (_domain.find(q) != _domain.end());
}


void
ReachabilityRestriction::
update(const HybridDenotableSet& set)
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
refine(int accuracy)
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


HybridDenotableSet
ReachabilityRestriction::
apply_to(const HybridDenotableSet& set) const
{
	ARIADNE_ASSERT_MSG(_grid == set.grid(), "You cannot apply a restriction to a set with a mismatching grid.");

	HybridDenotableSet result = set;

	for (HybridDenotableSet::locations_iterator result_it = result.locations_begin(); result_it != result.locations_end(); ++result_it) {
		const DiscreteLocation& loc = result_it->first;

		DenotableSetType& local_result = result_it->second;

		// If the domain covers the result, no restriction is required. Otherwise we discretise, if necessary, and restrict.
		if (!this->has_discretised(loc) && this->_bounding_box(loc).covers(local_result.bounding_box())) {

		} else {
			if (!this->has_discretised(loc))
				_insert_domain_discretisation(loc);

			local_result.restrict(_set[loc]);
		}
	}
	return result;
}


std::list<EnclosureType>
ReachabilityRestriction::
apply_to(const std::list<EnclosureType>& enclosures) const
{
	std::list<EnclosureType> result;

	for (std::list<EnclosureType>::const_iterator encl_it = enclosures.begin(); encl_it != enclosures.end(); ++encl_it) {
		const DiscreteLocation& loc = encl_it->first;
		ARIADNE_ASSERT_MSG(this->has_location(loc),
				"The location " << loc.name() << " is not present in the ReachabilityRestriction: cannot apply the restriction.");

		Box encl_bb = encl_it->second.bounding_box();

		// If the domain covers the local enclosure, no restriction is required. Otherwise we discretise, if necessary, and restrict.
		if (!this->has_discretised(loc) && this->_bounding_box(loc).covers(encl_bb)) {
			result.push_back(EnclosureType(loc,encl_it->second));
		} else {
			if (!this->has_discretised(loc))
				_insert_domain_discretisation(loc);

			if (possibly(_set[loc].overlaps(encl_bb)))
				result.push_back(EnclosureType(loc,encl_it->second));
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


HybridDenotableSet
ReachabilityRestriction::
forward_jump_set(
		const HybridDenotableSet& set,
		const HybridAutomatonInterface& sys) const
{
	/*
	 * \forall DiscreteLocation l1 in set locations {
	 * 		\forall cells c in set[l1] not outside invariants {
	 * 			\forall possibly active transition t to location l2 whose reset is not outside invariants {
	 * 				adjoin reset(c) to result
	 * 			}
	 * 		}
	 * 	}
	 *  apply the restriction on the result
	 */

	ARIADNE_ASSERT_MSG(set.grid() == _grid, "To create a forward jump set, the target set grid must match the restriction grid.");

	HybridDenotableSet result(_grid);

	for (HybridDenotableSet::locations_const_iterator src_loc_it = set.locations_begin();
			src_loc_it != set.locations_end(); ++src_loc_it) {

		const DiscreteLocation& src_location = src_loc_it->first;
		const DenotableSetType& src_cells = src_loc_it->second;

		for (DenotableSetType::const_iterator src_cell_it = src_cells.begin(); src_cell_it != src_cells.end(); ++src_cell_it) {

			Box src_cell_bx = src_cell_it->box();
			const ContinuousEnclosureType src_enclosure(src_cell_bx);

			if (!_is_outside_invariants(src_location,src_cell_bx,sys))
				_adjoin_forward_jump_sets(src_location,src_enclosure,sys,result);
		}
	}

	return this->apply_to(result);
}

HybridDenotableSet
ReachabilityRestriction::
backward_jump_set(
		const HybridDenotableSet& set,
		const HybridAutomatonInterface& sys) const
{
	/*
	 * \forall DiscreteLocation l1 in non-empty restriction locations {
	 * 		if \exists transition t1 to location l2 in targetCells and r_t1(bounding[l1]) intersects targetCells[l2] {
	 * 			if not discretised[l2]
	 * 				discretise in l2
	 * 			\forall cells c in restriction[l2] not outside invariants {
	 * 				if \exists possibly active transition t whose reset overlaps targetCells[l2] {
	 * 					adjoin c to result
	 * 					break cells scan
	 * 				}
	 * 			}
	 * 		}
	 * 	}
	 */

	ARIADNE_ASSERT_MSG(set.grid() == _grid, "To create a backward jump set, the target set grid must match the restriction grid.");

	HybridDenotableSet result(_grid);

	HybridSpace space = _grid.state_space();
	for (HybridSpace::const_iterator src_space_it = space.begin(); src_space_it != space.end(); ++src_space_it) {

		const DiscreteLocation& src_location = src_space_it->first;
		if (!(this->has_discretised(src_location) && _set[src_location].empty())) {

			if (_has_feasible_transition(src_location,set,sys)) {

				if (!this->has_discretised(src_location))
					_insert_domain_discretisation(src_location);

				DenotableSetType src_cells = _set[src_location];
				for (DenotableSetType::const_iterator src_cell_it = src_cells.begin(); src_cell_it != src_cells.end(); ++src_cell_it) {

					Box src_cell_bx = src_cell_it->box();
					const ContinuousEnclosureType src_enclosure(src_cell_bx);

					if (!_is_outside_invariants(src_location,src_cell_bx,sys))
						_adjoin_src_of_forward_jump_sets(src_location,src_enclosure,sys,set,result);
				}
			}
		}
	}

	return result;
}


HybridDenotableSet
ReachabilityRestriction::
possibly_feasible_projection(const HybridConstraintSet& constraint) const
{
	HybridDenotableSet result;

	for (HybridConstraintSet::const_iterator cons_it = constraint.begin(); cons_it != constraint.end(); ++cons_it) {
		const DiscreteLocation& loc = cons_it->first;
		ARIADNE_ASSERT_MSG(this->has_location(loc), "The constraint has location " << loc.name() << ", but the restriction does not.");
		if (!this->has_discretised(loc))
			_insert_domain_discretisation(loc);

		result.insert(make_pair(loc,possibly_overlapping_subset(_set[loc],cons_it->second)));
	}

	return result;
}


HybridDenotableSet
ReachabilityRestriction::
definitely_feasible_projection(const HybridConstraintSet& constraint) const
{
	HybridDenotableSet result;

	for (HybridConstraintSet::const_iterator cons_it = constraint.begin(); cons_it != constraint.end(); ++cons_it) {
		const DiscreteLocation& loc = cons_it->first;
		ARIADNE_ASSERT_MSG(this->has_location(loc), "The constraint has location " << loc.name() << ", but the restriction does not.");
		if (!this->has_discretised(loc))
			_insert_domain_discretisation(loc);

		result.insert(make_pair(loc,definitely_covered_subset(_set[loc],cons_it->second)));
	}

	return result;
}


HybridDenotableSet
ReachabilityRestriction::
possibly_infeasible_projection(const HybridConstraintSet& constraint) const
{
	HybridDenotableSet result;

	for (HybridConstraintSet::const_iterator cons_it = constraint.begin(); cons_it != constraint.end(); ++cons_it) {
		const DiscreteLocation& loc = cons_it->first;
		ARIADNE_ASSERT_MSG(this->has_location(loc), "The constraint has location " << loc.name() << ", but the restriction does not.");
		if (!this->has_discretised(loc))
			_insert_domain_discretisation(loc);

		DenotableSetType definitely_covered = definitely_covered_subset(_set[loc],cons_it->second);
		DenotableSetType possibly_infeasible = _set[loc];
		possibly_infeasible.remove(definitely_covered);

		result.insert(make_pair(loc,possibly_infeasible));
	}

	return result;
}


HybridDenotableSet
ReachabilityRestriction::
definitely_infeasible_projection(const HybridConstraintSet& constraint) const
{
	HybridDenotableSet result;

	for (HybridConstraintSet::const_iterator cons_it = constraint.begin(); cons_it != constraint.end(); ++cons_it) {
		const DiscreteLocation& loc = cons_it->first;
		ARIADNE_ASSERT_MSG(this->has_location(loc), "The constraint has location " << loc.name() << ", but the restriction does not.");
		if (!this->has_discretised(loc))
			_insert_domain_discretisation(loc);

		DenotableSetType possibly_overlapping = possibly_overlapping_subset(_set[loc],cons_it->second);
		DenotableSetType definitely_infeasible = _set[loc];
		definitely_infeasible.remove(possibly_overlapping);

		result.insert(make_pair(loc,definitely_infeasible));
	}

	return result;
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

    const Float accuracy_divider = (1<<_accuracy);

    if (this->has_discretised(q)) {
    	return _set[q].bounding_box();
    } else {
    	const Vector<Float>& grid_lengths = _grid.find(q)->second.lengths();
    	Box result = _domain.find(q)->second;
    	// Add the minimum cell size on each extreme of the domain box
    	for (uint i=0; i<result.size();++i) {
    		Float inaccuracy = grid_lengths[i]/accuracy_divider;
    		result[i] = Interval(result[i].lower()-inaccuracy,result[i].upper()+inaccuracy);
    	}
    	return result;
    }
}


bool
ReachabilityRestriction::
_has_feasible_transition(
		const DiscreteLocation& src_location,
		const HybridDenotableSet& trg_set,
		const HybridAutomatonInterface& sys) const
{
	const ContinuousEnclosureType src_enclosure(this->bounding_box().at(src_location));
	const Set<DiscreteEvent> events = sys.events(src_location);

	for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {

		const DiscreteEvent event = *event_it;
		EventKind kind = sys.event_kind(src_location,event);

		if (kind == URGENT || kind == PERMISSIVE) {

			DiscreteLocation trg_location = sys.target(src_location,event);

			if (!trg_set[trg_location].empty()) {

				const ContinuousEnclosureType target_encl = _calculus->reset_step(
						sys.reset_function(src_location,event),src_enclosure);

				const HybridBox target_hbounding(trg_location,target_encl.bounding_box());

				if (possibly(trg_set.overlaps(target_hbounding)))
					return true;
			}
		}
	}
	return false;
}

void
ReachabilityRestriction::
_adjoin_forward_jump_sets(
		const DiscreteLocation& src_location,
		const ContinuousEnclosureType& src_enclosure,
		const HybridAutomatonInterface& sys,
		HybridDenotableSet& result_set) const
{
	Set<DiscreteEvent> events = sys.events(src_location);
	for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {
		const DiscreteEvent event = *event_it;
		EventKind kind = sys.event_kind(src_location,event);
		if (kind == URGENT || kind == PERMISSIVE) {

			RealScalarFunction guard = sys.guard_function(src_location,event);
			RealVectorFunction dynamic = sys.dynamic_function(src_location);

			if (_is_transition_feasible(guard,kind,dynamic,src_enclosure)) {
				const DiscreteLocation& target_loc = sys.target(src_location,event);
				const ContinuousEnclosureType target_encl = _calculus->reset_step(
						sys.reset_function(src_location,event),src_enclosure);

				if (!_is_outside_invariants(target_loc,target_encl.bounding_box(),sys))
					result_set[target_loc].adjoin_outer_approximation(target_encl,_accuracy);
			}
		}
	}
}


void
ReachabilityRestriction::
_adjoin_src_of_forward_jump_sets(
		const DiscreteLocation& src_location,
		const ContinuousEnclosureType& src_enclosure,
		const HybridAutomatonInterface& sys,
		const HybridDenotableSet& trg_set,
		HybridDenotableSet& set_to_adjoin) const
{
	const Set<DiscreteEvent> events = sys.events(src_location);

	// In order to add a source hybrid box, just one possibly overlapping target enclosure suffices
	for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {
		const DiscreteEvent event = *event_it;
		EventKind kind = sys.event_kind(src_location,event);
		if (kind == URGENT || kind == PERMISSIVE) {

			RealScalarFunction guard = sys.guard_function(src_location,event);
			RealVectorFunction dynamic = sys.dynamic_function(src_location);

			if (_is_transition_feasible(guard,kind,dynamic,src_enclosure)) {
				const DiscreteLocation& target_loc = sys.target(src_location,event);
				const ContinuousEnclosureType image_encl = _calculus->reset_step(
						sys.reset_function(src_location,event),src_enclosure);
				const HybridBox image_hbx(target_loc,image_encl.bounding_box());

				if (possibly(trg_set.overlaps(image_hbx))) {
					set_to_adjoin[src_location].adjoin_outer_approximation(src_enclosure,_accuracy);
					break;
				}
			}
		}
	}

}


bool
ReachabilityRestriction::
_is_transition_feasible(
		const ScalarFunction& activation,
		EventKind event_kind,
		const VectorFunction& dynamic,
		const ContinuousEnclosureType& source) const
{
	bool result = false;

	const bool is_urgent = (event_kind == URGENT);

	tribool is_guard_active = _calculus->active(VectorFunction(1,activation),source);

	ARIADNE_LOG(6,"Guard activity: " << is_guard_active);

	/*
	 * a) If the guard is definitely active and the transition is urgent, then we are definitely outside the related invariant
	 * b) If the transition is not urgent, it suffices to have a possibly active guard: we then must perform the transition
	 * c) If the transition is urgent and the guard is only possibly active, we check the crossing:
	 *    i) If it is negative, then no transition is possible
	 *	 ii) If it is possibly positive, then we must take the transition
	 */

	if (definitely(is_guard_active) && is_urgent) {
		result = false;
	} else if (possibly(is_guard_active) && !is_urgent) {
		result = true;
	} else if (possibly(is_guard_active) && is_urgent) {
		tribool positive_crossing = is_positively_crossing(source.bounding_box(),dynamic,activation);
		result = possibly(positive_crossing);
	} else {
		result = false;
	}

	return result;
}


bool
ReachabilityRestriction::
_is_outside_invariants(
		const DiscreteLocation& location,
		const Box& bx,
		const HybridAutomatonInterface& sys) const
{
	Set<DiscreteEvent> events = sys.events(location);
	for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {
		const DiscreteEvent event = *event_it;
		EventKind kind = sys.event_kind(location,event);
		if (kind == INVARIANT) {
			const ScalarFunction& activation = sys.invariant_function(location,event);
			tribool is_active = _calculus->active(VectorFunction(1,activation),bx);
			if (definitely(is_active)) {
				return true;
			}
		}
	}

	return false;
}


tribool
is_positively_crossing(
		const Box& set_bounds,
		const RealVectorFunction& dynamic,
		const RealScalarFunction& activation)
{
    RealScalarFunction derivative=lie_derivative(activation,dynamic);
    Interval derivative_range = derivative.evaluate(set_bounds);

    if (derivative_range.lower() > 0)
    	return true;
    else if (derivative_range.upper() < 0)
    	return false;
    else return indeterminate;
}


} // namespace Ariadne
