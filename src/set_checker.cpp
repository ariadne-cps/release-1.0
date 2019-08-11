/***************************************************************************
 *            set_checker.cc
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
 
#include "set_checker.h"
#include "reachability_restriction.h"
#include "hybrid_automaton_interface.h"
#include "taylor_calculus.h"

namespace Ariadne {


ConstraintSetChecker::
ConstraintSetChecker(const ConstraintSet& constraint) :
	_constraint(constraint.clone())
{

}

tribool
ConstraintSetChecker::
check(const Box& bx) const
{
	if (definitely(_constraint->covers(bx)))
		return true;
	if (definitely(_constraint->disjoint(bx)))
		return false;
	return indeterminate;
}

DiscreteJumpSetCheckerBase::
DiscreteJumpSetCheckerBase(
		const DiscreteLocation& loc,
		const HybridAutomatonInterface& sys) :
		_location(loc),
		_sys(sys.clone()),
		_calculus(new TaylorCalculus())
{
	HybridSpace sys_space = sys.state_space();

	ARIADNE_ASSERT_MSG(sys_space.find(loc) != sys_space.end(),
			"The set location is not present in the system space.");
}


tribool
DiscreteJumpSetCheckerBase::
_is_transition_taken(
		const ScalarFunction& activation,
		EventKind event_kind,
		const VectorFunction& dynamic,
		const EnclosureType& src_encl) const
{
	tribool result;

	Box src_encl_bx = src_encl.box();

	const bool is_urgent = (event_kind == URGENT);

	tribool is_guard_active = _calculus->active(VectorFunction(1,activation),src_encl);

	/*
	 * a) If the guard is definitely active and the transition is urgent, then we are definitely outside the related invariant
	 * b) If the transition is not urgent, then taking the transition reflect the activity (this would disallow the transition for inner approximations)
	 * c) If the transition is urgent and the guard is only possibly active, we refine by checking the crossing:
	 *    - if it is definitely negative, then no transition is possible
	 *    - otherwise we leave the possibility
	 * d) If the activity is false, the result is obviously false
	 */

	if (definitely(is_guard_active) && is_urgent) {
		result = false;
	} else if (possibly(is_guard_active) && !is_urgent) {
		result = is_guard_active;
	} else if (possibly(is_guard_active) && is_urgent) {
		tribool positive_crossing = _is_positively_crossing(src_encl.bounding_box(),dynamic,activation);
		if (definitely(!positive_crossing))
			result = false;
		else
			result = indeterminate;
	} else {
		result = false;
	}

	return result;
}


tribool
DiscreteJumpSetCheckerBase::
_is_outside_any_invariant(
		const DiscreteLocation& location,
		const Box& bx) const
{
	tribool result = false;

	Set<DiscreteEvent> events = _sys->events(location);
	for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {
		const DiscreteEvent event = *event_it;
		EventKind kind = _sys->event_kind(location,event);
		if (kind == INVARIANT || kind == URGENT) {
			const ScalarFunction& activation = (kind == INVARIANT ?
					_sys->invariant_function(location,event) : _sys->guard_function(location,event));
			tribool is_active = _calculus->active(VectorFunction(1,activation),bx);
			if (definitely(is_active)) {
				return true;
			} else if (possibly(is_active)) {
				result = indeterminate;
			}
		}
	}

	return result;
}


tribool
DiscreteJumpSetCheckerBase::
_is_positively_crossing(
		const Box& set_bounds,
		const RealVectorFunction& dynamic,
		const RealScalarFunction& activation) const
{
    RealScalarFunction derivative= lie_derivative(activation,dynamic);
    Interval derivative_range = derivative.evaluate(set_bounds);

    if (derivative_range.lower() > 0)
    	return true;
    else if (derivative_range.upper() < 0)
    	return false;
    else
    	return indeterminate;
}


ForwardDiscreteJumpSetChecker::
ForwardDiscreteJumpSetChecker(
		const DiscreteLocation& loc,
		const HybridAutomatonInterface& sys,
		DiscreteEvent event) :
		DiscreteJumpSetCheckerBase(loc,sys),
		_event(event)
{
	Set<DiscreteEvent> events = sys.events(loc);
	ARIADNE_ASSERT_MSG(events.find(event) != events.end(),
			"The provided event " << event.name() << " is not present in the system at location " << loc.name());
	EventKind kind = sys.event_kind(loc,event);
	ARIADNE_ASSERT_MSG(kind == URGENT || kind == PERMISSIVE, "The event must be associated with a transition.");
}


tribool
ForwardDiscreteJumpSetChecker::
check(const Box& bx) const
{
	tribool is_outside_any_src_invariant = _is_outside_any_invariant(_location,bx);

	if (definitely(is_outside_any_src_invariant))
		return false;

	tribool result = false;

	EventKind kind = _sys->event_kind(_location,_event);

	RealScalarFunction guard = _sys->guard_function(_location,_event);
	RealVectorFunction dynamic = _sys->dynamic_function(_location);

	EnclosureType src_enclosure(bx);

	tribool is_transition_taken = _is_transition_taken(guard,kind,dynamic,src_enclosure);

	if (possibly(is_transition_taken)) {
		const DiscreteLocation& trg_location = _sys->target(_location,_event);
		const EnclosureType trg_enclosure = _calculus->reset_step(
				_sys->reset_function(_location,_event),src_enclosure);

		tribool is_outside_any_trg_invariant = _is_outside_any_invariant(trg_location,trg_enclosure.bounding_box());

		result = (!is_outside_any_src_invariant && is_transition_taken && !is_outside_any_trg_invariant);

	}

	return result;
}


DenotableSetType
ForwardDiscreteJumpSetChecker::
get_reset(
		const DenotableSetType& src_set,
		DiscreteEvent event,
		const Grid& trg_grid,
		int accuracy) const
{
	DiscreteLocation trg_location = _sys->target(_location,event);

	DenotableSetType result(trg_grid);

	for (DenotableSetType::const_iterator cell_it = src_set.begin(); cell_it != src_set.end(); ++cell_it) {
		EnclosureType src_enclosure(cell_it->box());
		EnclosureType trg_enclosure = _calculus->reset_step(
				_sys->reset_function(_location,event),src_enclosure);
		result.adjoin_outer_approximation(trg_enclosure,accuracy);
	}

	return result;
}


BackwardDiscreteJumpSetChecker::
BackwardDiscreteJumpSetChecker(
		const DiscreteLocation& loc,
		const HybridAutomatonInterface& sys,
		const HybridDenotableSet& reset_restriction) :
		DiscreteJumpSetCheckerBase(loc,sys),
		_starting_set(reset_restriction.clone())
{
	HybridSpace sys_space = sys.state_space();

	ARIADNE_ASSERT_MSG(sys_space == reset_restriction.space(),
			"The system and reset restriction have different state spaces.");
}


tribool
BackwardDiscreteJumpSetChecker::
check(const Box& bx) const
{
	// Remember that source and target are absolute and refer to a forward transition, so they appear
	// switched in backward evolution

	tribool is_outside_any_src_invariant = _is_outside_any_invariant(_location,bx);

	if (definitely(is_outside_any_src_invariant))
		return false;

	if (!_has_feasible_transitions(bx))
		return false;

	const Set<DiscreteEvent> events = _sys->events(_location);

	EnclosureType src_encl(bx);

	tribool result = false;

	for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {
		const DiscreteEvent event = *event_it;
		EventKind kind = _sys->event_kind(_location,event);
		if (kind == URGENT || kind == PERMISSIVE) {

			RealScalarFunction guard = _sys->guard_function(_location,event);
			RealVectorFunction dynamic = _sys->dynamic_function(_location);

			tribool is_transition_taken = _is_transition_taken(guard,kind,dynamic,src_encl);

			if (possibly(is_transition_taken)) {
				const DiscreteLocation trg_location = _sys->target(_location,event);
				const EnclosureType trg_enclosure = _calculus->reset_step(
						_sys->reset_function(_location,event),src_encl);
				const Box trg_bbx = trg_enclosure.bounding_box();

				tribool is_outside_any_trg_invariant = _is_outside_any_invariant(trg_location,trg_bbx);

				if (possibly(!is_outside_any_trg_invariant)) {
					const LocalisedBox trg_hbbx(trg_location,trg_enclosure.bounding_box());

					tribool is_covered_backward = _is_covered_backward(trg_hbbx);

					tribool check_successful = (
							!is_outside_any_src_invariant &&
							is_transition_taken &&
							is_covered_backward &&
							!is_outside_any_trg_invariant);

					if (definitely(check_successful))
						return true;

					if (possibly(check_successful))
						result = indeterminate;
				}
			}
		}
	}

	return result;
}


tribool
BackwardDiscreteJumpSetChecker::
_is_covered_backward(const LocalisedBox& src_hbx) const {

	if (definitely(_starting_set->superset(src_hbx)))
		return true;

	if (definitely(_starting_set->disjoint(src_hbx)))
		return false;

	return indeterminate;
}


bool
BackwardDiscreteJumpSetChecker::
_has_feasible_transitions(const Box& src_bx) const
{
	const EnclosureType src_enclosure(src_bx);
	const Set<DiscreteEvent> events = _sys->events(_location);

	for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {

		const DiscreteEvent event = *event_it;
		EventKind kind = _sys->event_kind(_location,event);

		if (kind == URGENT || kind == PERMISSIVE) {

			DiscreteLocation trg_location = _sys->target(_location,event);
			const DenotableSetType trg_restriction = _starting_set->find(trg_location)->second;

			if (!trg_restriction.empty()) {

				const EnclosureType trg_enclosure = _calculus->reset_step(
						_sys->reset_function(_location,event),src_enclosure);

				Box trg_enclosure_bx = trg_enclosure.bounding_box();
				Box trg_restriction_bx = trg_restriction.bounding_box();

				if (possibly(trg_restriction_bx.overlaps(trg_enclosure_bx)))
					return true;
			}
		}
	}
	return false;
}


} // namespace Ariadne
