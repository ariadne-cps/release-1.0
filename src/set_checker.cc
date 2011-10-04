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
		const HybridAutomatonInterface& sys,
		boost::shared_ptr<ReachabilityRestriction> restriction) :
		_location(loc),
		_sys(sys.clone()),
		_restriction(restriction),
		_calculus(new TaylorCalculus())
{
	HybridSpace sys_space = sys.state_space();

	ARIADNE_ASSERT_MSG(sys_space == restriction->grid().state_space(),
			"The system and restriction have different state spaces.");
	ARIADNE_ASSERT_MSG(sys_space.find(loc) != sys_space.end(),
			"The set location is not present in the system space.");
}


bool
DiscreteJumpSetCheckerBase::
_has_feasible_transitions(const Box& src_bx) const
{
	const ContinuousEnclosureType src_enclosure(src_bx);
	const Set<DiscreteEvent> events = _sys->events(_location);

	for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {

		const DiscreteEvent event = *event_it;
		EventKind kind = _sys->event_kind(_location,event);

		if (kind == URGENT || kind == PERMISSIVE) {

			DiscreteLocation trg_location = _sys->target(_location,event);

			if (!_restriction->has_empty(trg_location)) {

				const ContinuousEnclosureType target_encl = _calculus->reset_step(
						_sys->reset_function(_location,event),src_enclosure);

				Box trg_encl_bx = target_encl.bounding_box();
				Box trg_restr_bx = _restriction->bounding_box().at(trg_location);

				if (possibly(trg_restr_bx.overlaps(trg_encl_bx)))
					return true;
			}
		}
	}
	return false;
}


tribool
DiscreteJumpSetCheckerBase::
_is_transition_taken(
		const ScalarFunction& activation,
		EventKind event_kind,
		const VectorFunction& dynamic,
		const ContinuousEnclosureType& src_encl) const
{
	bool result = indeterminate;

	const bool is_urgent = (event_kind == URGENT);

	tribool is_guard_active = _calculus->active(VectorFunction(1,activation),src_encl);

	/*
	 * Definite cases:
	 *
	 * a) If the guard is definitely active and the transition is urgent, then we are definitely outside the related invariant
	 * b) If the transition is not urgent, we must perform the transition if the guard is possibly active
	 * c) If the transition is urgent and the guard is only possibly active, we refine by checking the crossing:
	 *    if it is definitely negative, then no transition is possible
	 */

	if (definitely(is_guard_active) && is_urgent) {
		result = false;
	} else if (possibly(is_guard_active) && !is_urgent) {
		result = true;
	} else if (possibly(is_guard_active) && is_urgent) {
		tribool positive_crossing = is_positively_crossing(src_encl.bounding_box(),dynamic,activation);
		if (definitely(!positive_crossing))
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
		if (kind == INVARIANT) {
			const ScalarFunction& activation = _sys->invariant_function(location,event);
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
    RealScalarFunction derivative=lie_derivative(activation,dynamic);
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
		boost::shared_ptr<ReachabilityRestriction> restriction,
		DiscreteEvent event) :
		DiscreteJumpSetCheckerBase(loc,sys,restriction),
		_event(event)
{
	Set<DiscreteEvent> events = sys.events(loc);
	ARIADNE_ASSERT_MSG(events.find(event) != events.end(),
			"The provided event " << event.name() << " is not present in the system at location " << loc.name());
}

tribool
ForwardDiscreteJumpSetChecker::
check(const Box& bx) const
{
	ARIADNE_NOT_IMPLEMENTED;
}

tribool
BackwardDiscreteJumpSetChecker::
check(const Box& bx) const
{
	tribool result = indeterminate;

	if (_is_outside_any_invariant(_location,bx))
		return false;

	if (!_has_feasible_transitions(bx))
		return false;

	const Set<DiscreteEvent> events = _sys->events(_location);

	ContinuousEnclosureType src_encl(bx);

	for (Set<DiscreteEvent>::const_iterator event_it = events.begin(); event_it != events.end(); ++event_it) {
		const DiscreteEvent event = *event_it;
		EventKind kind = _sys->event_kind(_location,event);
		if (kind == URGENT || kind == PERMISSIVE) {

			RealScalarFunction guard = _sys->guard_function(_location,event);
			RealVectorFunction dynamic = _sys->dynamic_function(_location);

			if (_is_transition_taken(guard,kind,dynamic,src_encl)) {
				const DiscreteLocation trg_location = _sys->target(_location,event);
				const ContinuousEnclosureType trg_encl = _calculus->reset_step(
						_sys->reset_function(_location,event),src_encl);
				const Box trg_bx = trg_encl.bounding_box();

				if (!_is_outside_any_invariant(trg_location,trg_bx)) {
					const HybridBox trg_hbx(trg_location,trg_encl.bounding_box());

					if (definitely(_restriction->superset(trg_hbx)))
						return true;
				}
			}
		}
	}

	return result;
}

} // namespace Ariadne
