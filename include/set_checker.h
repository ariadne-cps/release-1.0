/***************************************************************************
 *            set_checker.h
 *
 *  Copyright  2011  Luca Geretti
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

/*! \file set_checker.h
 *  \brief Classes for performing checks on sets.
 */

#ifndef ARIADNE_SET_CHECKER_H
#define ARIADNE_SET_CHECKER_H

#include "tribool.h"
#include <boost/shared_ptr.hpp>
#include "hybrid_automaton_interface.h"

namespace Ariadne {

class Box;
class TaylorModel;
template<class T> class CalculusInterface;
class ReachabilityRestriction;
class HybridDenotableSet;
class TaylorSet;
template<class T> class HybridBasicSet;
class ConstraintSet;

typedef TaylorSet ContinuousEnclosureType;
typedef HybridBasicSet<ContinuousEnclosureType> EnclosureType;


//! \brief Interface for performing checks on sets.
class SetCheckerInterface
{
  public:
	//! \brief Checks whether the bounding box of a set satisfies some properties.
	//! \details Returns true if all points in \a bx pass the check, false if no point
	//! in \a bx passes the check, indeterminate otherwise.
	virtual tribool check(const Box& bx) const = 0;
};


//! \brief Constraint set satisfaction.
class ConstraintSetChecker : public SetCheckerInterface
{
  protected:

	//! The constraint set to check against.
	boost::shared_ptr<ConstraintSet> _constraint;

  public:

	//! \brief Assign the constraint.
    ConstraintSetChecker(const ConstraintSet& constraint);

    //! \brief Virtual destructor.
    virtual ~ConstraintSetChecker() { }

    //! \brief Check whether the image of the constraint function on \a bx is inside the constraint codomain.
    virtual tribool check(const Box& bx) const;
};

//! \brief Class for checking against forward jump sets of discrete evolution jumps.
class DiscreteJumpSetCheckerBase : public SetCheckerInterface
{
  protected:

	//! The source location for the jump
	DiscreteLocation _location;
	//! The system defining the transitions
    boost::shared_ptr<HybridAutomatonInterface> _sys;
    //! A reachability restriction
    boost::shared_ptr<ReachabilityRestriction> _restriction;
    //! A calculus for performing activation and reset operations on transitions
    boost::shared_ptr<CalculusInterface<TaylorModel> > _calculus;

  public:

    //! \brief Assign the basic fields.
    DiscreteJumpSetCheckerBase(
    		const DiscreteLocation& loc,
    		const HybridAutomatonInterface& sys,
    		boost::shared_ptr<ReachabilityRestriction> restriction);

    //! \brief Virtual destructor.
    virtual ~DiscreteJumpSetCheckerBase() { }

  protected:

	//! \brief Whether there is a feasible transition from \a src_bx at the \a _src_location.
	//! \details This is a rough check using only the bounding domain of the restriction.
	bool _has_feasible_transitions(const Box& src_bx) const;

	//! \brief Whether the transition specified by \a activation would be possible from \a source.
	//! \details It depends on the \a event_kind and the direction of the \a dynamic.
	tribool _is_transition_taken(
			const ScalarFunction& activation,
			EventKind event_kind,
			const VectorFunction& dynamic,
			const ContinuousEnclosureType& source) const;

	//! \brief Whether the box \bx is outside any invariant in the given \a location.
	tribool _is_outside_any_invariant(
			const DiscreteLocation& location,
			const Box& bx) const;

	//! \brief Whether the crossing of \a activation under \a dynamic is positive for a given \a set_bounds.
	tribool _is_positively_crossing(
			const Box& set_bounds,
			const RealVectorFunction& dynamic,
			const RealScalarFunction& activation) const;

};


//! \brief Handles the case of forward jumps.
class ForwardDiscreteJumpSetChecker : public DiscreteJumpSetCheckerBase
{
  protected:

	//! The event to take into account
	DiscreteEvent _event;

  public:

	//! \brief Constructor accepting an \a event too.
	ForwardDiscreteJumpSetChecker(
    		const DiscreteLocation& loc,
    		const HybridAutomatonInterface& sys,
    		boost::shared_ptr<ReachabilityRestriction> restriction,
			DiscreteEvent event);

	//! \brief Check whether \a bx allows a target set for a specific event.
	virtual tribool check(const Box& bx) const;
};


//! \brief Handles the case of backward jumps.
class BackwardDiscreteJumpSetChecker : public DiscreteJumpSetCheckerBase
{
  public:

	//! \brief Check whether \a bx is reached by backward transitions from the location in the object.
	virtual tribool check(const Box& bx) const;
};


} // namespace Ariadne




#endif // ARIADNE_SET_CHECKER_H
