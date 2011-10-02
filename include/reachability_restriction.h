/***************************************************************************
 *            reachability_restriction.h
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file reachability_restriction.h
 *  \brief Allows a computationally flexible domain discretisation and set restricting operations.
 */

#ifndef REACHABILITY_RESTRICTION_H
#define REACHABILITY_RESTRICTION_H

#include "taylor_set.h"
#include "hybrid_set.h"

using namespace std;
using namespace Ariadne;

namespace Ariadne {

typedef TaylorSet ContinuousEnclosureType;
typedef HybridBasicSet<ContinuousEnclosureType> EnclosureType;

//! A structure for holding a restriction to reachability, with facilities to avoid unnecessary discretisations.
class ReachabilityRestriction
{
  public:

	//@{
	//! \name Constructors and destructors

	//! \brief Default constructor.
	ReachabilityRestriction(
			const HybridBoxes& domain,
			const HybridGrid& grid,
			int accuracy);

	//! \brief Copy constructor.
	ReachabilityRestriction(const ReachabilityRestriction& other);

	//! \brief Cloner.
	ReachabilityRestriction* clone() const { return new ReachabilityRestriction(*this); }

	virtual ~ReachabilityRestriction() { }

	//@}

    //@{
    //! \name Getters

	//! \brief Return the grid as exposed by the given domain.
	const HybridGrid& grid() const { return _grid; }

	//! \brief Return the hybrid space where the restriction is defined.
	HybridSpace space() const { return _grid.state_space(); }

	//! \brief Return the accuracy.
	//! \details This is not mandatory to expose this information, but useful
	//! for adopting the %HybridRestrictionSet as the accuracy reference.
	int accuracy() const { return _accuracy; }

	//! \brief Return an overapproximation of the bounding box of the restriction.
	HybridBoxes bounding_box() const;

	//@}

	//@{
	//! \name Predicates

	//! \brief Whether location \a q has already been discretised.
	bool has_discretised(DiscreteLocation q) const;

	//! \brief Whether location \a q is present in the space of the set.
	bool has_location(DiscreteLocation q) const;

	//! \brief Whether applying the restriction to \a set would affect it.
	bool restricts(const HybridDenotableSet& set) const;

	//@}

	//@{
	//! \name Setters

    //! \brief Update the content with the given \a set.
    //! \details Only the locations present in \a set are updated.
    void update(const HybridDenotableSet& set);

	//! \brief Refine the content at the given accuracy.
	//! \details Does not perform the operation on non-discretised locations or locations where the
	//! restriction set does not cross the domain border.
	void refine(int accuracy);

	//@}

	//@{
	//! \name Restricting operations

	//! \brief Apply the restriction to \a set.
	HybridDenotableSet apply_to(const HybridDenotableSet& set) const;

	//! \brief Apply the restriction to a list of enclosures
	std::list<EnclosureType> apply_to(const std::list<EnclosureType>& enclosures) const;

	//@}

	//@{
	//! \name Set computations

	//! \brief Compute the set obtained performing forward jumps from \a set in respect to the transitions in \a sys.
	HybridDenotableSet forward_jump_set(
			const HybridDenotableSet& set,
			const HybridAutomatonInterface& sys) const;

	//! \brief Compute the set obtained performing backward jumps from \a set in respect to the transitions in \a sys.
	HybridDenotableSet backward_jump_set(
			const HybridDenotableSet& set,
			const HybridAutomatonInterface& sys) const;

	//! \brief Computes the subset that is possibly feasible in respect to \a constraint.
	//! \details The result in general is a projection since \a constraint is allowed to be on a subspace of the hybrid
	//! space of the restriction.
	HybridDenotableSet possibly_feasible_projection(const HybridConstraintSet& constraint) const;

	//! \brief Computes the subset that is definitely feasible in respect to \a constraint.
	//! \details The result in general is a projection since \a constraint is allowed to be on a subspace of the hybrid
	//! space of the restriction.
	HybridDenotableSet definitely_feasible_projection(const HybridConstraintSet& constraint) const;

	//! \brief Computes the subset that is possibly infeasible in respect to \a constraint.
	//! \details The result in general is a projection since \a constraint is allowed to be on a subspace of the hybrid
	//! space of the restriction.
	HybridDenotableSet possibly_infeasible_projection(const HybridConstraintSet& constraint) const;

	//! \brief Computes the subset that is definitely infeasible in respect to \a constraint.
	//! \details The result in general is a projection since \a constraint is allowed to be on a subspace of the hybrid
	//! space of the restriction.
	HybridDenotableSet definitely_infeasible_projection(const HybridConstraintSet& constraint) const;

	//@}

  private:

    //! \brief Insert the discretisation at location q.
    //! \details It is assumed, not checked, that in location q no discretisation has already been performed.
    //! This operation will overwrite the existing content, in that case. It is const for compatibility with
    //! the operator[] const method, and it const-casts this to work the problem around.
    void _insert_domain_discretisation(DiscreteLocation q) const;

    //! \brief Get the bounding box for the location \a q.
	//! \details Does not discretise if not necessary: hence the result for non-discretised locations
	//! is a reasonable overapproximation.
	Box _bounding_box(DiscreteLocation q) const;

	//! \brief Whether there is a feasible transition from the restriction set at \a src_location towards \a trg_set.
	//! \details This is a rough check using only the bounding domain of the restriction.
	bool
	_has_feasible_transition(
			const DiscreteLocation& src_location,
			const HybridDenotableSet& trg_set,
			const HybridAutomatonInterface& sys) const;

	//! \brief Checks \a src_enclosure's jump sets from \a src_location. Adjoins all found jump sets to \a set_to_adjoin.
	void
	_adjoin_forward_jump_sets(
			const DiscreteLocation& src_location,
			const ContinuousEnclosureType& src_enclosure,
			const HybridAutomatonInterface& sys,
			HybridDenotableSet& result_set) const;

	//! \brief Checks if \a src_enclosure's jump set from \a src_location intersects \a trg_set. If it does, adjoins it to \a set_to_adjoin.
	void
	_adjoin_src_of_forward_jump_sets(
			const DiscreteLocation& src_location,
			const ContinuousEnclosureType& src_enclosure,
			const HybridAutomatonInterface& sys,
			const HybridDenotableSet& trg_set,
			HybridDenotableSet& set_to_adjoin) const;

	//! \brief Whether the transition specified by \a activation would be possible from \a source.
	//! \details It depends on the \a event_kind and the direction of the \a dynamic.
	bool
	_is_transition_feasible(
			const ScalarFunction& activation,
			EventKind event_kind,
			const VectorFunction& dynamic,
			const ContinuousEnclosureType& source) const;

	//! \brief Whether the box \bx is definitely outside the invariants in the given \a location according to \a sys.
	bool
	_is_outside_invariants(
			const DiscreteLocation& location,
			const Box& bx,
			const HybridAutomatonInterface& sys) const;

  private:

	// The domain that would be discretised into the related denotable set
	HybridBoxes _domain;
	// The grid used
	HybridGrid _grid;
	// The accuracy to be used when operating with the set
	int _accuracy;
	// The internal set
	HybridDenotableSet _set;
	// A calculus for activation/reset operations on sets
	boost::shared_ptr<CalculusInterface<TaylorModel> > _calculus;
};


//! \brief Whether the crossing of \a activation under \a dynamic is positive for a given \a set_bounds.
tribool
is_positively_crossing(
		const Box& set_bounds,
		const RealVectorFunction& dynamic,
		const RealScalarFunction& activation);

}

#endif /* REACHABILITY_RESTRICTION_H */
