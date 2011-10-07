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

#ifndef ARIADNE_REACHABILITY_RESTRICTION_H
#define ARIADNE_REACHABILITY_RESTRICTION_H

#include "taylor_set.h"
#include "hybrid_set.h"

using namespace std;
using namespace Ariadne;

namespace Ariadne {

class HybridAutomatonInterface;

typedef TaylorSet EnclosureType;
typedef HybridBasicSet<EnclosureType> LocalisedEnclosureType;

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

	//! \brief Return an over-approximation of the bounding box of the restriction.
	HybridBoxes bounding_box() const;

	//! \brief Return an outer approximation of the domain box of the restriction.
	HybridBoxes outer_domain_box() const;

	//@}

	//@{
	//! \name Predicates

	//! \brief Whether location \a q is present in the space of the set.
	bool has_location(DiscreteLocation q) const;

	//! \brief Whether location \a q has already been discretised.
	bool has_discretised(DiscreteLocation q) const;

	//! \brief Whether the restriction is empty (i.e., it restricts to the empty set) at location \a q.
	bool has_empty(DiscreteLocation q) const;

	//! \brief Whether applying the restriction to \a set would affect it.
	bool restricts(const HybridDenotableSet& set) const;

	//@}

	//@{
	//! \name Setters

    //! \brief Update the content with the given \a set.
    //! \details Only the locations present in \a set are updated.
    void update_with(const HybridDenotableSet& set);

	//! \brief Refine the content at the given accuracy.
	//! \details Does not perform the operation on non-discretised locations or locations where the
	//! restriction set does not cross the domain border.
	void refine_at(int accuracy);

	//@}

	//@{
	//! \name Restricting operations

	//! \brief Apply the restriction to \a set.
	void apply_to(HybridDenotableSet& set) const;

	//! \brief Apply the restriction to a list of enclosures.
	//! \details This method differs in signature from apply_to since we need to produce a new list set.
	std::list<LocalisedEnclosureType> filter(const std::list<LocalisedEnclosureType>& enclosures) const;

	//@}

	//! \name Hybrid box checks

    tribool disjoint(const LocalisedBox& hbx) const;
    tribool overlaps(const LocalisedBox& hbx) const;
    tribool superset(const LocalisedBox& hbx) const;

    //@}

    //@{

	//! \name Set intersection and difference for constraint sets.

	//! \brief An outer approximation of the intersection with \a constraint_set.
	HybridDenotableSet outer_intersection_with(const HybridConstraintSet& constraint_set) const;

	//! \brief An inner approximation of the intersection with \a constraint_set.
	HybridDenotableSet inner_intersection_with(const HybridConstraintSet& constraint_set) const;

	//! \brief An outer approximation of the difference from \a constraint_set.
	HybridDenotableSet outer_difference_from(const HybridConstraintSet& constraint_set) const;

	//! \brief An inner approximation of the difference from \a constraint_set.
	HybridDenotableSet inner_difference_from(const HybridConstraintSet& constraint_set) const;

	//@{

	//! \name Feasibility projection for bounded constraint sets.

	//! \brief An outer approximation of the intersection with \a constraint_set.
	HybridDenotableSet outer_intersection_with(const HybridBoundedConstraintSet& constraint_set) const;

	//! \brief An inner approximation of the intersection with \a constraint_set.
	HybridDenotableSet inner_intersection_with(const HybridBoundedConstraintSet& constraint_set) const;

	//! \brief An outer approximation of the difference from \a constraint_set.
	HybridDenotableSet outer_difference_from(const HybridBoundedConstraintSet& constraint_set) const;

	//! \brief An inner approximation of the difference from \a constraint_set.
	HybridDenotableSet inner_difference_from(const HybridBoundedConstraintSet& constraint_set) const;

	//@}

	//@{
	//! \name Discrete jump set computations

	//! \brief Compute the set obtained performing forward jumps from \a set in respect to the transitions in \a sys.
	HybridDenotableSet forward_jump_set(
			const HybridDenotableSet& set,
			const HybridAutomatonInterface& sys) const;

	//! \brief Compute the set obtained performing backward jumps from \a set in respect to the transitions in \a sys.
	HybridDenotableSet backward_jump_set(
			const HybridDenotableSet& set,
			const HybridAutomatonInterface& sys) const;

	//@}

	std::ostream& write(std::ostream&) const;

  private:

    //! \brief Insert the discretisation at \a q.
    //! \details It is assumed, not checked, that in location q no discretisation has already been performed.
    //! This operation will overwrite the existing content, in that case. It is const for compatibility with
    //! the operator[] const method, and it const-casts this to work the problem around.
    void _insert_domain_discretisation(DiscreteLocation q) const;

    //! \brief Get the bounding box for the location \a q.
	//! \details Does not discretise if not necessary: hence the result for non-discretised locations
	//! is a reasonable overapproximation.
	Box _bounding_box(DiscreteLocation q) const;

	//! \brief Get an outer approximation of the domain box at \a q.
	//! \details This is a rough over-approximation based on the grid lengths.
	Box _outer_domain_box(DiscreteLocation q) const;

  private:

	// The domain that would be discretised into the related denotable set
	HybridBoxes _domain;
	// The grid used
	HybridGrid _grid;
	// The accuracy to be used when operating with the set
	int _accuracy;
	// The internal set
	HybridDenotableSet _set;
};


}

#endif /* ARIADNE_REACHABILITY_RESTRICTION_H */
