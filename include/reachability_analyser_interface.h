/***************************************************************************
 *            reachability_analyser_interface.h
 *
 *  Copyright  2006-8  Pieter Collins
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

/*! \file reachability_analyser_interface.h
 *  \brief Interface for performing reachability analysis.
 */

#ifndef ARIADNE_REACHABILITY_ANALYSER_INTERFACE_H
#define ARIADNE_REACHABILITY_ANALYSER_INTERFACE_H

#include <boost/smart_ptr.hpp>

#include "set_interface.h"
#include "hybrid_set_interface.h"
#include "logging.h"

namespace Ariadne {

/*! \brief Interface for computing reachable sets of a dynamic system.
 */
template<class SYS> class ReachabilityAnalyserInterface
    : public Loggable
{
  public:  
    //! \brief The type of the system.
    typedef SYS SystemType;
    //! \brief The type used to define the elapsed evolution time for the system type.
    typedef typename SystemType::TimeType TimeType;
    //! \brief The type used to describe the state space of the system evolution.
    typedef typename SystemType::StateSpaceType StateSpaceType;
    //! \brief The type used to describe the type used for concrete approximations to sets.
    typedef typename StateSpaceType::SetApproximationType SetApproximationType;
    //! \brief The type used to pass around references to set approximations.
    typedef typename boost::shared_ptr<SetApproximationType> SetApproximationPtrType;
  public:
    //! \brief Virtual destructor.
    virtual ~ReachabilityAnalyserInterface() { }
    
    //@{
    //! \name Evaluation of maps on abstract sets

    //! \name Evaluation of systems on abstract sets
    /*! \brief Compute a lower-approximation to the set obtained by evolving \a system for \a time starting in \a initial_set. */
    virtual SetApproximationType lower_evolve(
            const HybridImageSet& initial_set,
            const TimeType& time) const = 0;

    /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time (discrete part only). */
    virtual SetApproximationType lower_reach(
            const HybridImageSet& initial_set,
            const TimeType& time) const = 0;

    /*! \brief Compute a lower-approximation to the reachable and evolved sets of \a system starting in \a initial_set up to \a time. */
    virtual std::pair<SetApproximationType,SetApproximationType> lower_reach_evolve(
            const HybridImageSet& initial_set,
            const TimeType& time) const = 0;

    /*! \brief Compute an approximation to the set obtained by iterating \a time times \a system starting in \a initial_set. */
    virtual SetApproximationType upper_evolve(
            const HybridImageSet& initial_set,
            const TimeType& time) const = 0;

    /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
    virtual SetApproximationType upper_reach(
            const HybridImageSet& initial_set,
            const TimeType& timeType) const = 0;

    /*! \brief Compute an approximation to the reachable and evolved sets of \a system starting in \a initial_set iterating at most \a time times. */
    virtual std::pair<SetApproximationType,SetApproximationType> upper_reach_evolve(
            const HybridImageSet& initial_set,
            const TimeType& time) const = 0;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set with a given \a direction, using
     * upper semantics.
     * \return The reach set. */
    virtual SetApproximationType outer_chain_reach(
            const HybridImageSet& initial_set,
            ContinuousEvolutionDirection direction = DIRECTION_FORWARD) const = 0;

    virtual SetApproximationType outer_chain_reach(
            const SetApproximationType& initial_set,
            ContinuousEvolutionDirection direction = DIRECTION_FORWARD) const = 0;

    /*! \brief Compute the epsilon lower bounds of \a system starting in \a initial_set.
     * \return The reach and the epsilon values. */
    virtual std::pair<SetApproximationType,HybridFloatVector> lower_chain_reach_and_epsilon(
            const HybridImageSet& initial_set) const = 0;

    //@}
    
    //@{
    //! \name Tuning and utilities
    
    /*! \brief Tune the settings.
     * \details Much of the arguments do not need to be defined, being either empty or null.
     */
    virtual void tune_settings(
            const Set<Identifier>& locked_params_ids,
            const HybridConstraintSet& constraint_set,
            unsigned free_cores,
            bool enable_lower_reach_restriction_check,
            Semantics semantics) = 0;

    //@}
    
};


} // namespace Ariadne




#endif // ARIADNE_ANALYSER_INTERFACE_H
