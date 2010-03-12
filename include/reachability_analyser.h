/***************************************************************************
 *            reachability_analyser.h
 *
 *  Copyright  2006-10  Alberto Casagrande, Pieter Collins, Luca Geretti
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
 
/*! \file reachability_analyser.h
 *  \brief Methods for computing abstract reachable sets.
 */

#ifndef ARIADNE_REACHABILITY_ANALYSER_H
#define ARIADNE_REACHABILITY_ANALYSER_H

#include <boost/smart_ptr.hpp>

#include "hybrid_set_interface.h"
#include "evolver_interface.h"
#include "discretiser_interface.h"
#include "reachability_analyser_interface.h"

#include "orbit.h"
#include "grid_set.h"
#include "hybrid_set.h"

#include "discretiser.h"
#include "hybrid_evolver.h"

#include "logging.h"

namespace Ariadne {
 
template<class ES> class Orbit;

class DiscreteState;

template<class BS> class HybridBasicSet;
typedef HybridBasicSet<Box> HybridBox;

class HybridGrid;
class HybridGridCell;
class HybridGridTreeSet;

template<class ES> class HybridListSet;
template<class ES> class HybridDiscretiser;

/*! \brief A class for performing reachability analysis on a hybrid system.
	\details Log levels and the methods that are allowed to log output at such a level.
	<br>	1: Methods using the verify_iterative() method; 
	<br>	2: The verify_iterative() method;								
	<br>	3: The verify() and viable() methods;
	<br>	4: Remaining analysis methods implemented from reachability_analyser_interface.h;
	<br>	5: Internal analysis methods.
 */
class HybridReachabilityAnalyser
    : public Loggable
{
  private:
    boost::shared_ptr< DiscreteEvolutionParameters > _parameters;
	boost::shared_ptr< DiscreteEvolutionStatistics > _statistics;
    boost::shared_ptr< HybridDiscretiser<HybridEvolver::ContinuousEnclosureType> > _discretiser;
  public:
    typedef DiscreteEvolutionParameters EvolutionParametersType;
	typedef DiscreteEvolutionStatistics	EvolutionStatisticsType;
    typedef HybridAutomaton SystemType;
    typedef SystemType::StateSpaceType StateSpaceType;
    typedef SystemType::TimeType TimeType;
    typedef HybridGridTreeSet SetApproximationType;
    typedef HybridEvolver::EnclosureType EnclosureType;
    typedef HybridEvolver::ContinuousEnclosureType ContinuousEnclosureType;
  public:
    //@{
    //! \name Constructors and destructors
    /*! \brief Virtual destructor */
    virtual ~HybridReachabilityAnalyser();

    /*! \brief Construct from a method for evolving basic sets. */
    HybridReachabilityAnalyser(const HybridDiscretiser<HybridEvolver::ContinuousEnclosureType>& discretiser);

    /*! \brief Construct from evolution parameters and a method for evolving basic sets. */
    template<class HybridEnclosureType>
    HybridReachabilityAnalyser(const EvolutionParametersType& parameters,
                               const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver);

    template<class HybridEnclosureType>
    HybridReachabilityAnalyser(const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver);

    /*! \brief Make a dynamically-allocated copy. */
    virtual HybridReachabilityAnalyser* clone() const { return new HybridReachabilityAnalyser(*this); }
    //@}
  
    //@{ 
    //! \name Methods to set and get the parameters controlling the accuracy
    /*! \brief The parameters controlling the accuracy. */
    const EvolutionParametersType& parameters() const { return *this->_parameters; }
    /*! \brief A reference to the parameters controlling the accuracy. */
    EvolutionParametersType& parameters() { return *this->_parameters; }
    //@}

    //@{ 
    //! \name Methods to set and get the statistics related to the analyses
    /*! \brief The statistics stemming from the analyses. */
    const EvolutionStatisticsType& statistics() const { return *this->_statistics; }
    /*! \brief A reference to the statistics stemming from the analyses. */
    EvolutionStatisticsType& statistics() { return *this->_statistics; }
    //@}
  
    //@{
    //! \name Evaluation of systems on abstract sets
    /*! \brief Compute a lower-approximation to the set obtained by evolving \a system for \a time starting in \a initial_set. */
    virtual SetApproximationType lower_evolve(const SystemType& system,
                                              const HybridImageSet& initial_set, 
                                              const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time. */
    virtual SetApproximationType
    lower_reach(const SystemType& system, 
                const HybridImageSet& initial_set, 
                const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable and evolved sets of \a system starting in \a initial_set up to \a time. */
    virtual std::pair<SetApproximationType,SetApproximationType>
    lower_reach_evolve(const SystemType& system, 
                       const HybridImageSet& initial_set, 
                       const TimeType& time) const;
  
    /*! \brief Compute an approximation to the set obtained by iterating \a time times \a system starting in \a initial_set. */
    virtual SetApproximationType upper_evolve(const SystemType& system,
                                              const HybridImageSet& initial_set, 
                                              const TimeType& time) const;
  
    /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
    virtual SetApproximationType upper_reach(const SystemType& system,
                                             const HybridImageSet& initial_set, 
                                             const TimeType& timeType) const;
  
    /*! \brief Compute an approximation to the reachable and evolved sets of \a system starting in \a initial_set iterating at most \a time times. */
    virtual std::pair<SetApproximationType,SetApproximationType>
    upper_reach_evolve(const SystemType& system, 
                       const HybridImageSet& initial_set, 
                       const TimeType& time) const;
  
    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set. */
    virtual SetApproximationType chain_reach(const SystemType& system,
                                             const HybridImageSet& initial_set) const;
  
    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set and remaining in \a bounding_domain. \deprecated */
    virtual SetApproximationType chain_reach(const SystemType& system,
                                             const HybridImageSet& initial_set, 
                                             const HybridBoxes& bounding_domain) const;
  
    /*! \brief Compute an outer-approximation to the viability kernel of \a system within \a bounding_set. */
    virtual SetApproximationType viable(const HybridAutomaton& system,
                                        const HybridImageSet& bounding_set) const;
  
    /*! \brief Attempt to verify that the reachable set of \a system starting in \a initial_set remains in \a safe_set. */
    virtual tribool verify(const HybridAutomaton& system, 
                           const HybridImageSet& initial_set, 
                           const HybridImageSet& safe_set) const;
    //@}
  
  public:
    typedef HybridTime T;
    typedef HybridAutomaton Sys;
    typedef HybridListSet<Box> BxLS;
    typedef HybridGrid Gr;
    typedef HybridGridCell GC;
    typedef HybridGridTreeSet GCLS;
    typedef HybridGridTreeSet GTS;
    typedef HybridOpenSetInterface OpSI;
    typedef HybridOvertSetInterface OvSI;
    typedef HybridCompactSetInterface CoSI;
  public:
    // Helper functions for operators on lists of sets.
    GTS _upper_reach(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    GTS _upper_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
    std::pair<GTS,GTS> _upper_reach_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
  private:
    // Helper functions for approximating sets
};


template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver)
    : _parameters(new EvolutionParametersType())
	, _statistics(new EvolutionStatisticsType())
    , _discretiser(new HybridDiscretiser<typename HybridEnclosureType::ContinuousStateSetType>(evolver))
{
}


template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolutionParametersType& parameters,
                           const EvolverInterface<HybridAutomaton,HybridEnclosureType>& evolver)
    : _parameters(new EvolutionParametersType(parameters))
	, _statistics(new EvolutionStatisticsType())
    , _discretiser(new HybridDiscretiser<typename HybridEnclosureType::ContinuousStateSetType>(evolver))
{
}


} // namespace Ariadne


#endif // ARIADNE_REACHABILITY_ANALYSER_H
