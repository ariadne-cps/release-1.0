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
#include <cstdarg>
#include <config.h>

#include "hybrid_set_interface.h"
#include "evolver_interface.h"
#include "reachability_analyser_interface.h"

#include "orbit.h"
#include "grid_set.h"
#include "hybrid_set.h"
#include "graphics.h"

#include "hybrid_evolver.h"

#include "logging.h"
#include "parametric.h"

namespace Ariadne {

enum SplittingHalf { LEFT_HALF, RIGHT_HALF };
 
template<class ES> class Orbit;

class DiscreteLocation;
class HybridReachabilityAnalyserSettings;

template<class BS> class HybridBasicSet;
typedef HybridBasicSet<Box> HybridBox;
typedef std::map<DiscreteLocation,Vector<Float> > HybridFloatVector;
typedef std::map<Identifier,int> ParameterIdIntMap;
typedef HybridEvolver::EnclosureType EnclosureType;
typedef HybridEvolver::ContinuousEnclosureType ContinuousEnclosureType;
typedef boost::shared_ptr<HybridDenotableSet> HybridDenotableSetPtr;

class HybridGrid;
class HybridGridCell;
class HybridDenotableSet;
class ReachabilityRestriction;

template<class ES> class HybridListSet;
template<class ES> class HybridDiscretiser;

/*! \brief A class for performing reachability analysis on a hybrid system
 */
class HybridReachabilityAnalyser
    : public ReachabilityAnalyserInterface<HybridAutomatonInterface>
{
  public:
    typedef HybridReachabilityAnalyserSettings SettingsType;
    typedef ImageSetHybridEvolver::EnclosureType EnclosureType;
    typedef ImageSetHybridEvolver::ContinuousEnclosureType ContinuousEnclosureType;
  private:
    typedef boost::shared_ptr<EvolverInterface<SystemType,EnclosureType> > EvolverPtrType;
  private:
    boost::shared_ptr< SettingsType > _settings;
    boost::shared_ptr< ReachabilityRestriction > _restriction;
    mutable boost::shared_ptr< SystemType > _system;
  public:

    // The reduction in the number of logical cores used in multithreading (down from the maximum concurrency of the machine) (zero by default)
    uint free_cores;
  public:
    //@{
    //! \name Constructors and destructors

    /*! \brief Virtual destructor */
    virtual ~HybridReachabilityAnalyser();

    HybridReachabilityAnalyser(
    		const SystemType& system,
    		const ReachabilityRestriction& restriction);

    HybridReachabilityAnalyser(
    		const SystemType& system,
    		const HybridBoxes& domain,
    		int accuracy);

    /*! \brief Make a dynamically-allocated copy. */
    virtual HybridReachabilityAnalyser* clone() const { return new HybridReachabilityAnalyser(*this); }

    //@}

    //@{ 
    //! \name Methods to set and get the settings controlling the accuracy
    /*! \brief The settings controlling the accuracy. */
    const SettingsType& settings() const { return *this->_settings; }
    /*! \brief A reference to the settings controlling the accuracy. */
    SettingsType& settings() { return *this->_settings; }
    //@}
  
    //@{

    //! \name Evaluation of systems on abstract sets
    /*! \brief Compute a lower-approximation to the set obtained by evolving \a system for \a time starting in \a initial_set. */
    virtual SetApproximationType lower_evolve(
            const HybridBoundedConstraintSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time (discrete part only). */
    virtual SetApproximationType lower_reach(
            const HybridBoundedConstraintSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable and evolved sets of \a system starting in \a initial_set up to \a time. */
    virtual std::pair<SetApproximationType,SetApproximationType> lower_reach_evolve(
            const HybridBoundedConstraintSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute an approximation to the set obtained by iterating \a time times \a system starting in \a initial_set. */
    virtual SetApproximationType upper_evolve(
            const HybridBoundedConstraintSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
    virtual SetApproximationType upper_reach(
            const HybridBoundedConstraintSet& initial_set,
            const TimeType& timeType) const;
  
    /*! \brief Compute an approximation to the reachable and evolved sets of \a system starting in \a initial_set iterating at most \a time times. */
    virtual std::pair<SetApproximationType,SetApproximationType> upper_reach_evolve(
            const HybridBoundedConstraintSet& initial_set,
            const TimeType& time) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set with a given \a direction, using
     * upper semantics.
     * \return The reach set. */
    virtual SetApproximationType outer_chain_reach(
			const HybridBoundedConstraintSet& initial_set,
			ContinuousEvolutionDirection direction = DIRECTION_FORWARD) const;

    virtual SetApproximationType outer_chain_reach(
    		const SetApproximationType& initial_set,
    		ContinuousEvolutionDirection direction = DIRECTION_FORWARD) const;

    /*! \brief Compute the epsilon lower bounds of \a system starting in \a initial_set.
     * \return The reach and the epsilon values. */
    virtual std::pair<SetApproximationType,HybridFloatVector> lower_chain_reach_and_epsilon(
            const HybridBoundedConstraintSet& initial_set) const;

    //@}

    /*! \brief Tune the settings. */
    virtual void tune_settings(
            const Set<Identifier>& locked_params_ids,
            const HybridConstraintSet& constraint_set,
            unsigned free_cores,
            bool enable_lower_reach_restriction_check,
            Semantics semantics);

  public:

    typedef HybridTime T;
    typedef HybridListSet<Box> BxLS;
    typedef HybridGrid Gr;
    typedef HybridGridCell GC;
    typedef SetApproximationType GCLS;
    typedef SetApproximationType HDS;
    typedef HybridOpenSetInterface OpSI;
    typedef HybridOvertSetInterface OvSI;
    typedef HybridCompactSetInterface CoSI;

  private:

    /*! \brief Returns the accuracy, grid and domain used, taken from the internal restriction */
    int _accuracy() const;
    const HybridGrid& _grid() const;
    HybridBoxes _domain() const;

    /*! \brief Plots \a reach in \a plot_dirpath directory, where \a name_prefix as a prefix to the filename */
    void _plot_reach(
    		const SetApproximationType& reach,
    		string plot_dirpath,
    		string name_prefix) const;

    /*! \brief Obtains an evolver from the system, already tuned in respect to accuracy and semantics */
    EvolverPtrType _get_tuned_evolver(
            const SystemType& sys,
            unsigned ADD_TAB_OFFSET,
            Semantics semantics) const;

    std::pair<HDS,HDS> _upper_reach_evolve(
    		const SystemType& sys,
    		const SetApproximationType& initial_enclosures,
    		const T& time,
    		bool enable_premature_termination_on_blocking_event = false,
    		ContinuousEvolutionDirection direction = DIRECTION_FORWARD) const;

    /*! \brief Performs outer chain reach calculation, where the constants of the \a system are assumed to be already splitted. */
    SetApproximationType _outer_chain_reach_splitted(
    		const SystemType& system,
    		const SetApproximationType& initial,
    		ContinuousEvolutionDirection direction) const;

    /*! \brief Gets the lower reach and the epsilon for the \a system.
     * \details The \a constraint_set is checked: if not empty and its epsilon relaxation is not satisfied
     * for the current lower reach, an exception is raised. */
    std::pair<SetApproximationType,HybridFloatVector> _lower_chain_reach_and_epsilon(
    		const SystemType& system,
    		const HybridBoundedConstraintSet& initial_set) const;

    /*! \brief Checks whether \a reach with \a epsilon satisfies the constraint.
     * \details Throws ReachUnsatisfiesConstraintException if it doesn't.
     */
    void _lower_chain_reach_and_epsilon_constraint_check(
    		const SystemType& system,
    		const HDS& reach,
    		const HybridFloatVector& epsilon) const;

    /*! \brief Gets whether \a reach has a subset definitely infeasible in respect to the
     * constraint enlarged by \eps. */
    bool _has_eps_definitely_infeasible_subset(
    		const HybridDenotableSet& reach,
    		const HybridFloatVector& eps,
    		const HybridSpace& space) const;

    /*! \brief Filters \a final_enclosures into \a initial_enclosures.
     * \details The procedure prunes a percentage of the enclosures based on \a adjoined_evolve_sizes and \a superposed_evolve_sizes. */
    void _filter_enclosures(
    		std::list<EnclosureType>& final_enclosures,
    		std::list<EnclosureType>& initial_enclosures,
    		const std::map<DiscreteLocation,uint>& adjoined_evolve_sizes,
    		const std::map<DiscreteLocation,uint>& superposed_evolve_sizes) const;

    /*! \brief Gets the set of all the split intervals from the parameters of the system.*/
    std::list<RealParameterSet> _getSplitParameterSetList() const;

    /*! \brief Gets the maximum score corresponding to the derivative widths with no splitting. */
    Float _getDerivativeWidthsScore(
            const HybridFloatVector& hmad) const;

    /*! \brief Gets the scores corresponding to the derivative widths with \a param splitted. */
    std::pair<Float,Float> _getSplitDerivativeWidthsScores(
            const RealParameter& param,
            const HybridFloatVector& hmad) const;

    /*! \brief Updates the set \a working_scored_parameter_set_list by splitting onto the "best" parameter, or updates \a result_parameter_set_list if no splitting is possible. */
    /*! \details The \a hmad are the maximum derivative values with no splitting applied, hence they stay constant for successive calls of
     * this method.
     */
    void _updateSplitParameterSetLists(
            std::list<std::pair<Float,RealParameterSet> >& working_scored_parameter_set_list,
            std::list<RealParameterSet>& result_parameter_set_list,
            const Float& initial_score,
            const HybridFloatVector& hmad) const;

    /*! \brief Gets the ratio of derivative widths when substituting the given \a half of \a param.
     * \details The \a max_der_widths and mid_der_widths are the the maximum derivative widths with no splitting applied, and
     * the derivative widths for the midpoint of some parameter set (\a param included). */
    Float _getDerivativeWidthsRatio(
            const RealParameter& param,
            SplittingHalf half,
            const HybridFloatVector& max_der_widths,
            const HybridFloatVector& mid_der_widths) const;

    //! \brief Creates enclosures from the midpoints of the discretisation of \a initial_set.
    list<EnclosureType>
    _enclosures_from_discretised_initial_set_midpoints(const HybridBoundedConstraintSet initial_set) const;
};


//! \brief Settings for controlling the accuracy of discretised evolution methods and reachability analysis.
class HybridReachabilityAnalyserSettings {
    friend class HybridReachabilityAnalyser;
  public:
    typedef int IntType;
    typedef uint UnsignedIntType;
    typedef double RealType;
    typedef HybridAutomatonInterface SystemType;

  private:

    //! \brief Default constructor based on a system and some bounds.
    HybridReachabilityAnalyserSettings(
    		const SystemType& sys,
    		const HybridBoxes& domain);

  public:

    //! \brief The time after which an upper evolution or reachability analysis routine
    //! may approximate computed sets on a grid, in order to use previously cached
    //! integration results for the grid.
    //! \details
    //! Increasing this parameter improves the accuracy of the computations.
    //! Setting this parameter too low usually results in meaningless computations.
    //! As a rule of thumb, a typical system trajectory should move at least four
    //! times the grid size between locking to the grid. <br>
    //! For forced oscillators, this parameter should be set to the forcing time,
    //! or a multiple or fraction thereof.
    //! <br>
    //! This parameter is only used for continuous-time computation.
    RealType lock_to_grid_time;

    //! \brief The number of transitions before a discretisation.
    IntType lock_to_grid_steps;

    //! \brief Set the constraint set for reachability.
    //! \details Used for early termination of lower chain reachability. An empty constraint set implies no constraint at all.
    HybridConstraintSet constraint_set;

    //! \brief The parameters that must not be automatically split inside a system.
    Set<Identifier> locked_parameters_ids;

    //! \brief The target ratio of derivatives width to obtain when splitting parameters (must be > 0).
    RealType splitting_parameters_target_ratio;

    //! \brief Checks whether the restriction would affect the lower reach and issues an error if it does.
    //! \details Useful to validate the restriction used.
    bool enable_lower_reach_restriction_check;

    //! \brief Enable the pruning of the trajectories in lower semantics when too many.
    //! \details The pruning is done probabilistically.
    bool enable_lower_pruning;

};


std::ostream& operator<<(std::ostream& os, const HybridReachabilityAnalyserSettings& s);


/*! \brief Gets the minimum allowed widths of the cells of the \a grid under \a maximum_grid_depth */
HybridFloatVector min_cell_widths(
		const HybridGrid& grid,
		int maximum_grid_depth);

/*! \brief Removes from \a params those parameters having identifier in \a locked_params_ids. */
void remove_nonlocked_parameters(
        RealParameterSet& params,
        const Set<Identifier>& locked_params_ids);

/*! \brief Helper function to get the hybrid widths of the derivatives from the \a system. */
HybridFloatVector getHybridDerivativeWidths(
		const HybridReachabilityAnalyser::SystemType& system,
		const HybridBoxes& domain);

/*! \brief Helper function to get the widths of the derivatives from the \a system in a given location. */
Vector<Float> getDerivativeWidths(
        const HybridReachabilityAnalyser::SystemType& system,
        const DiscreteLocation& loc,
        const Box& bx);

/*! \brief Gets the set of all the midpoints of the split intervals in \a intervals_set. */
std::list<RealParameterSet> getMidpointsSet(const std::list<RealParameterSet>& intervals_set);

/*! \brief Get the hybrid grid given the maximum derivative \a hmad and the \a domain parameter, where the grid is chosen differently for each location.
 * \details The grid is chosen so that each cell is included into the domains for all locations. The \a equal_for_all_locations flag decides whether the
 * grid is the same for all locations. */
HybridGrid getHybridGrid(
		const HybridFloatVector& hmad,
		const HybridBoxes& domain,
		bool equal_for_all_locations);

/*! \brief Get the hybrid grid given the maximum derivative \a hmad and the \a domain parameter, where the grid is chosen differently for each location.
 * \details The grid is chosen so that each cell is included into the domains for all locations. */
HybridGrid getHybridGrid(
		const HybridFloatVector& hmad,
		const HybridBoxes& domain);

/*! \brief Get the grid given the maximum derivative \a hmad and the \a domain parameter, where the grid is chosen differently for each location.
 * \details The grid is chosen so that each cell is included into the domains for all locations. */
Grid getGrid(
		const HybridFloatVector& hmad,
		const HybridBoxes& domain);

/*! \brief Set the lock to grid time of \system.
	\details The value is taken as the maximum over the times required by any variable on any location to cover a distance equal to
	the domain width of the location, moving at the maximum absolute derivative.
	ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
Float getLockToGridTime(
		const HybridReachabilityAnalyser::SystemType& system,
		const HybridBoxes& domain);

/*! \brief Get the hybrid midpoint absolute derivatives of \system given a domain \a domain_constraint box.
 * \details ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
HybridFloatVector
getHybridMidpointAbsoluteDerivatives(
		const HybridReachabilityAnalyser::SystemType& sys,
		const HybridBoxes& bounding_domain);

/*! \brief Get the hybrid maximum absolute derivatives of \system given a previously computed outer approximation
 *  \a outer_approximation and a domain \a domain_constraint.
 * \details ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
HybridFloatVector
getHybridMaximumAbsoluteDerivatives(
		const HybridReachabilityAnalyser::SystemType& system,
		const HybridDenotableSetPtr& outer_approximation,
		const HybridBoxes& domain_constraint);

//! \brief Copies the \a reach set into enclosures.
std::list<EnclosureType> to_enclosures(const HybridDenotableSet& reach);

} // namespace Ariadne

#endif // ARIADNE_REACHABILITY_ANALYSER_H
