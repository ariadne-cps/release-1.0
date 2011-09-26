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
    mutable boost::shared_ptr< SystemType > _system;
  public:

    // The reduction in the number of logical cores used in multithreading (down from the maximum concurrency of the machine) (zero by default)
    uint free_cores;
  public:
    //@{
    //! \name Constructors and destructors

    /*! \brief Virtual destructor */
    virtual ~HybridReachabilityAnalyser();

    HybridReachabilityAnalyser(const SystemType& system);

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
            const HybridImageSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time (discrete part only). */
    virtual SetApproximationType lower_reach(
            const HybridImageSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable and evolved sets of \a system starting in \a initial_set up to \a time. */
    virtual std::pair<SetApproximationType,SetApproximationType> lower_reach_evolve(
            const HybridImageSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute an approximation to the set obtained by iterating \a time times \a system starting in \a initial_set. */
    virtual SetApproximationType upper_evolve(
            const HybridImageSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
    virtual SetApproximationType upper_reach(
            const HybridImageSet& initial_set,
            const TimeType& timeType) const;
  
    /*! \brief Compute an approximation to the reachable and evolved sets of \a system starting in \a initial_set iterating at most \a time times. */
    virtual std::pair<SetApproximationType,SetApproximationType> upper_reach_evolve(
            const HybridImageSet& initial_set,
            const TimeType& time) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set with a given \a direction, using
     * upper semantics.
     * \return The reach set. */
    virtual SetApproximationType outer_chain_reach(
			const HybridImageSet& initial_set,
			ContinuousEvolutionDirection direction = DIRECTION_FORWARD) const;

    virtual SetApproximationType outer_chain_reach(
    		const SetApproximationType& initial_set,
    		ContinuousEvolutionDirection direction = DIRECTION_FORWARD) const;

    /*! \brief Compute the epsilon lower bounds of \a system starting in \a initial_set.
     * \return The reach and the epsilon values. */
    virtual std::pair<SetApproximationType,HybridFloatVector> lower_chain_reach_and_epsilon(
            const HybridImageSet& initial_set) const;

    //@}

    /*! \brief Tune the settings. */
    virtual void tune_settings(
            const HybridBoxes& domain,
            const Set<Identifier>& locked_params_ids,
            const SetApproximationPtrType& reach_outer_approx,
            const SetApproximationPtrType& reach_restriction,
            const HybridConstraintSet& constraint_set,
            bool EQUAL_GRID_FOR_ALL_LOCATIONS,
            int accuracy,
            unsigned free_cores,
            Semantics semantics);

    /*! \brief Produces the starting cells from an enclosure set. */
    virtual SetApproximationType initial_cells_set(const HybridImageSet& initial_enclosure_set) const;

    /*! \brief Produces the starting cells from a constraint set. */
    virtual SetApproximationType initial_cells_set(const HybridConstraintSet& initial_constraint_set) const;
  
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

    /*! \brief Plots \a reach in \a plot_dirpath directory, where \a name_prefix as a prefix to the filename */
    void _plot_reach(
    		const SetApproximationType& reach,
    		string plot_dirpath,
    		string name_prefix) const;

    /*! \brief Obtains an evolver from the system, already tuned in respect to accuracy and semantics */
    EvolverPtrType _get_tuned_evolver(
            const SystemType& sys,
            int accuracy,
            unsigned ADD_TAB_OFFSET,
            Semantics semantics) const;

    std::pair<HDS,HDS> _upper_reach_evolve(
    		const SystemType& sys,
    		const HDS& set,
    		const T& time,
    		const int accuracy) const;

    std::pair<HDS,HDS> _upper_reach_evolve(
    		const SystemType& sys,
    		const list<EnclosureType>& initial_enclosures,
    		const T& time,
    		bool enable_premature_termination_on_blocking_event,
    		ContinuousEvolutionDirection direction,
    		int accuracy) const;

    /*! \brief Performs outer chain reach calculation. */
    SetApproximationType _outer_chain_reach(
    		const std::list<EnclosureType>& initial_enclosures,
    		ContinuousEvolutionDirection direction) const;

    /*! \brief Performs outer chain reach calculation, where the constants of the \a system are assumed to be already splitted. */
    SetApproximationType _outer_chain_reach_splitted(
    		const SystemType& system,
    		const std::list<EnclosureType>& initial_enclosures,
    		ContinuousEvolutionDirection direction) const;

    /*! \brief Pushes into \a result_enclosures the enclosures from \a reachCells.
     * \details Ignores enclosures that lie outside the domain.
     */
    void _outer_chain_reach_forward_pushTargetCells(
    		const SystemType& system,
    		const SetApproximationType& reachCells,
    		std::list<EnclosureType>& result_enclosures,
    		bool use_domain_checking) const;

    /*! \brief Pushes into \a result_enclosures the source enclosures from \a sourceCellsOverapprox that reach \a reachCells. */
    void _outer_chain_reach_backward_pushSourceCells(
    		const SystemType& system,
    		const SetApproximationType& reachCells,
    		std::list<EnclosureType>& result_enclosures) const;

    /*! \brief Checks whether a box \a bx is outside any invariant from \a invariants. */
    bool _outer_chain_reach_isOutsideInvariants(
    		const SystemType& system,
    		const DiscreteLocation& location,
    		const Box& bx) const;

    /*! \brief Checks whether the transition of kind \a event_kind, with a given \a activation, is feasible for enclosure \a source under dynamic \a dynamic.
     * \details By feasibility we mean that, under upper semantics, it would be taken by \a source.
     */
    bool _is_transition_feasible(
    		const ScalarFunction& activation,
    		EventKind event_kind,
    		const VectorFunction& dynamic,
    		const ContinuousEnclosureType& source,
    		Semantics semantics) const;

    /*! \brief Pushes the enclosures from the \a source enclosure into the \a destination enclosure list, for all \a transitions.
     */
    void _outer_chain_reach_forward_pushTargetEnclosures(
    		const SystemType& system,
    		const DiscreteLocation& sourceLocation,
    		const ContinuousEnclosureType& sourceEnclosure,
			const HybridGrid& grid,
			std::list<EnclosureType>& result_enclosures,
			bool use_domain_checking) const;

    void _outer_chain_reach_backward_pushSourceEnclosures(
    		const SystemType& system,
    		const DiscreteLocation& sourceLocation,
    		const ContinuousEnclosureType& sourceEnclosure,
			const SetApproximationType& targetCells,
			const HybridGrid& grid,
			std::list<EnclosureType>& result_enclosures) const;

    /*! \brief Pushes the enclosures from the \a finalCells tree set into the \a result_enclosures list.
     */
    void _outer_chain_reach_pushLocalFinalCells(
    		const SetApproximationType& finalCells,
    		std::list<EnclosureType>& result_enclosures,
    		bool use_domain_checking) const;

    /*! \brief Gets the lower reach and the epsilon for the \a system.
     * \details The \a constraint_set is checked: if not empty and its epsilon relaxation is not satisfied
     * for the current lower reach, an exception is raised. */
    std::pair<SetApproximationType,HybridFloatVector> _lower_chain_reach_and_epsilon(
    		const SystemType& system,
    		const HybridImageSet& initial_set) const;

    /*! \brief Filters \a final_enclosures into \a initial_enclosures.
     * \details The procedure prunes a percentage of the enclosures based on \a adjoined_evolve_sizes and \a superposed_evolve_sizes. */
    void _filter_enclosures(
    		std::list<EnclosureType>& final_enclosures,
    		std::list<EnclosureType>& initial_enclosures,
    		const std::map<DiscreteLocation,uint>& adjoined_evolve_sizes,
    		const std::map<DiscreteLocation,uint>& superposed_evolve_sizes,
    		bool use_domain_checking) const;

    /*! \brief Gets the set of all the split intervals from the parameters of the system.*/
    std::list<RealParameterSet> _getSplitParameterSetList() const;

    /*! \brief Gets the maximum score corresponding to the derivative widths with no splitting. */
    Float _getDerivativeWidthsScore(
            const HybridFloatVector& hmad,
            const SetApproximationType& discretised_domain) const;

    /*! \brief Gets the scores corresponding to the derivative widths with \a param splitted. */
    std::pair<Float,Float> _getSplitDerivativeWidthsScores(
            const RealParameter& param,
            const HybridFloatVector& hmad,
            const SetApproximationType& discretised_domain) const;

    /*! \brief Updates the set \a working_scored_parameter_set_list by splitting onto the "best" parameter, or updates \a result_parameter_set_list if no splitting is possible. */
    /*! \details The \a hmad are the maximum derivative values with no splitting applied, hence they stay constant for successive calls of
     * this method. The \a discretised_domain is used to evaluate the derivatives in different cells of the discretised domain.
     */
    void _updateSplitParameterSetLists(
            std::list<std::pair<Float,RealParameterSet> >& working_scored_parameter_set_list,
            std::list<RealParameterSet>& result_parameter_set_list,
            const Float& initial_score,
            const HybridFloatVector& hmad,
            const SetApproximationType& discretised_domain) const;

    /*! \brief Gets the ratio of derivative widths when substituting the given \a half of \a param.
     * \details The \a max_der_widths and mid_der_widths are the the maximum derivative widths with no splitting applied, and
     * the derivative widths for the midpoint of some parameter set (\a param included). */
    Float _getDerivativeWidthsRatio(
            const RealParameter& param,
            SplittingHalf half,
            const HybridFloatVector& max_der_widths,
            const HybridFloatVector& mid_der_widths) const;
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

    //! \brief Default constructor based on a system.
    HybridReachabilityAnalyserSettings(const SystemType& sys);

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

    //! \brief The time after which an evolver may approximate computed sets on a grid,
    //! in order to use previously cached results for the grid.
    //! \details
    //! Increasing this parameter may improve the accuracy of the computations.
    //! If there is recurrence in the system, then this parameter should be set to
    //! the average recurrence time, if known.
    //!  <br>
    //! This parameter is only used for discrete-time computation.
    UnsignedIntType lock_to_grid_steps;

    //! \brief Set the depth used for approximation on a grid for computations using upper semantics.
    //! \details
    //! Increasing this value increases the accuracy of the computation.
    //!  <br>
    //! This parameter is only used in upper_evolve(), upper_reach() and chain_reach() routines.
    IntType maximum_grid_depth;

    //! \brief Set the allowed bounding domain for chain reachability computations.
    //! \details Defaults to an unbounded box. Since it is also used to tune the evolver, it could be necessary to provide an explicit bounded value for it.
    HybridBoxes domain_bounds;

    //! \brief Set the constraint set for reachability.
    //! \details Used for early termination of lower chain reachability. An empty constraint set implies no constraint at all.
    HybridConstraintSet constraint_set;

    //! \brief Set the restriction for reachability.
    //! \details Be aware that such restriction is used for upper semantics strictly (tree set intersection), while for
    //! lower semantics the bounding boxes of the cells are employed instead. This permits to have different grids for different semantics.
    //! Hence, the analyser grid must match the reachability restriction grid for upper semantics only.
    boost::shared_ptr<HybridDenotableSet> reachability_restriction;

    //! \brief The grid to use.
    HybridGrid grid;

    //! \brief The parameters that must not be automatically split inside a system.
    Set<Identifier> locked_parameters_ids;

    //! \brief The target ratio of derivatives widthd to obtain when splitting parameters (must be > 0).
    RealType splitting_parameters_target_ratio;

    //! \brief Enable the pruning of the trajectories when too large.
    //! \details The pruning is done probabilistically.
    //! <br>
    //! This parameter is used only under lower semantics.
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

/*! \brief Generates a list of hybrid enclosures from the domain midpoints of a splitting of \a img_set,
 * where the image of each part has bounding box widths lower than the corresponding \a max_cell_widths */
list<EnclosureType> enclosures_from_split_domain_midpoints(
		const HybridImageSet img_set,
		const HybridFloatVector max_cell_widths);

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

/*! \brief Splits \a target_encl for location \a target_loc, storing the result in \a initial_enclosures.
 * \details The function is recursive.
 */
void pushSplitTargetEnclosures(
		std::list<EnclosureType>& initial_enclosures,
		const DiscreteLocation& target_loc,
		const ContinuousEnclosureType& target_encl,
	    const Vector<Float>& minTargetCellWidths,
		const Box& target_domain_constraint,
		bool use_domain_checking);

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

/*! \brief Get the hybrid maximum absolute derivatives of \system given a previously computed outer approximation
 *  \a outer_approximation and a domain \a domain_constraint.
 * \details ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
HybridFloatVector getHybridMaximumAbsoluteDerivatives(
		const HybridReachabilityAnalyser::SystemType& system,
		const HybridDenotableSetPtr& outer_approximation,
		const HybridBoxes& domain_constraint);

/*! \brief Checks whether \a new_restriction is equal to \a old_restriction.
 * \details The routine assumes that \a new_restriction has been obtained from \a old_restriction, thus
 * eliminating the need for any check of the two tree sets.
 */
bool new_reachability_restriction_equals(
		const HybridDenotableSet& new_restriction,
		const HybridDenotableSet& old_restriction);

/*! \brief Turns the cells from the grid tree set \a cells into enclosures of smallest size in respect to a given
 * \a maximum_grid_depth. */
std::list<EnclosureType> cells_to_smallest_enclosures(
		HybridDenotableSet cells,
		int maximum_grid_depth);

/*! \brief Restricts the \a enclosures to those possibly overlapping \a restriction */
std::list<EnclosureType> restrict_enclosures(
		const std::list<EnclosureType> enclosures,
		const HybridDenotableSet& restriction);

/*! \brief Gets the cells from \a reach (inside \a reachability_restriction) that satisfy \a constraint
 * under a relaxation of the reachability given by \a eps. */
HybridDenotableSet possibly_feasible_cells(
		const HybridDenotableSet& reach,
		const HybridConstraintSet& constraint,
		const HybridFloatVector eps,
		HybridDenotableSet reachability_restriction,
		int accuracy);


} // namespace Ariadne

#endif // ARIADNE_REACHABILITY_ANALYSER_H
