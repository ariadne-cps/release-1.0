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
#include "discretiser_interface.h"
#include "reachability_analyser_interface.h"

#include "orbit.h"
#include "grid_set.h"
#include "hybrid_set.h"
#include "graphics.h"

#include "discretiser.h"
#include "hybrid_evolver.h"

#include "logging.h"
#include "parametric.h"

namespace Ariadne {
 
template<class ES> class Orbit;

class DiscreteLocation;

template<class BS> class HybridBasicSet;
typedef HybridBasicSet<Box> HybridBox;
typedef std::map<DiscreteLocation,Vector<Float> > HybridFloatVector;
typedef std::map<Identifier,int> ParameterIdIntMap;
typedef ParameterizableHybridAutomatonInterface SystemType;
typedef HybridEvolver::EnclosureType EnclosureType;
typedef HybridEvolver::ContinuousEnclosureType ContinuousEnclosureType;

class HybridGrid;
class HybridGridCell;
class HybridGridTreeSet;

template<class ES> class HybridListSet;
template<class ES> class HybridDiscretiser;

// Keep this value in sync with the maximum verbosity level on the Analyser methods
const unsigned analyser_max_verbosity_level_used = 7;

/*! \brief A class for performing reachability analysis on a hybrid system-
 */
class HybridReachabilityAnalyser
    : public Loggable
{
  private:
    boost::shared_ptr< DiscreteEvolutionSettings > _settings;
    boost::shared_ptr< HybridDiscretiser<HybridEvolver::ContinuousEnclosureType> > _discretiser;
  public:
    typedef DiscreteEvolutionSettings EvolutionSettingsType;
    typedef ParameterizableHybridAutomatonInterface SystemType;
    typedef SystemType::StateSpaceType StateSpaceType;
    typedef SystemType::TimeType TimeType;
    typedef HybridGridTreeSet SetApproximationType;
    typedef HybridEvolver::EnclosureType EnclosureType;
    typedef HybridEvolver::ContinuousEnclosureType ContinuousEnclosureType;
    typedef SetApproximationType (HybridReachabilityAnalyser::*OuterChainReachFuncPtr)(const SystemType&, const HybridImageSet&) const;
  public:
    //@{
    //! \name Constructors and destructors
    /*! \brief Virtual destructor */
    virtual ~HybridReachabilityAnalyser();

    /*! \brief Construct from a method for evolving basic sets. */
    HybridReachabilityAnalyser(const HybridDiscretiser<HybridEvolver::ContinuousEnclosureType>& discretiser);

    /*! \brief Construct from evolution parameters and a method for evolving basic sets. */
    template<class HybridEnclosureType>
    HybridReachabilityAnalyser(
    		const EvolutionSettingsType& parameters,
            const EvolverInterface<HybridAutomatonInterface,HybridEnclosureType>& evolver);

    template<class HybridEnclosureType>
    HybridReachabilityAnalyser(const EvolverInterface<HybridAutomatonInterface,HybridEnclosureType>& evolver);

    /*! \brief Make a dynamically-allocated copy. */
    virtual HybridReachabilityAnalyser* clone() const { return new HybridReachabilityAnalyser(*this); }
    //@}
  
    //@{ 
    //! \name Methods to set and get the settings controlling the accuracy
    /*! \brief The settings controlling the accuracy. */
    const EvolutionSettingsType& settings() const { return *this->_settings; }
    /*! \brief A reference to the settings controlling the accuracy. */
    EvolutionSettingsType& settings() { return *this->_settings; }
    //@}
  
    //@{

    //! \name Evaluation of systems on abstract sets
    /*! \brief Compute a lower-approximation to the set obtained by evolving \a system for \a time starting in \a initial_set. */
    virtual SetApproximationType lower_evolve(
    		const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time (discrete part only). */
    virtual SetApproximationType lower_reach(
    		const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute a lower-approximation to the reachable and evolved sets of \a system starting in \a initial_set up to \a time. */
    virtual std::pair<SetApproximationType,SetApproximationType> lower_reach_evolve(
    		const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute an approximation to the set obtained by iterating \a time times \a system starting in \a initial_set. */
    virtual SetApproximationType upper_evolve(
    		const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
    virtual SetApproximationType upper_reach(
    		const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& timeType) const;
  
    /*! \brief Compute an approximation to the reachable and evolved sets of \a system starting in \a initial_set iterating at most \a time times. */
    virtual std::pair<SetApproximationType,SetApproximationType> upper_reach_evolve(
    		const SystemType& system,
            const HybridImageSet& initial_set,
            const TimeType& time) const;
  
    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set. */
    virtual SetApproximationType chain_reach(
    		const SystemType& system,
            const HybridImageSet& initial_set) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set and remaining in \a bounding_domain. \deprecated */
    virtual SetApproximationType chain_reach(
    		const SystemType& system,
            const HybridImageSet& initial_set,
            const HybridBoxes& bounding_domain) const;

    /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set with a given \a direction, using
     * upper semantics; the method performs discretisation before transitions, then checks activations on the discretised cells.
     * \return The reach set. */
    virtual SetApproximationType outer_chain_reach(
    		SystemType& system,
			const HybridImageSet& initial_set,
			EvolutionDirection direction,
			const HybridGridTreeSet& reachability_restriction) const;

    virtual SetApproximationType outer_chain_reach(
    		SystemType& system,
    		const HybridGridTreeSet& initial,
    		EvolutionDirection direction,
    		const HybridGridTreeSet& reachability_restriction) const;

    /*! \brief Compute the epsilon lower bounds of \a system starting in \a initial_set.
     * \return The reach and the epsilon values. */
    virtual std::pair<HybridGridTreeSet,HybridFloatVector> lower_reach_and_epsilon(
    		SystemType& system,
			const HybridImageSet& initial_set,
			const HybridGridTreeSet& reachability_restriction) const;

    /*! \brief Compute the epsilon lower bounds of \a system starting in \a initial_set.
     * \details If \a constraint_set is not empty, it checks whether the reachable area does not satisfy the constraint,
     * raising an exception if this is the case.
     * \return The reach and the epsilon values. */
    virtual std::pair<HybridGridTreeSet,HybridFloatVector> lower_reach_and_epsilon(
    		SystemType& system,
			const HybridImageSet& initial_set,
			const HybridConstraintSet& constraint_set,
			const HybridGridTreeSet& reachability_restriction) const;
  
    /*! \brief Tunes the settings of the internal evolver. */
    void tuneEvolverSettings(
    		const SystemType& system,
    		const HybridFloatVector& hmad,
			uint maximum_grid_depth,
			Semantics semantics);

    //@}

  public:

	// The reduction in the number of logical cores used in multithreading (down from the maximum concurrency of the machine) (zero by default)
	uint free_cores;
  
  public:

    typedef HybridTime T;
    typedef HybridListSet<Box> BxLS;
    typedef HybridGrid Gr;
    typedef HybridGridCell GC;
    typedef HybridGridTreeSet GCLS;
    typedef HybridGridTreeSet GTS;
    typedef HybridOpenSetInterface OpSI;
    typedef HybridOvertSetInterface OvSI;
    typedef HybridCompactSetInterface CoSI;

  private:

    /*! \brief Gets the calculus interface from the hybrid evolver. */
    const CalculusInterface<TaylorModel>& _getCalculusInterface(Semantics semantics) const;

    /*! \brief Plots \a reach in \a plot_dirpath directory, where \a name_prefix as a prefix to the filename */
    void _plot_reach(
    		const HybridGridTreeSet& reach,
    		string plot_dirpath,
    		string name_prefix) const;

    // Helper functions for operators on lists of sets.
    GTS _upper_reach(
    		const SystemType& sys,
    		const GTS& set,
    		const T& time,
    		const int accuracy) const;

    std::pair<GTS,GTS> _upper_reach_evolve(
    		const SystemType& sys,
    		const GTS& set,
    		const T& time,
    		const int accuracy) const;

    std::pair<GTS,GTS> _upper_reach_evolve(
    		const SystemType& sys,
    		const list<EnclosureType>& initial_enclosures,
    		const T& time,
    		EvolutionDirection direction,
    		bool enable_premature_termination_on_blocking_event,
    		int accuracy) const;

    /*! \brief Performs outer chain reach calculation. */
    SetApproximationType _outer_chain_reach(
    		SystemType& system,
    		const std::list<EnclosureType>& initial_enclosures,
    		EvolutionDirection direction,
    		const HybridGridTreeSet& reachability_restriction) const;

    /*! \brief Performs outer chain reach calculation, where the constants of the \a system are assumed to be already splitted. */
    SetApproximationType _outer_chain_reach_splitted(
    		const SystemType& system,
    		const std::list<EnclosureType>& initial_enclosures,
    		EvolutionDirection direction,
    		const HybridGridTreeSet& reachability_restriction) const;

    /*! \brief Pushes into \a result_enclosures the enclosures from \a reachCells.
     * \details Ignores enclosures that lie outside the domain.
     */
    void _outer_chain_reach_forward_pushTargetCells(
    		const SystemType& system,
    		const HybridGridTreeSet& reachCells,
    		std::list<EnclosureType>& result_enclosures,
    		bool use_domain_checking) const;

    /*! \brief Pushes into \a result_enclosures the source enclosures from \a sourceCellsOverapprox that reach \a reachCells. */
    void _outer_chain_reach_backward_pushSourceCells(
    		const SystemType& system,
    		const HybridGridTreeSet& reachCells,
    		std::list<EnclosureType>& result_enclosures,
    		HybridGridTreeSet sourceCellsOverapprox) const;

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
			const HybridGridTreeSet& targetCells,
			const HybridGrid& grid,
			std::list<EnclosureType>& result_enclosures) const;

    /*! \brief Pushes the enclosures from the \a finalCells tree set into the \a result_enclosures list.
     */
    void _outer_chain_reach_pushLocalFinalCells(
    		const HybridGridTreeSet& finalCells,
    		std::list<EnclosureType>& result_enclosures,
    		bool use_domain_checking) const;

    /*! \brief Gets the lower reach and the epsilon for the \a system.
     * \details The \a constraint_set is checked: if not empty and its epsilon relaxation is not satisfied
     * for the current lower reach, an exception is raised. */
    std::pair<HybridGridTreeSet,HybridFloatVector> _lower_reach_and_epsilon(
    		const SystemType& system,
    		const HybridImageSet& initial_set,
			const HybridConstraintSet& constraint_set,
			const HybridGridTreeSet& reachability_restriction) const;

    /*! \brief Filters \a final_enclosures into \a initial_enclosures.
     * \details The procedure prunes a percentage of the enclosures based on \a adjoined_evolve_sizes and \a superposed_evolve_sizes. */
    void _filter_enclosures(
    		std::list<EnclosureType>& final_enclosures,
    		std::list<EnclosureType>& initial_enclosures,
    		const std::map<DiscreteLocation,uint>& adjoined_evolve_sizes,
    		const std::map<DiscreteLocation,uint>& superposed_evolve_sizes,
    		bool use_domain_checking) const;

    /*! \brief Gets the set of all the split intervals from the parameters of the \a system with a given \a tolerance.
     *  \details The calculation is performed over the domain with a splitting limit controlled by the \a tolerance, excluding
     *  those parameters present in the locked_parameters.
     *  Orders the list elements by first picking the leftmost subintervals, followed by the rightmost and then
     *  all the remaining from right to left. If no parameters to split are available, returns the original parameters.
     */
    std::list<RealParameterSet> _getSplitParametersIntervalsSet(
    		SystemType& system,
    		float tolerance) const;
};

template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(const EvolverInterface<HybridAutomatonInterface,HybridEnclosureType>& evolver)
	: _settings(new EvolutionSettingsType())
	, _discretiser(new HybridDiscretiser<typename HybridEnclosureType::ContinuousStateSetType>(evolver))
	, free_cores(0)

{
	ImageSetHybridEvolver& internal_evolver = dynamic_cast<ImageSetHybridEvolver&>(*this->_discretiser->evolver());
	internal_evolver.verb_tab_prefix = verb_tab_prefix + analyser_max_verbosity_level_used;
}


template<class HybridEnclosureType>
HybridReachabilityAnalyser::
HybridReachabilityAnalyser(
		const EvolutionSettingsType& parameters,
		const EvolverInterface<HybridAutomatonInterface,HybridEnclosureType>& evolver)
	: _settings(new EvolutionSettingsType(parameters))
	, _discretiser(new HybridDiscretiser<typename HybridEnclosureType::ContinuousStateSetType>(evolver))
	, free_cores(0)
{
	ImageSetHybridEvolver& internal_evolver = dynamic_cast<ImageSetHybridEvolver&>(*this->_discretiser->evolver());
	internal_evolver.verb_tab_prefix = verb_tab_prefix + analyser_max_verbosity_level_used;
}

/*! \brief Gets the minimum allowed widths of the cells of the \a grid under \a maximum_grid_depth */
HybridFloatVector min_cell_widths(
		const HybridGrid& grid,
		int maximum_grid_depth);

/*! \brief Generates a list of hybrid enclosures from the domain midpoints of a splitting of \a img_set,
 * where the image of each part has bounding box widths lower than the corresponding \a max_cell_widths */
list<EnclosureType> enclosures_from_split_domain_midpoints(
		const HybridImageSet img_set,
		const HybridFloatVector max_cell_widths);

/*! \brief Gets for each non-singleton parameter the factor determining the number of chunks its interval should be split into.
 *
 * \details Splits until the deviation of the derivatives is reasonably low in respect to the deviation calculated at the midpoint. This
 * limit value is expressed as a percentage using \a targetRatioPerc.
 *
 * @param system The system to get the parameters from.
 * @param targetRatioPerc The derivative widths ratio percentage to reach before termination.
 *
 * @return A split factor for each non-singleton parameter name of the \a system.
 */
ParameterIdIntMap getSplitFactorsOfParameters(
		SystemType& system,
		const Set<Identifier>& locked_params_ids,
		const Float& targetRatioPerc,
		const HybridBoxes& bounding_domain);

/*! \brief Gets the best parameter among the \a working_parameters of the \a system to split, in terms of
 * relative reduction of derivative widths compared to some \a referenceWidths.
 */
RealParameter getBestParameterToSplit(
		SystemType& system,
		const RealParameterSet& working_parameters,
		const HybridFloatVector& referenceWidths,
		const HybridBoxes& domain);

/*! \brief Helper function to get the maximum value of the derivative width ratio \f$ (w-w^r)/w_r \f$, where the \f$ w^r \f$ values
 * are stored in \a referenceWidths and the \f$ w \f$ values are obtained from the \a system.
 */
Float getMaxDerivativeWidthRatio(
		const SystemType& system,
		const HybridFloatVector& referenceWidths,
		const HybridBoxes& domain);

/*! \brief Helper function to get the widths of the derivatives from the \a system */
HybridFloatVector getDerivativeWidths(
		const SystemType& system,
		const HybridBoxes& domain);

/*! \brief Gets the set of all the midpoints of the split intervals in \a intervals_set. */
std::list<RealParameterSet> getSplitParametersMidpointsSet(const std::list<RealParameterSet>& intervals_set);

/*! \brief Splits the constant \a con into \a numParts parts. */
std::vector<RealParameter> split(
		const RealParameter& con,
		uint numParts);

/*! \brief Creates a set \a dest of all the possible combinations of split values from \a src. */
void fillSplitSet(
		const std::vector<std::vector<RealParameter> >& src,
		std::vector<std::vector<RealParameter> >::iterator col_it,
		std::vector<RealParameter>::iterator row_it,
		RealParameterSet s,
		std::list<RealParameterSet>& dest);

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

/*! \brief Set the maximum enclosure cell from the hybrid grid \a hgrid and the \a maximum_grid_depth. */
Vector<Float> getMaximumEnclosureCell(
		const HybridGrid& hgrid,
		int maximum_grid_depth);

/*! \brief Get the hybrid maximum integration step size, under the assumption that given the maximum derivatives \a hmad,
	all variables in a step must cover a length greater than a length determined by the \a hgrid with a given \a maximum_grid_depth.
	\details The actual result is scaled based on the \a semantics. */
std::map<DiscreteLocation,Float> getHybridMaximumStepSize(
		const HybridFloatVector& hmad,
		const HybridGrid& hgrid,
		int maximum_grid_depth,
		Semantics semantics);

/*! \brief Set the lock to grid time of \system.
	\details The value is taken as the maximum over the times required by any variable on any location to cover a distance equal to
	the domain width of the location, moving at the maximum absolute derivative.
	ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
Float getLockToGridTime(
		const SystemType& system,
		const HybridBoxes& domain);

/*! \brief Get the hybrid maximum absolute derivatives of \system given a previously computed outer approximation
 *  \a outer_approx_constraint and a domain \a domain_constraint.
 * \details ASSUMPTION: the continuous variables are preserved in order and quantity between discrete states. */
HybridFloatVector getHybridMaximumAbsoluteDerivatives(
		const SystemType& system,
		const HybridGridTreeSet& outer_approx_constraint,
		const HybridBoxes& domain_constraint);

/*! \brief Checks whether \a new_restriction is equal to \a old_restriction.
 * \details The routine assumes that \a new_restriction has been obtained from \a old_restriction, thus
 * eliminating the need for any check of the two tree sets.
 */
bool new_reachability_restriction_equals(
		const HybridGridTreeSet& new_restriction,
		const HybridGridTreeSet& old_restriction);

/*! \brief Turns the cells from the grid tree set \a cells into enclosures of smallest size in respect to a given
 * \a maximum_grid_depth. */
std::list<EnclosureType> cells_to_smallest_enclosures(
		HybridGridTreeSet cells,
		int maximum_grid_depth);

/*! \brief Restricts the \a enclosures to those possibly overlapping \a restriction */
std::list<EnclosureType> restrict_enclosures(
		const std::list<EnclosureType> enclosures,
		const HybridGridTreeSet& restriction);

/*! \brief Gets the cells from \a reach (inside \a reachability_restriction) that satisfy \a constraint
 * under a relaxation of the reachability given by \a eps. */
HybridGridTreeSet possibly_feasible_cells(
		const HybridGridTreeSet& reach,
		const HybridConstraintSet& constraint,
		const HybridFloatVector eps,
		HybridGridTreeSet reachability_restriction,
		int accuracy);


} // namespace Ariadne

#endif // ARIADNE_REACHABILITY_ANALYSER_H
