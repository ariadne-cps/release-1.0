/***************************************************************************
 *            verifier.h
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

/*! \file verifier.h
 *  \brief Methods for performing verification over hybrid systems.
 */

#ifndef ARIADNE_VERIFIER_H_
#define ARIADNE_VERIFIER_H_

#include "reachability_analyser.h"
#include "logging.h"
#include "verification_input.h"

namespace Ariadne {

// Keep this value in sync with the maximum verbosity level on the Verifier methods
const unsigned verifier_max_verbosity_level_used = 6;

class HybridReachabilityAnalyser;

enum DominanceSystem { DOMINATING_SYSTEM, DOMINATED_SYSTEM };

/** \brief Performs verification over reachable sets information. */
class Verifier
    : public Loggable
{
  public:

    typedef HybridAutomatonInterface SystemType;

  private:

    typedef boost::shared_ptr<HybridReachabilityAnalyser> AnalyserPtrType;

  private:

    /*! \brief Holds the cached result of an outer approximation of a system. */
    struct OuterApproximationCache
    {
      private:
    	bool _is_set;
    	HybridGridTreeSet _value;
      public:
    	OuterApproximationCache() : _is_set(false), _value(HybridGridTreeSet()) {}
    	virtual ~OuterApproximationCache() {}
    	bool is_set() const { return _is_set; }
    	HybridGridTreeSet get() const { return _value; }
    	void set(const HybridGridTreeSet& hgts) { ARIADNE_ASSERT(!_is_set); _is_set = true; _value = hgts; }
    	void reset() { _is_set = false; _value = HybridGridTreeSet(); }
    };

  private:
    mutable boost::shared_ptr< SystemType > _system;
    mutable boost::shared_ptr< HybridReachabilityAnalyser > _analyser;
    boost::shared_ptr< VerificationSettings > _settings;
    mutable std::string _plot_dirpath;

	/*! \brief Fields for caching outer approximations obtained during safety or dominance methods.
	 * \details Their presence is useful to refine the derivatives and consequently the grid on successive iterations. Set to mutable since their value is
	 * valid only transitively during one (iterative) verification method, therefore they do not add external state to the verifier. */
	mutable boost::shared_ptr< OuterApproximationCache > _safety_coarse_outer_approximation;
	mutable boost::shared_ptr< OuterApproximationCache > _dominating_coarse_outer_approximation;
	mutable boost::shared_ptr< OuterApproximationCache > _dominated_coarse_outer_approximation;

	/*! \brief Fields for holding the reachability restriction for reachability analyses.
	 * \details The latter two are mandatory since an analyser alternatively works on the dominating and dominated systems. */
	mutable HybridGridTreeSet _safety_reachability_restriction;
	mutable HybridGridTreeSet _dominating_reachability_restriction;
	mutable HybridGridTreeSet _dominated_reachability_restriction;

  public:

    //@{
    //! \name Constructors and destructors

    /*! \brief Virtual destructor */
    virtual ~Verifier();

    /*! \brief Default constructor. */
    Verifier();

    //@}

    //@{
    //! \name Methods to set and get the settings controlling the verification

    const VerificationSettings& settings() const { return *this->_settings; }
    VerificationSettings& settings() { return *this->_settings; }
    //@}

    //! \brief Overrides the Loggable method.
    //! \details Necessary to percolate to the internal analyser.
    void set_verbosity(int verbosity);

    //@{
    //! \name Safety methods

	/*! \brief Attempt to verify that the reachable set of a system starting in an initial_set remains in a safe region inside a domain.
	 * \details This is done in an iterative way by tuning the evolution/analysis parameters. The \a verInput contains all the information
	 * necessary for verification.
	 */
	tribool safety(SafetyVerificationInput& verInput) const;

	/**
	 * \brief Performs a parametric verification on a set of parameters \a params, by partitioning the parameters space.
	 * \details The \a logNumIntervalsPerParam variable determines how many times any parameter is split in two.
	 * The values in \a params are substituted into the system, the latter
	 * being restored to its initial conditions by the end of the method.
	 */
	std::list<ParametricOutcome> parametric_safety(
			SafetyVerificationInput& verInput,
			const RealParameterSet& params) const;

	//@}

	//@{
	//! \name Dominance methods

	/**
	 * \brief Performs dominance checking.
	 * \details Verifies if the \a dominating system dominates the \a dominated system. Dominance is intended in terms of
	 * the outer reachability of a dominating system being inside the BOUNDING BOX of the lower reachability of a dominated
	 * system minus its epsilon.
	 */
	tribool dominance(
			DominanceVerificationInput& dominating,
			DominanceVerificationInput& dominated) const;

	/**
	 * \brief Performs a parametric dominance checking on a set of parameters \a dominating_params of the \a dominating system, by partitioning
	 * the parameters space.
	 * \details The \a logNumIntervalsPerParam variable determines how many times any parameter is split in two.
     * The values in \a dominating_params are substituted into the \a dominating
	 * system alone, the latter being restored to its initial conditions by the end of the method.
	 */
	std::list<ParametricOutcome> parametric_dominance(
			DominanceVerificationInput& dominating,
			DominanceVerificationInput& dominated,
			const RealParameterSet& dominating_params) const;

	//@}

  private:

	//@{
	//! \name Safety methods

	/*! \brief Prove (once, i.e. for a given grid depth) that the the reachable set of \a system starting in \a initial_set
	 * definitely respects the \a safety_constraint.
	 * \details The \a constants are substituted into the system. */
	bool _safety_proving_once(
	        SafetyVerificationInput& verInput,
			const unsigned int& accuracy,
			const RealParameterSet& params) const;

    /*! \brief Returns the new initial grid tree set for a chain reach. */
    HybridGridTreeSet _reachability_refinement_starting_set(
            const HybridReachabilityAnalyser& analyser,
    		SystemType& system,
    		const HybridImageSet& initial_set,
    		const HybridConstraintSet& constraint_set,
    		const HybridGridTreeSet& reachability_restriction,
    		ContinuousEvolutionDirection direction) const;

	/*! \brief Prove (once, i.e. for a given grid depth) that the reachable set of \a system starting in \a initial_set
	 * does definitely DOES NOT respect the \a safety_constraint.
	 * \details The \a params are substituted into the system. */
	bool _safety_disproving_once(
	        SafetyVerificationInput& verInput,
	        const unsigned int& accuracy,
			const RealParameterSet& params) const;

    /*! \brief Attempt (once, i.e. for a given grid depth) to verify that the reachable set of \a system starting in \a initial_set
     * respects the \a safety_constraint.
     * \details The \a params are substituted into the system. */
    tribool _safety_once(
            SafetyVerificationInput& verInput,
            const unsigned int& accuracy,
			const RealParameterSet& constants) const;

	/*! \brief Performs iterative safety verification where \a parameter is substituted into the system.
	 */
	tribool _safety(
			SafetyVerificationInput& verInput,
			const RealParameter& parameter) const;

	/*! \brief Performs iterative safety verification where the singleton \a value is substituted into the system for the given \a constant.
	 */
	tribool _safety(
			SafetyVerificationInput& verInput,
			const RealParameter& parameter,
			const Float& value) const;

	/*! \brief Performs iterative safety verification, with \a params_to_substitute substituted into the system.
	 * \details The \a params are substituted in the system and are not allowed to be split */
	tribool _safety_nosplitting(
			SafetyVerificationInput& verInput,
			const RealParameterSet& params) const;

	//@}

	//@{
	//! \name Dominance methods

	/**
	 * \brief Performs dominance checking with \a constant substituted into the \a dominating system.
	 * \details Verifies if the \a dominating system dominates the \a dominated system.
	 */
	tribool _dominance(
			DominanceVerificationInput& dominating,
			DominanceVerificationInput& dominated,
			const RealParameter& param) const;

	/**
	 * \brief Performs dominance checking with \a constant substituted into the \a dominating system with a value of \a value.
	 * \details Verifies if the \a dominating system dominates the \a dominated system.
	 */
	tribool _dominance(
			DominanceVerificationInput& dominating,
			DominanceVerificationInput& dominated,
			const RealParameter& param,
			const Float& value) const;

	/*! \brief Helper function to perform dominance in the more general case when some \a constants are substituted into the dominating system. */
	tribool _dominance(
			DominanceVerificationInput& dominating,
			DominanceVerificationInput& dominated,
			const RealParameterSet& params) const;

	/*! \brief Performs the proving part of dominance checking.
	 * \details Tries proving only once, in respect to the current grid depth. */
	bool _dominance_proving_once(
			DominanceVerificationInput& dominating,
			DominanceVerificationInput& dominated,
			const RealParameterSet& params,
			const unsigned int& accuracy) const;

	/*! \brief Performs the disproving part of dominance checking.
	 * \details Tries disproving only once, in respect to the current grid depth. */
	bool _dominance_disproving_once(
			DominanceVerificationInput& dominating,
			DominanceVerificationInput& dominated,
			const RealParameterSet& params,
			const unsigned int& accuracy) const;

	/*! \brief Gets the flattened lower reach and the epsilon of the \a dominanceSystem type. */
	std::pair<GridTreeSet,Vector<Float> > _dominance_flattened_lower_reach_and_epsilon(
			DominanceVerificationInput& verInput,
			const RealParameterSet& params,
			DominanceSystem dominanceSystem,
			const unsigned int& accuracy) const;

	/*! \brief Gets the flattened outer reach of the \a dominanceSystem type. */
	GridTreeSet _dominance_flattened_outer_reach(
			DominanceVerificationInput& verInput,
			const RealParameterSet& params,
			DominanceSystem dominanceSystem,
			const unsigned int& accuracy) const;

	//@}

	//@{
	//! \name Other helper methods

    /*! \brief Obtains an analyser from the verification input, already tuned in respect to the other arguments. */
    AnalyserPtrType _get_tuned_analyser(
            const VerificationInput& verInput,
            const Set<Identifier>& locked_params_ids,
            const HybridGridTreeSet& outer_approx,
            int accuracy,
            bool EQUAL_GRID_FOR_ALL_LOCATIONS,
            Semantics semantics) const;

	/*! \brief Resets cached verification state information, for safety. */
	void _reset_safety_state() const;

	/*! \brief Resets cached verification state information, for dominance. */
	void _reset_dominance_state() const;

	/*! \brief Updates with \a reach the outer approximation or the reachability restriction.
	 * \details Which field is set depends on the current state of such variables: if no coarse outer
	 * approximation is set, then it is set, otherwise the reachability restriction is set (updated).
	 */
	void _update_safety_cached_reachability_with(const HybridGridTreeSet& reach) const;

	// Reached region plotting methods
	void _plot_dirpath_init(std::string basename) const;
	void _plot_reach(
			const HybridGridTreeSet& reach,
			string base_filename,
			int accuracy) const;
	void _plot_dominance(
			const HybridGridTreeSet& reach,
			DominanceSystem dominanceSystem,
			int accuracy,
			Semantics semantics) const;

	//@}

};

/*! \brief Splits the parameters to the maximum based on the \a tolerance
 *  \details The \a numIntervalsPerParam is the number of intervals to split for each parameter.
 *  \return The resulting split parameters sets.
 */
std::list<RealParameterSet> maximally_split_parameters(
		const RealParameterSet& params,
		const uint& maximum_parameter_depth);

}

#endif /* ARIADNE_VERIFIER_H_ */
