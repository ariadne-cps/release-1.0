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
#include "reachability_restriction.h"
#include "verification_input.h"

namespace Ariadne {

class VerifierSettings;
class Loggable;
class Interruptible;

enum DominanceSystem { DOMINATING_SYSTEM, DOMINATED_SYSTEM };

/** \brief Performs verification over reachable sets information. */
class Verifier
    : public Loggable
    , public Interruptible
{
  public:

    typedef HybridAutomatonInterface SystemType;
    typedef VerifierSettings SettingsType;
    typedef ReachabilityAnalyserInterface<SystemType>::SetApproximationType SetApproximationType;

  private:

    typedef boost::shared_ptr<ReachabilityAnalyserInterface<SystemType> > AnalyserPtrType;
    typedef boost::shared_ptr<ReachabilityRestriction> ReachabilityRestrictionPtr;

  private:
    boost::shared_ptr<SettingsType> _settings;
    mutable std::string _plot_dirpath;
    mutable time_t _start_time;

	/*! \brief "Stateless" fields for holding reachability restrictions between successive internal calls. */
	mutable ReachabilityRestrictionPtr _safety_restriction;
	mutable ReachabilityRestrictionPtr _dominating_restriction;
	mutable ReachabilityRestrictionPtr _dominated_restriction;

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

    const VerifierSettings& settings() const { return *this->_settings; }
    VerifierSettings& settings() { return *this->_settings; }

    //@}

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
			const RealParameterSet& params) const;

	/*! \brief Prove (once, i.e. for a given grid depth) that the reachable set of \a system starting in \a initial_set
	 * does definitely DOES NOT respect the \a safety_constraint.
	 * \details The \a params are substituted into the system. */
	bool _safety_disproving_once(
	        SafetyVerificationInput& verInput,
			const RealParameterSet& params) const;

    /*! \brief Attempt (once, i.e. for a given grid depth) to verify that the reachable set of \a system starting in \a initial_set
     * respects the \a safety_constraint.
     * \details The \a params are substituted into the system. */
    tribool _safety_once(
            SafetyVerificationInput& verInput,
			const RealParameterSet& params) const;

	/*! \brief Performs iterative safety verification, with \a params_to_substitute substituted into the system.
	 * \details The \a params are substituted in the system and are not allowed to be split */
	tribool _safety_nosplitting(
			SafetyVerificationInput& verInput,
			const RealParameterSet& params) const;

	/*! \brief Performs backward refinement of the safety reachability restriction. */
	void _safety_proving_once_backward_refinement(
	        const SystemType& sys,
	        const SetApproximationType& initial_set,
	        const HybridConstraintSet& constraint,
	        const RealParameterSet& params) const;

	/*! \brief Performs forward analysis for safety proving. */
	bool _safety_proving_once_forward_analysis(
	        const SystemType& sys,
	        const SetApproximationType& initial_set,
	        const HybridConstraintSet& constraint,
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
	 * \details Tries proving only once, in respect to the current restriction accuracy. */
	bool _dominance_proving_once(
			DominanceVerificationInput& dominating,
			DominanceVerificationInput& dominated,
			const RealParameterSet& params) const;

	/*! \brief Performs the disproving part of dominance checking.
	 * \details Tries disproving only once, in respect to the current restriction accuracy. */
	bool _dominance_disproving_once(
			DominanceVerificationInput& dominating,
			DominanceVerificationInput& dominated,
			const RealParameterSet& params) const;

	/*! \brief Gets the flattened lower reach and the epsilon of the \a dominanceSystem type. */
	std::pair<DenotableSetType,Vector<Float> > _dominance_flattened_lower_reach_and_epsilon(
			DominanceVerificationInput& verInput,
			const RealParameterSet& params,
			DominanceSystem dominanceSystem) const;

	/*! \brief Gets the flattened outer reach of the \a dominanceSystem type. */
	DenotableSetType _dominance_flattened_outer_reach(
			DominanceVerificationInput& verInput,
			const RealParameterSet& params,
			DominanceSystem dominanceSystem) const;

	//@}

	//@{
	//! \name Other helper methods

    /*! \brief Obtains an analyser from the system, already tuned in respect to the other arguments. */
    AnalyserPtrType _get_tuned_analyser(
            const SystemType& sys,
            const Set<Identifier>& locked_params_ids,
            const ReachabilityRestrictionPtr& restriction,
            const HybridConstraintSet& constraint_set,
            unsigned ADD_TAB_OFFSET) const;

	/*! \brief Initialises the safety restriction. */
	void _init_safety_restriction(const SafetyVerificationInput& verInput) const;

	/*! \brief Initialises the dominance restriction for the given \a dominanceSystem. */
	void _init_dominance_restriction(
			const DominanceVerificationInput& verInput,
			DominanceSystem dominanceSystem) const;

	/*! \brief Utility function generalising the safety and dominance cases. */
	void _init_restriction(
			const VerificationInput& verInput,
			ReachabilityRestrictionPtr& restriction,
			bool EQUAL_GRID_FOR_ALL_LOCATIONS) const;

	// Reached region plotting methods
	void _plot_dirpath_init(std::string basename) const;
	void _plot_reach(
			const SetApproximationType& reach,
			string base_filename,
			int accuracy) const;
	void _plot_dominance(
			const SetApproximationType& reach,
			DominanceSystem dominanceSystem,
			int accuracy,
			Semantics semantics) const;
	//@}

};


//! \brief Settings for controlling the verification flow.
class VerifierSettings {

    friend class Verifier;

  public:
    typedef uint UnsignedIntType;
    typedef double RealType;

  private:

    //! \brief Default constructor gives reasonable values.
    VerifierSettings();


  public:
    /*! \brief Whether the analysis results must be plotted. */
    bool plot_results;

    /*! \brief The maximum depth of parameter range splitting.
     * \details A value of zero means that the parameter space is not splitted at all. The total verification time
     * in this case is multiplied by 2^(number_of_parameters * maximum_parameter_depth). */
    uint maximum_parameter_depth;

    /*! \brief Whether to substitute midpoints of parameter boxes when proving.
     * \details Defaults to false. A value of true would not yield a formal result for the parameter box
     * but would be useful for quick pre-analysis. */
    bool use_param_midpoints_for_proving;

    /*! \brief Whether to substitute midpoints of parameter boxes when disproving.
     * \details Defaults to true. Indeed, if we use a value of false and successfully disprove, we gain no additional insight.
     * Choosing false has the benefit of exploring the whole parameter box, but the drawback of possibly be unable to successfully disprove at all
     * due to error radii. */
    bool use_param_midpoints_for_disproving;

    /*/ \brief Enable backward refinement of reachability when safety proving.
     * \details If disabled, it is possible for safety disproving to check whether the domain is incorrect. */
    bool enable_backward_refinement_for_safety_proving;
};

std::ostream& operator<<(std::ostream& os, const VerifierSettings& s);

/*! \brief Splits the parameters to the maximum based on the \a tolerance
 *  \details The \a numIntervalsPerParam is the number of intervals to split for each parameter.
 *  \return The resulting split parameters sets.
 */
std::list<RealParameterSet> maximally_split_parameters(
		const RealParameterSet& params,
		const uint& maximum_parameter_depth);

}

#endif /* ARIADNE_VERIFIER_H_ */
