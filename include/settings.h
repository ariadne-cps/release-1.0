/***************************************************************************
 *            settings.h
 *
 *  Copyright  2007-11  Davide Bresolin, Alberto Casagrande, Pieter Collins, 
 *                      Luca Geretti
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
 
/*! \file settings.h
 *  \brief Settings for controlling the accuracy of evaluation methods.
 */

#ifndef ARIADNE_SETTINGS_H
#define ARIADNE_SETTINGS_H

#include <cstddef>
#include <boost/smart_ptr.hpp>

#include "grid_set.h"
#include "hybrid_set.h"
#include "variables.h"
#include "hybrid_automaton_interface.h"

namespace Ariadne {

class DiscreteLocation;

//! \brief Settings for controlling the accuracy of evolution methods on enclosure sets.
class EnclosedEvolutionSettings {
	friend class ImageSetHybridEvolver;
  public:
    typedef uint UnsignedIntType;
    typedef double RealType;
    typedef HybridAutomatonInterface SystemType;

  protected:

    //! \brief Default constructor gives reasonable values.
    EnclosedEvolutionSettings(const SystemType& sys);

  public:

    //! \brief The maximum allowable step size for integration, different for each location.
    //! Decreasing the values increases the accuracy of the computation.
    std::map<DiscreteLocation,RealType> hybrid_maximum_step_size;

    //! \brief The maximum allowable cell of a basic set during integration. 
    //! Decreasing the volume of the cell increases the accuracy of the computation of an over-approximation. 
    Vector<RealType> maximum_enclosure_cell;
    
    //! \brief Enable subdivision of basic sets (false by default).
    bool enable_subdivisions;

    //! \brief Terminate evolution if basic sets became too large (true by default).
	//! \details In the case of upper semantics, if true and no subdivisions are present, the set is put into the final sets. In the case of lower semantics, the set is discarded.
    bool enable_premature_termination_on_enclosure_size;
};


//! \brief Settings for controlling the accuracy of discretised evolution methods and reachability analysis.
class DiscretisedEvolutionSettings {
	friend class HybridReachabilityAnalyser;
  public:
    typedef int IntType;
    typedef uint UnsignedIntType;
    typedef double RealType;
    typedef HybridAutomatonInterface SystemType;
  
  protected:

    //! \brief Default constructor based on a system.
    DiscretisedEvolutionSettings(const SystemType& sys);

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
    //! \details Applied only to the chain reach routines. Assumed as not used if not assigned.
    //! (while, on the contrary, an empty reachable set would restrict any set to the empty set).
    boost::shared_ptr<HybridGridTreeSet> reachability_restriction;

    //! \brief The grid to use.
    HybridGrid grid;

    //! \brief The parameters that must not be automatically split inside a system.
    Set<Identifier> locked_parameters_ids;

    //! \brief The target ratio of derivatives width to obtain when splitting constants.
    RealType splitting_constants_target_ratio;

	//! \brief Enable the pruning of the trajectories when too large (false by default).
    //! \details The pruning is done probabilistically.
	//! <br>
    //! This parameter is used only under lower semantics.
	bool enable_lower_pruning;

};


//! \brief Settings for controlling the accuracy of continuous evolution methods.
class VerificationSettings {
  public:
	typedef uint UnsignedIntType;
	typedef double RealType;

	//! \brief Default constructor gives reasonable values.
	VerificationSettings();

    /*! \brief Whether the analysis results must be plotted. */
	bool plot_results;

	/*! \brief The time (in seconds) under which a verification outcome should be obtained.
	 * \details The verifier would stop iterative refinement (yielding indeterminate) as soon as
	 * it surpasses this value. */
	uint time_limit_for_outcome;

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

    /*/ \brief Enable backward refinement of reachability when testing inclusion in a constraint.
 	 * \details The refinement itself is done only after a reachability restriction is available. */
 	bool enable_backward_refinement_for_testing_inclusion;
};

inline
EnclosedEvolutionSettings::EnclosedEvolutionSettings(const SystemType& sys)
    : maximum_enclosure_cell(Vector<RealType>(sys.state_space().begin()->second,2.0)),
      enable_subdivisions(false),
      enable_premature_termination_on_enclosure_size(true)
{
	HybridSpace hspace(sys.state_space());
	for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
		hybrid_maximum_step_size[hs_it->first] = 1.0;
	}
}

inline
DiscretisedEvolutionSettings::DiscretisedEvolutionSettings(const SystemType& sys)
    : lock_to_grid_time(1.0),
      lock_to_grid_steps(1),
      maximum_grid_depth(6),
      domain_bounds(unbounded_hybrid_boxes(sys.state_space())),
      constraint_set(),
      reachability_restriction(),
      grid(HybridGrid(sys.state_space())),
      splitting_constants_target_ratio(0.1),
	  enable_lower_pruning(true)
{
}

inline
VerificationSettings::VerificationSettings() :
		plot_results(false),
		time_limit_for_outcome(10),
    	maximum_parameter_depth(3),
    	use_param_midpoints_for_proving(false),
    	use_param_midpoints_for_disproving(true),
    	enable_backward_refinement_for_testing_inclusion(true)
{ }


inline
std::ostream& 
operator<<(std::ostream& os, const EnclosedEvolutionSettings& s) 
{
    os << "ContinuousEvolutionSettings"
       << ",\n  hybrid_maximum_step_size=" << s.hybrid_maximum_step_size
       << ",\n  maximum_enclosure_cell=" << s.maximum_enclosure_cell
       << ",\n  enable_subdivisions=" << s.enable_subdivisions
       << ",\n  enable_premature_termination_on_enclosure_size=" << s.enable_premature_termination_on_enclosure_size
       << "\n)\n";
    return os;
}


inline
std::ostream& 
operator<<(std::ostream& os, const DiscretisedEvolutionSettings& s) 
{
    os << "DiscreteEvolutionSettings"
       << "(\n  lock_to_grid_steps=" << s.lock_to_grid_steps
       << ",\n  lock_to_grid_time=" << s.lock_to_grid_time
       << ",\n  maximum_grid_depth=" << s.maximum_grid_depth
       << ",\n  bounding_domain=" << s.domain_bounds
       << ",\n  grid=" << s.grid
       << ",\n  splitting_constants_target_ratio=" << s.splitting_constants_target_ratio
       << ",\n  enable_lower_pruning=" << s.enable_lower_pruning
       << "\n)\n";
    return os;
}


inline
std::ostream&
operator<<(std::ostream& os, const VerificationSettings& s)
{
    os << "VerificationSettings"
       << "(\n  plot_results=" << s.plot_results
       << ",\n  time_limit_for_outcome" << s.time_limit_for_outcome
       << ",\n  maximum_parameter_depth=" << s.maximum_parameter_depth
       << ",\n  use_param_midpoints_for_proving=" << s.use_param_midpoints_for_proving
       << ",\n  use_param_midpoints_for_disproving=" << s.use_param_midpoints_for_disproving
       << ",\n  enable_backward_refinement_for_testing_inclusion=" << s.enable_backward_refinement_for_testing_inclusion
       << "\n)\n";
    return os;
}

} //!namespace Ariadne

#endif //!ARIADNE_SETTINGS_H
