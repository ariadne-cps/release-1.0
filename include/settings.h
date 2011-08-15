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
  public:
    typedef uint UnsignedIntType;
    typedef double RealType;
    typedef HybridAutomatonInterface SystemType;

    //! \brief Default constructor gives reasonable values.
    EnclosedEvolutionSettings(const SystemType& sys);

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
  public:
    //! \brief The integer type.
    typedef int IntType;
    //! \brief The unsigned integer type.
    typedef uint UnsignedIntType;
    //! \brief The real type.
    typedef double RealType;
  
    //! \brief Default constructor gives reasonable values.
    DiscretisedEvolutionSettings();

    //! \brief The time after which infinite-time upper-evolution routines
    //! may approximate computed sets on a grid. 
    //! \details
    //! This value should be set to the time after which the transient
    //! behaviour has mostly died out. If there are no transients (i.e. the system evolves
    //! for a certain time and then stops), then this parameter should be set to a value 
    //! slightly higher than the maximum running time.
    //! <br> 
    //! Setting this value to the time at which transients die out can improve the
    //! speed and accuracy of the computations. 
    //!  <br> 
    //! This parameter is only used by chain_reach routine.
    RealType transient_time;

    //! \brief The number of discrete steps after which infinite-time upper evolution 
    //! routines may approximate computed sets on the grid. 
    //! \details
    //! Note that the transients are assumed to have died out after <em>either</em> 
    //! transient_time or transient_steps has been reached. 
    //! <br> 
    //! For example, if a hybrid system makes at most three steps before settling down
    //! into its limiting behaviour, then this parameter should be set to three. 
    //! <br> 
    //! Setting this value to the number of steps at which transients die out can improve the
    //! speed and accuracy of the computations. 
    //! <br> 
    //! This parameter is only used by chain_reach() routine.
    //! \sa #transient_time
    UnsignedIntType transient_steps;
    // (See the #transient_time parameter.) 

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

    //! \brief Set the depth used for approximation on a grid for the initial set in computations using lower semantics.
    //! \details
    //! Setting this value to \a d will on a grid with lengths \a l will result in the use of initial boxes
    //! with sides of length \f$l/2^{d}\f$.
    //! If the initial set is an open set, then this parameter may be unused; instead the initial sets are points,
    //! spaced according to the initial_grid_density.
    //!  <br> 
    //! Increasing this value increases the accuracy of the computation of lower evolution.
    //!  <br> 
    //! This parameter is only used in the lower_evolve() and lower_reach() routines.
    IntType initial_grid_depth;

    //! \brief Set the density of initial values on a grid for the initial set in computations using lower semantics.
    //! \details
    //! Setting this value to \a g will on a grid with lengths \a l will result in the use of one initial box
    //! (one simulation) for each cell of size \f$l/2^g\f$. If the #initial_grid_depth parameter is higher, then
    //! the initial sets will be smaller than the specified cell. 
    //!  <br> 
    //! Increasing this value increases the number of initial sets used in the computation of lower_evolve() and lower_reach().
    //! and decreases their spacing.
    //!  <br> 
    //! This parameter is only used in the lower_evolve() and lower_reach() routines.
    //! \internal Pieter: I don't like the name of this parameter very much, any other suggestions?
    IntType initial_grid_density;

    //! \brief Set the depth used for approximation on a grid for computations using upper semantics.
    //! \details
    //! Increasing this value increases the accuracy of the computation. 
    //!  <br> 
    //! This parameter is only used in upper_evolve(), upper_reach() and chain_reach() routines.
    IntType maximum_grid_depth;

    //! \brief Set the highest allowed value of maximum_grid_depth to be used by an analysis method that iteratively refines the grid depth.
    //! \details
    //! Increasing this value (greater or equal to zero) reduces the computation time when very high accuracy is required for verification,
	//! but can increase the computation time when very low accuracy is sufficient instead.
    //!  <br> 
    //! This parameter is only used in the iterative routines such as verify_iterative() and the *_parametric() routines.
	IntType lowest_maximum_grid_depth;

    //! \brief Set the highest allowed value of maximum_grid_depth to be used by an analysis method that iteratively refines the grid depth.
    //! \details
    //! Increasing this value increases the maximum accuracy of the computation for iterative methods. 
    //!  <br> 
    //! This parameter is only used in the verify_iterative() routine.
	IntType highest_maximum_grid_depth;

    //! \brief Set the maximum height used for approximation on a grid for chain reachability computations.
    //! \details
    //! Increasing this value increases domain over which computation is performed. 
    //!  <br> 
    //! This parameter is only used in the chain_reach() routine.
    IntType maximum_grid_height;

    //! \brief Set the allowed bounding domain for chain reachability computations.
    //! <br>
	//! This parameters is only used in the chain_reach() routines.
    HybridBoxes domain_bounds;

    //! \brief The grid to use.
    boost::shared_ptr<HybridGrid> grid;

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

/** 
 * Returns an appropriated grid for the system specified.
 *
 * This method returns a hybrid grid for the specified system. 
 * If \a *settings.grid is not empty, the method returns it, 
 * otherwise, it returns a hybrid grid for the state space 
 * of the provided system.
 *
 * @param system 
 *    The system for which we want a hybrid grid.
 * @param settings 
 *    The settings provided for the computation.
 * @return 
 *    An appropriate hybrid grid for the specified system. 
 */
inline
HybridGrid grid_for(const HybridAutomatonInterface& system,
                    const DiscretisedEvolutionSettings& settings)
{
  return ((settings.grid)->empty()?
            HybridGrid(system.state_space()):
            *(settings.grid));
}


//! \brief Settings for controlling the accuracy of continuous evolution methods.
class VerificationSettings {
  public:
	typedef uint UnsignedIntType;
	typedef double RealType;

	//! \brief Default constructor gives reasonable values.
	VerificationSettings();

    /*! \brief Whether the analysis results must be plotted. */
	bool plot_results;
	/*! \brief The maximum depth of parameter range splitting.
	 * \details A value of zero means that the parameter space is not splitted at all. */
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

    //! \brief Enable backward refinement of reachability when testing inclusion in a constraint.
 	//! \details The refinement itself is done only after a reachability restriction is available.
 	bool enable_backward_refinement_for_testing_inclusion;

 	/*! \brief Enable domain enforcing.
 	 *  \details If enforcing is disabled, the domain is only used to terminate reachability as soon as the domain is
 	 *  not respected. If enforcing is enabled, then reachability is strictly restricted inside the domain. In this case,
 	 *  there is an implicit assumption that the reachability is guaranteed to be inside the domain (the user must provide
 	 *  such guarantee by properly choosing the domain, otherwise the results are valid only on the domain itself).
 	 */
 	bool enable_domain_enforcing;
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
DiscretisedEvolutionSettings::DiscretisedEvolutionSettings() 
    : transient_time(0.0),
      transient_steps(0),
      lock_to_grid_time(1.0),
      lock_to_grid_steps(1),
      initial_grid_depth(10),
      initial_grid_density(8),
      maximum_grid_depth(6),
	  lowest_maximum_grid_depth(0),
	  highest_maximum_grid_depth(9),
      maximum_grid_height(16),
      grid(new HybridGrid()),
      splitting_constants_target_ratio(0.1),
	  enable_lower_pruning(true)
{ }

inline
VerificationSettings::VerificationSettings() :
		plot_results(false),
    	maximum_parameter_depth(3),
    	use_param_midpoints_for_proving(false),
    	use_param_midpoints_for_disproving(true),
    	enable_backward_refinement_for_testing_inclusion(true),
		enable_domain_enforcing(false)
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

       << ",\n  transient_time=" << s.transient_time
       << ",\n  transient_steps=" << s.transient_steps

       << ",\n  initial_grid_depth=" << s.initial_grid_depth
       << ",\n  initial_grid_density=" << s.initial_grid_density
       << ",\n  maximum_grid_depth=" << s.maximum_grid_depth
       << ",\n  lowest_maximum_grid_depth=" << s.lowest_maximum_grid_depth
       << ",\n  highest_maximum_grid_depth=" << s.highest_maximum_grid_depth
       << ",\n  maximum_grid_height=" << s.maximum_grid_height
       << ",\n  bounding_domain=" << s.domain_bounds
       << ",\n  grid=" << *s.grid
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
       << ",\n  maximum_parameter_depth=" << s.maximum_parameter_depth
       << ",\n  use_param_midpoints_for_proving=" << s.use_param_midpoints_for_proving
       << ",\n  use_param_midpoints_for_disproving=" << s.use_param_midpoints_for_disproving
       << ",\n  enable_fb_refinement_for_safety_proving=" << s.enable_backward_refinement_for_testing_inclusion
       << ",\n  enable_domain_enforcing=" << s.enable_domain_enforcing
       << "\n)\n";
    return os;
}

} //!namespace Ariadne

#endif //!ARIADNE_SETTINGS_H
