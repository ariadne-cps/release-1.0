/***************************************************************************
 *            hybrid_automaton.h
 *
 *  Copyright  2004-8  Alberto Casagrande, Pieter Collins
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

/*! \file hybrid_automaton.h
 *  \brief Main hybrid system class.
 */

#ifndef ARIADNE_HYBRID_AUTOMATON_H
#define ARIADNE_HYBRID_AUTOMATON_H

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "function.h"
#include "discrete_location.h"
#include "discrete_event.h"
#include "variables.h"
#include "hybrid_automaton_interface.h"

namespace Ariadne {


class HybridTime;
class HybridSpace;
class HybridSet;
class HybridGrid;

class DiscreteMode;
class DiscreteTransition;
class HybridAutomaton;

class ScalarFunction;
class VectorFunction;
class Grid;


/*! \brief A discrete mode of a hybrid automaton, comprising continuous evolution given by a vector field
 * within and invariant constraint set.
 *
 * A %DiscreteMode can only be created using the new_mode() method in
 * the %HybridAutomaton class.
 *
 * \sa \link Ariadne::HybridAutomaton \c HybridAutomaton \endlink, \link Ariadne::DiscreteTransition \c DiscreteTransition \endlink
 */
class DiscreteMode {
    friend class HybridAutomaton;
  private:

    // The discrete mode's discrete location.
    DiscreteLocation _location;

    // The discrete mode's vector field.
    VectorFunction _dynamic;
    // The discrete mode's invariants.
    std::map< DiscreteEvent, VectorFunction > _invariants;

  public:
    //! \brief The mode's discrete location.
    DiscreteLocation location() const {
        return this->_location; }

    //! \brief The discrete mode's dynamic (a vector field).
    const VectorFunction& dynamic() const {
        return this->_dynamic; }

    //! \brief The discrete mode's invariants.
    const std::map< DiscreteEvent, VectorFunction >& invariants() const {
        return this->_invariants; }

    //! \brief The discrete mode's invariants, converted to scalar functions.
    std::map< DiscreteEvent, ScalarFunction > scalar_invariants() const {
        std::map<DiscreteEvent,ScalarFunction> result;
        for(std::map<DiscreteEvent,VectorFunction>::const_iterator iter=this->_invariants.begin();
            iter!=this->_invariants.end(); ++iter)
        {
            ARIADNE_ASSERT_MSG(iter->second.result_size()==1,"Invariant "<<*iter<<" is not scalar.");
            result[iter->first]=iter->second[0];
        }
        return result;
    }

	/*! \brief Substitute the parameter \a param, if present, on the invariants and dynamic functions. */
	void substitute(const RealParameter& param) {
		this->_dynamic.substitute(param);
		for (std::map<DiscreteEvent,VectorFunction>::iterator it=this->_invariants.begin();it!=this->_invariants.end();it++)
			it->second.substitute(param);
	}

	/*! \brief Get the parameters from the dynamics and invariants */
	RealParameterSet parameters() const {
		RealParameterSet result = this->_dynamic.parameters();
		for (std::map<DiscreteEvent,VectorFunction>::const_iterator it=this->_invariants.begin();it!=this->_invariants.end();it++) {
			RealParameterSet new_parameters = it->second.parameters();
			result.insert(new_parameters.begin(),new_parameters.end());
		}
		return result;
	}

    //! \brief The dimension of the discrete mode.
    uint dimension() const;

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream& os) const;

  private:
    // Construct discrete mode.
    //
    // \param id is the identifier of the mode.
    // \param dynamic is the mode's vector field.
    // \param invariants is the mode's invariants.
    DiscreteMode(DiscreteLocation location,
                 const VectorFunction& dynamic);

    // Construct from objects managed by shared pointers (for internal use)
    DiscreteMode(DiscreteLocation location,
                 const VectorFunction dynamic,
                 const std::vector< VectorFunction >& invariants);

};


std::ostream& operator<<(std::ostream& os, const DiscreteMode& dm);

inline bool operator<(const DiscreteMode& mode1, const DiscreteMode& mode2) {
    return mode1.location() < mode2.location(); }




/*! \brief A discrete transition of a hybrid automaton, representing an instantaneous
 * jump from one discrete mode to another, governed by an activation set and a reset map.
 *
 * A %DiscreteTransition can only be created using the new_transition() method in
 * the %HybridAutomaton class.
 *
 * An invariant is modelled by a discrete transition with negative event id and null reset pointer.
 *
 * \sa \link Ariadne::HybridAutomaton \c HybridAutomaton \endlink, \link Ariadne::DiscreteMode \c DiscreteMode \endlink
 */
class DiscreteTransition
{
    friend class HybridAutomaton;
  private:
    // \brief The discrete transition's identificator.
    DiscreteEvent _event;

    // \brief The source of the discrete transition.
    DiscreteLocation _source;

    // \brief The target of the discrete transition.
    DiscreteLocation _target;

    // \brief The activation region of the discrete transition.
    VectorFunction _activation;

    // \brief The reset of the discrete transition.
    VectorFunction _reset;

    EventKind _kind;

    // \brief Whether or not the transition is forced.
    bool _forced;

  public:

    //! \brief The discrete event associated with the discrete transition.
    DiscreteEvent event() const {
        return this->_event; }

    //! \brief The source mode of the discrete transition.
    DiscreteLocation source() const {
        return this->_source; }

    //! \brief The target of the discrete transition.
    DiscreteLocation target() const {
        return this->_target; }

	/*! \brief Substitute the parameter \a param, if present, on the reset and activation functions. */
	void substitute(const RealParameter& param) {
		this->_activation.substitute(param);
		this->_reset.substitute(param);
	}

	/*! \brief Get the parameters (i.e., the RealConstant whose name starts with a letter) from the transition and reset dynamics. */
	RealParameterSet parameters() const {
		RealParameterSet result = this->_activation.parameters();
		RealParameterSet reset_parameters = this->_reset.parameters();
		result.insert(reset_parameters.begin(),reset_parameters.end());
		return result;
	}

    //! \brief The activation region of the discrete transition.
    const VectorFunction& activation() const {
        return this->_activation;
    }

    //! \brief The activation region of the discrete transition.
    const ScalarFunction scalar_activation() const {
        ARIADNE_ASSERT_MSG(this->_activation.result_size()==1,"Constraint "<<this->_activation<<" is not scalar.");
        return this->_activation[0];
    }

    //! \brief The reset map of the discrete transition.
    const VectorFunction& reset() const {
        return this->_reset;
    }

    EventKind kind() const {
        return this->_kind;
    }

  private:


    // Construct from shared pointers (for internal use).
    DiscreteTransition(DiscreteEvent event,
                       const DiscreteMode& source,
                       const DiscreteMode& target,
                       const VectorFunction& reset,
                       const VectorFunction& activation,
                       EventKind kind);

};

std::ostream& operator<<(std::ostream& os, const DiscreteTransition& dt);


inline bool operator<(const DiscreteTransition& transition1, const DiscreteTransition& transition2) {
    return transition1.event() < transition2.event()
        || (transition1.event() == transition2.event()
            && transition1.source() < transition2.source());
}





/*! \brief A hybrid automaton, comprising continuous-time behaviour
 *  at each discrete mode, coupled by instantaneous discrete transitions.
 *  The state space is given by a hybrid set.
 *
 * A hybrid automaton is a dynamic system with evolution in both
 * continuous time and discrete time.
 * The state space is a product \f$X=\bigcup\{q\}\times X_q\f$
 * where \f$q\f$ is the <em>discrete location</em> and \f$X_q\f$
 * is the <em>continuous state space</em> of corresponding to
 * each discrete location.
 *
 * For each %DiscreteMode, the dynamics is given by a
 * %VectorField describing the continuous dynamics,
 * and a %Set giving an invariants which must be satisified at
 * all times.
 *
 * The discrete time behaviour is specified by %DiscreteTransition
 * objects.
 * Each discrete transition represents an jump from a \a source
 * mode to a \a target mode.
 * There can be at most one discrete transition in an automaton
 * with the same event and source.
 *
 * A discrete transision can either be \em forced or \em unforced.
 * A forced transition much occur as soon as it is activated.
 * An unforced transition may occur at any time it is activated,
 * but is only forced to occur if the continuous evolution is
 * blocked by an invariant.
 *
 * \sa \link Ariadne::DiscreteMode \c DiscreteMode \endlink, \link Ariadne::DiscreteTransition \c DiscreteTransition \endlink

 */
class HybridAutomaton
{
  public:
    //! \brief The type used to represent time.
    typedef HybridTime TimeType;
    //! \brief The type used to represent real numbers.
    typedef double RealType ;
    //! \brief The type used to describe the state space.
    typedef HybridSpace StateSpaceType;


    typedef std::map<DiscreteEvent,VectorFunction>::const_iterator invariant_const_iterator;
    typedef std::list<DiscreteTransition>::const_iterator discrete_transition_const_iterator;
    typedef std::list<DiscreteMode>::const_iterator discrete_mode_const_iterator;
  private:
    //! \brief The hybrid automaton's name.
    std::string _name;

    //! \brief The list of the hybrid automaton's discrete modes.
    std::list< DiscreteMode > _modes;

    //! \brief The hybrid automaton's transitions.
    std::list< DiscreteTransition > _transitions;

  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Construct an empty automaton with no name
    HybridAutomaton();

    //! \brief Construct an empty automaton with the given name
    HybridAutomaton(const std::string& name);

    //! \brief Construct dynamically-allocated copy. (Not currently implemented)
    HybridAutomaton* clone() const;

    //! \brief  Destructor.
    ~HybridAutomaton();
    //@}

    //@{
    //! \name Methods for building the automaton.

    //! \brief Adds a discrete mode to the automaton.
    //!
    //!   \param location is the mode's discrete location.
    //!   \param dynamic is the mode's vector field.
    const DiscreteMode& new_mode(DiscreteLocation location,
                                 const VectorFunction& dynamic);

    //! \brief Adds an invariant to a mode of the automaton.
    //!
    //!   \param location is the mode's discrete location.
    //!   \param invariant is the new invariant condition, in the form \f$g(x)<0\f$.

    const DiscreteMode& new_invariant(DiscreteLocation location,
                                      const ScalarFunction& invariant);

    //! \brief Adds an invariants to a mode of the automaton.
    //!
    //!   \param location is the mode's discrete location.
    //!   \param invariants is the new invariants condition.

    const DiscreteMode& new_invariant(DiscreteLocation location,
                                      const VectorFunction& invariants);

    //! \brief Adds an invariants to a mode of the automaton.
    //!
    //!    \param mode is the discrete mode.
    //!    \param invariants is the new invariants condition.

    const DiscreteMode& new_invariant(const DiscreteMode& mode,
                                      const VectorFunction& invariants);


    //! \brief Adds a discrete transition to the automaton using the discrete locations to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    //!    \param kind determines whether the transision is forced (URGENT) or unforced (PERMISSIVE).
    const DiscreteTransition& new_transition(DiscreteEvent event,
                                             DiscreteLocation source,
                                             DiscreteLocation target,
                                             const VectorFunction& reset,
                                             const ScalarFunction& activation,
                                             EventKind kind);

    //! \brief Adds a discrete transition to the automaton using the discrete locations to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    //!    \param kind determines whether the transision is forced (URGENT) or unforced (PERMISSIVE).
    const DiscreteTransition& new_transition(DiscreteEvent event,
                                             DiscreteLocation source,
                                             DiscreteLocation target,
                                             const VectorFunction& reset,
                                             const VectorFunction& activation,
                                             EventKind kind);

    //! \brief Adds a discrete transition to the automaton using the discrete modes to specify the source and target.
    //!
    //!    \param event is the discrete transition's discrete event.
    //!    \param source is the discrete transition's source mode.
    //!    \param target is the discrete transition's target mode.
    //!    \param reset is the discrete transition's reset.
    //!    \param activation is the discrete transition's activation region.
    //!    \param kind determines whether the transision is forced (URGENT) or unforced (PERMISSIVE).
    const DiscreteTransition& new_transition(DiscreteEvent event,
                                             const DiscreteMode& source,
                                             const DiscreteMode& target,
                                             const VectorFunction& reset,
                                             const ScalarFunction& activation,
                                             EventKind kind);

    //! \brief Adds a discrete transition to the automaton using the discrete modes to specify the source and target.
    //!
    //!    \param event is the discrete transition's discrete event.
    //!    \param source is the discrete transition's source mode.
    //!    \param target is the discrete transition's target mode.
    //!    \param reset is the discrete transition's reset.
    //!    \param activation is the discrete transition's activation region.
    //!    \param kind determines whether the transision is forced (URGENT) or unforced (PERMISSIVE).
    const DiscreteTransition& new_transition(DiscreteEvent event,
                                             const DiscreteMode& source,
                                             const DiscreteMode& target,
                                             const VectorFunction& reset,
                                             const VectorFunction& activation,
                                             EventKind kind);

    //! \brief Adds a forced (urgent) discrete transition to the automaton
    //! using the discrete locations to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const DiscreteTransition& new_forced_transition(DiscreteEvent event,
                                                    DiscreteLocation source,
                                                    DiscreteLocation target,
                                                    const VectorFunction& reset,
                                                    const ScalarFunction& activation);

    //! \brief Adds a forced (urgent) discrete transition to the automaton
    //! using the discrete locations to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const DiscreteTransition& new_forced_transition(DiscreteEvent event,
                                                    DiscreteLocation source,
                                                    DiscreteLocation target,
                                                    const VectorFunction& reset,
                                                    const VectorFunction& activation);

    //! \brief Adds a forced (urgent) discrete transition to the automaton
    //! using the discrete modes to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const DiscreteTransition& new_forced_transition(DiscreteEvent event,
                                                    const DiscreteMode& source,
                                                    const DiscreteMode& target,
                                                    const VectorFunction& reset,
                                                    const ScalarFunction& activation);

    //! \brief Adds a forced (urgent) discrete transition to the automaton
    //! using the discrete modes to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const DiscreteTransition& new_forced_transition(DiscreteEvent event,
                                                    const DiscreteMode& source,
                                                    const DiscreteMode& target,
                                                    const VectorFunction& reset,
                                                    const VectorFunction& activation);

    //! \brief Adds an unforced (non-urgent) discrete transition to the automaton
    //! using the discrete locations to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const DiscreteTransition& new_unforced_transition(DiscreteEvent event,
                                                      DiscreteLocation source,
                                                      DiscreteLocation target,
                                                      const VectorFunction& reset,
                                                      const ScalarFunction& activation);

    //! \brief Adds an unforced (non-urgent) discrete transition to the automaton
    //! using the discrete locations to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const DiscreteTransition& new_unforced_transition(DiscreteEvent event,
                                                      DiscreteLocation source,
                                                      DiscreteLocation target,
                                                      const VectorFunction& reset,
                                                      const VectorFunction& activation);

    //! \brief Adds an unforced (non-urgent) discrete transition to the automaton
    //! using the discrete locations to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const DiscreteTransition& new_unforced_transition(DiscreteEvent event,
                                                      const DiscreteMode& source,
                                                      const DiscreteMode& target,
                                                      const VectorFunction& reset,
                                                      const ScalarFunction& activation);

    //! \brief Adds an unforced (non-urgent) discrete transition to the automaton
    //! using the discrete locations to specify the source and target modes.
    //!
    //!    \param event is the transition's event.
    //!    \param source is the transition's source location.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param activation is the transition's activation region.
    const DiscreteTransition& new_unforced_transition(DiscreteEvent event,
                                                      const DiscreteMode& source,
                                                      const DiscreteMode& target,
                                                      const VectorFunction& reset,
                                                      const VectorFunction& activation);

    //! \brief Gets the parameters (i.e., any RealParameter whose name starts from a letter) from the modes and transitions. */
    RealParameterSet parameters() const;

    //! \brief Get the value of a parameter.
    Real parameter_value(String name) const;

	/*! \brief Substitute the parameter \a param, if present, on all the functions of modes and transitions. */
	void substitute(RealParameter param);

	/*! \brief Substitute parameters from a set \a params. */
	void substitute(const RealParameterSet& params);

	/*! \brief Substitute parameter values from a set \a params, using the midpoint if \a use_midpoint is set. */
	void substitute(const RealParameterSet& params, bool use_midpoint);


	//@}

    //@{
    //! \name Data access and queries.

    //! \brief Returns the hybrid automaton's name.
    const std::string& name() const;

    //! \brief Test if the hybrid automaton has a discrete mode with discrete location \a location.
    bool has_mode(DiscreteLocation location) const;

    //! \brief The kind (permissive, urgent etc) of the event.
    virtual EventKind event_kind(DiscreteLocation location, DiscreteEvent event) const;

    //! \brief Test if the hybrid automaton has a discrete transition with \a event_id and \a source_id.
    bool has_transition(DiscreteEvent event, DiscreteLocation source) const;

    //! \brief The discrete mode with given discrete location.
    const DiscreteMode& mode(DiscreteLocation location) const;

    //! \brief The discrete transition with given \a event and \a source location.
    const DiscreteTransition& transition(DiscreteEvent event, DiscreteLocation source) const;

    //! \brief The set of discrete modes. 
    const std::list< DiscreteMode >& modes() const;

    //! \brief The set of discrete transitions. 
    const std::list< DiscreteTransition >& transitions() const;

    //! \brief The discrete transitions from location \a source.
    std::list< DiscreteTransition > transitions(DiscreteLocation source) const;

    //! \brief The blocking events (invariants and urgent transitions) in \a location.
    std::map<DiscreteEvent,VectorFunction> blocking_guards(DiscreteLocation location) const;

    //! \brief The permissive events (invariants and urgent transitions) in \a location.
    std::map<DiscreteEvent,VectorFunction> permissive_guards(DiscreteLocation location) const;

    //! \brief The state space of the system.
    HybridSpace state_space() const;

    //@}
    
    //! \brief non-standard assignment operator from a pair (HybridAutomaton, RealSpace).
    HybridAutomaton& operator=(const std::pair< HybridAutomaton, RealSpace >& pair);

};

std::ostream& operator<<(std::ostream& os, const HybridAutomaton& ha);

/*! \brief Returns the subset of the parameters whose values are not singletons. */
RealParameterSet nonsingleton_parameters(const RealParameterSet& parameters);

/*! \brief Returns the set of identifiers of the \a parameters. */
Set<Identifier> parameters_identifiers(const RealParameterSet& parameters);


} // namespace Ariadne

#endif // ARIADNE_HYBRID_AUTOMATON_H
