/***************************************************************************
 *            hybrid_automaton.cc
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

#include <map>

#include "macros.h"
#include "stlio.h"
#include "function.h"
#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_automaton.h"
#include "discrete_event.h"
#include "space.h"


namespace Ariadne {

typedef uint DimensionType;

class HybridSet {};


uint
DiscreteMode::
dimension() const
{
    return this->_dynamic.argument_size();
}


DiscreteMode::
DiscreteMode(DiscreteLocation location,
             const VectorFunction& dynamic)
    :  _location(location), _dynamic(dynamic), _invariants()
{
}


/*
DiscreteMode::
DiscreteMode(DiscreteLocation location,
             const boost::shared_ptr< const VectorFunction > dynamic,
             const std::map< DiscreteEvent, boost::shared_ptr< const VectorFunction > >& invariants)
    :  _location(location), _dynamic(dynamic), _invariants(invariants), _grid(new Grid(dynamic->argument_size()))
{
    ARIADNE_ASSERT(dynamic->result_size()==dynamic->argument_size());
    for(uint i=0; i!=invariants.size(); ++i) {
        ARIADNE_ASSERT(invariants[i]->argument_size()==dynamic->argument_size());
        ARIADNE_ASSERT(invariants[i]->result_size()==1u);
    }
}
*/


std::ostream&
operator<<(std::ostream& os, const DiscreteMode& mode)
{
    return os << "DiscreteMode( "
              << "location=" << mode.location() << ", "
              << "dynamic=" << mode.dynamic() << ", "
              << "invariants=" << mode.invariants() << " )";
}





DiscreteTransition::
DiscreteTransition(DiscreteEvent event,
                   const DiscreteMode& source,
                   const DiscreteMode& target,
                   const VectorFunction& reset,
                   const VectorFunction& activation,
                   EventKind kind)
    : _event(event), _source(source.location()), _target(target.location()),
      _activation(activation), _reset(reset), _kind(kind)
{
    ARIADNE_ASSERT(activation.result_size()==1);
    ARIADNE_ASSERT(activation.argument_size()==source.dimension());
    ARIADNE_ASSERT(reset.argument_size()==source.dimension());
    ARIADNE_ASSERT(reset.result_size()==target.dimension());
}




std::ostream&
operator<<(std::ostream& os, const DiscreteTransition& transition)
{
    return os << "DiscreteTransition( "
              << "event=" << transition.event() << ", "
              << "source=" << transition.source() << ", "
              << "target=" << transition.target() << ", "
              << "reset=" << transition.reset() << ", "
              << "activation=" << transition.activation() << ", "
              << "kind=" << transition.kind() << " )";
}




HybridAutomaton::~HybridAutomaton()
{
}

HybridAutomaton::HybridAutomaton()
{
}

HybridAutomaton::HybridAutomaton(const std::string& name)
    : _name(name)
{
}





const DiscreteMode&
HybridAutomaton::new_mode(DiscreteLocation location,
                          const VectorFunction& dynamic)
{
    if(this->has_mode(location)) {
        throw std::runtime_error("The hybrid automaton already has a mode with the given id");
    }
    if(dynamic.result_size()!=dynamic.argument_size()) {
        ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_mode(location,dynamic)",
            "The dynamic has argument size " << dynamic.argument_size()
                << " and result size " << dynamic.result_size() << ", so does not define a vector field.");
    }
    this->_modes.push_back(DiscreteMode(location,dynamic));
    return this->mode(location);
}


const DiscreteMode&
HybridAutomaton::new_invariant(DiscreteLocation location,
                               const ScalarFunction& invariant)
{
    return this->new_invariant(location, VectorFunction(1u,invariant));
}


const DiscreteMode&
HybridAutomaton::new_invariant(DiscreteLocation location,
                               const VectorFunction& invariant)
{
    if(!this->has_mode(location)) {
        throw std::runtime_error("The location of the invariant must be in the automaton.");
    }
    DiscreteMode& mode=const_cast<DiscreteMode&>(this->mode(location));
	DiscreteEvent invariant_event(location.name()+"-inv"+to_str(mode._invariants.size()));

    if(this->has_transition(location,invariant_event)) {
        throw std::runtime_error("The automaton already has a transition with the would-be id and the given location as source id.");
    }

    if(invariant.argument_size()!=mode.dimension()) {
        ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_invariant(location,invariant)",
            "The invariant has argument size " << invariant.argument_size()
                << " but the mode has state-space dimension " << mode.dimension());
    }
    if(invariant.result_size()!=1u) {
        ARIADNE_THROW(std::runtime_error,"HybridAutomaton::new_invariant(location,invariant)",
            "The invariant has result size " << invariant.result_size()
                << " but only scalar invariants are currently supported.");
    }
    mode._invariants[invariant_event]=invariant;
    return mode;
}

const DiscreteTransition&
HybridAutomaton::new_transition(DiscreteEvent event,
                                DiscreteLocation source,
                                DiscreteLocation target,
                                const VectorFunction &reset,
                                const VectorFunction &activation,
                                EventKind kind)
{
    if(this->has_transition(source,event)) {
        throw std::runtime_error("The automaton already has a transition with the given id and source id.");
    }
    if(this->has_invariant(source,event)) {
    	throw std::runtime_error("The automaton already has an invariant with the given id and location id.");
    }
    if(!this->has_mode(source)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given source id");
    }
    if(!this->has_mode(target)) {
        throw std::runtime_error("The automaton does not contain a mode with ths given desitination id");
    }

    const DiscreteMode& source_mode=this->mode(source);
    const DiscreteMode& target_mode=this->mode(target);
    this->_transitions.push_back(DiscreteTransition(event,source_mode,target_mode,reset,activation,kind));
    return this->transition(event,source);
}

const DiscreteTransition&
HybridAutomaton::
new_transition(DiscreteEvent event,
               DiscreteLocation source,
               DiscreteLocation target,
               const VectorFunction &reset,
               const ScalarFunction &activation,
               EventKind kind)
{
    return this->new_transition(event,source,target,reset,VectorFunction(1u,activation),kind);
}

const DiscreteTransition&
HybridAutomaton::
new_transition(DiscreteEvent event,
               const DiscreteMode &source,
               const DiscreteMode &target,
               const VectorFunction &reset,
               const VectorFunction &activation,
               EventKind kind)
{
    DiscreteLocation source_id=source.location();
    DiscreteLocation target_id=target.location();
    return this->new_transition(event,source_id,target_id,reset,activation,kind);
}

const DiscreteTransition&
HybridAutomaton::
new_transition(DiscreteEvent event,
               const DiscreteMode &source,
               const DiscreteMode &target,
               const VectorFunction &reset,
               const ScalarFunction &activation,
               EventKind kind)
{
    DiscreteLocation source_id=source.location();
    DiscreteLocation target_id=target.location();
    return this->new_transition(event,source_id,target_id,reset,VectorFunction(1u,activation),kind);
}

const DiscreteTransition&
HybridAutomaton::
new_forced_transition(DiscreteEvent event,
                      DiscreteLocation source,
                      DiscreteLocation target,
                      const VectorFunction &reset,
                      const ScalarFunction &activation)
{
    return this->new_transition(event,source,target,reset,VectorFunction(1u,activation),URGENT);
}

const DiscreteTransition&
HybridAutomaton::
new_forced_transition(DiscreteEvent event,
                      DiscreteLocation source,
                      DiscreteLocation target,
                      const VectorFunction &reset,
                      const VectorFunction &activation)
{
    return this->new_transition(event,source,target,reset,activation,URGENT);
}

const DiscreteTransition&
HybridAutomaton::
new_forced_transition(DiscreteEvent event,
                      const DiscreteMode &source,
                      const DiscreteMode &target,
                      const VectorFunction &reset,
                      const VectorFunction &activation)
{
    DiscreteLocation source_id=source.location();
    DiscreteLocation target_id=target.location();
    return this->new_transition(event,source_id,target_id,reset,activation,URGENT);
}

const DiscreteTransition&
HybridAutomaton::
new_forced_transition(DiscreteEvent event,
                      const DiscreteMode &source,
                      const DiscreteMode &target,
                      const VectorFunction &reset,
                      const ScalarFunction &activation)
{
    DiscreteLocation source_id=source.location();
    DiscreteLocation target_id=target.location();
    return this->new_transition(event,source_id,target_id,reset,VectorFunction(1u,activation),URGENT);
}


const DiscreteTransition&
HybridAutomaton::
new_unforced_transition(DiscreteEvent event,
                        DiscreteLocation source,
                        DiscreteLocation target,
                        const VectorFunction &reset,
                        const ScalarFunction &activation)
{
    return this->new_transition(event,source,target,reset,VectorFunction(1u,activation),PERMISSIVE);
}

const DiscreteTransition&
HybridAutomaton::
new_unforced_transition(DiscreteEvent event,
                        DiscreteLocation source,
                        DiscreteLocation target,
                        const VectorFunction &reset,
                        const VectorFunction &activation)
{
    return this->new_transition(event,source,target,reset,activation,PERMISSIVE);

}

const DiscreteTransition&
HybridAutomaton::
new_unforced_transition(DiscreteEvent event,
                        const DiscreteMode &source,
                        const DiscreteMode &target,
                        const VectorFunction &reset,
                        const VectorFunction &activation)
{
    DiscreteLocation source_id=source.location();
    DiscreteLocation target_id=target.location();
    return this->new_transition(event,source_id,target_id,reset,activation,PERMISSIVE);
}

const DiscreteTransition&
HybridAutomaton::
new_unforced_transition(DiscreteEvent event,
                        const DiscreteMode &source,
                        const DiscreteMode &target,
                        const VectorFunction &reset,
                        const ScalarFunction &activation)
{
    DiscreteLocation source_id=source.location();
    DiscreteLocation target_id=target.location();
    return this->new_transition(event,source_id,target_id,reset,VectorFunction(1u,activation),PERMISSIVE);
}


EventKind
HybridAutomaton::event_kind(DiscreteLocation location, DiscreteEvent event) const
{
	const DiscreteMode& mode = this->mode(location);

	std::list< DiscreteTransition > trans = transitions(location);

	for (std::list<DiscreteTransition>::const_iterator trans_it = trans.begin(); trans_it != trans.end(); ++trans_it) {
		if (trans_it->event() == event)
			return trans_it->kind();
	}

	for (std::map<DiscreteEvent,VectorFunction>::const_iterator inv_it = mode.invariants().begin();
			inv_it != mode.invariants().end(); ++inv_it) {
		if (inv_it->first == event)
			return INVARIANT;
	}

	ARIADNE_FAIL_MSG("The event '" << event.name() << "' for location '" << location.name() << "' is not present.");
}


bool
HybridAutomaton::has_mode(DiscreteLocation location) const
{
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter) {
        if(mode_iter->location()==location)
            return true;

    }
    return false;
}


Set<DiscreteEvent>
HybridAutomaton::events(DiscreteLocation location) const
{
	Set<DiscreteEvent> result;

	const DiscreteMode& mode = this->mode(location);

	std::list< DiscreteTransition > trans = this->transitions(location);

	for (std::list<DiscreteTransition>::const_iterator trans_it = trans.begin(); trans_it != trans.end(); ++trans_it) {
		result.insert(trans_it->event());
	}

	for (std::map<DiscreteEvent,VectorFunction>::const_iterator inv_it = mode.invariants().begin();
			inv_it != mode.invariants().end(); ++inv_it) {
		result.insert(inv_it->first);
	}

	return result;
}


bool
HybridAutomaton::has_guard(DiscreteLocation location, DiscreteEvent event) const
{
	return this->has_transition(location,event) || this->has_invariant(location,event);
}


bool
HybridAutomaton::has_transition(DiscreteLocation source, DiscreteEvent event) const
{
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->event()==event && transition_iter->source()==source) {
                return true;
            }
        }
    return false;
}


bool
HybridAutomaton::has_invariant(DiscreteLocation location, DiscreteEvent event) const
{
	DiscreteMode mode = this->mode(location);
    for(std::map<DiscreteEvent,VectorFunction>::const_iterator inv_it=mode.invariants().begin();
        inv_it!=mode.invariants().end(); ++inv_it)
        {
            if(inv_it->first==event)
                return true;
        }
    return false;
}


DiscreteLocation
HybridAutomaton::target(DiscreteLocation source, DiscreteEvent event) const {
    if(this->has_transition(source,event)) {
        return this->transition(event,source).target();
    } else {
        return source;
    }
}


uint
HybridAutomaton::dimension(DiscreteLocation location) const
{
	uint dim = 0;
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter) {
    	if (mode_iter->location() == location) {
    		dim = mode_iter->dimension();
    		break;
    	}
    }
    ARIADNE_ASSERT_MSG(dim > 0, "The location is not present into the automaton.");
    return dim;
}


HybridSpace
HybridAutomaton::state_space() const
{
    HybridSpace result;
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            result[mode_iter->location()]=mode_iter->dimension();
        }
    return result;
}


RealSpace
HybridAutomaton::continuous_state_space(DiscreteLocation location) const
{
	ARIADNE_NOT_IMPLEMENTED;
}


RealVectorFunction
HybridAutomaton::dynamic_function(DiscreteLocation location) const
{
	RealVectorFunction func;

	bool found = false;
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter) {
    	if (mode_iter->location() == location) {
    		func = mode_iter->_dynamic;
    		found = true;
    		break;
    	}
    }
    ARIADNE_ASSERT_MSG(found, "The location is not present into the automaton.");
    return func;
}


RealScalarFunction
HybridAutomaton::invariant_function(DiscreteLocation location, DiscreteEvent event) const
{
	RealScalarFunction func;

	bool found = false;
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter) {
    	if (mode_iter->location() == location) {
    		invariant_const_iterator inv_it = mode_iter->_invariants.find(event);
    		ARIADNE_ASSERT_MSG(inv_it != mode_iter->_invariants.end(),
    				"The invariant with event '" << event.name() << "' is not present into the automaton.");
    		func = inv_it->second[0];
    		found = true;
    		break;
    	}
    }
    ARIADNE_ASSERT_MSG(found, "The location is not present into the automaton.");
    return func;
}


RealScalarFunction
HybridAutomaton::guard_function(DiscreteLocation location, DiscreteEvent event) const
{
	RealScalarFunction func;

	bool found = false;
    for(discrete_transition_const_iterator trans_iter=this->_transitions.begin();
    		trans_iter!=this->_transitions.end(); ++trans_iter) {
    	if (trans_iter->source() == location && trans_iter->event() == event) {
    		func = trans_iter->_activation[0];
    		found = true;
    		break;
    	}
    }
    ARIADNE_ASSERT_MSG(found, "No transition with the given location and event is present into the automaton.");
    return func;
}


RealVectorFunction
HybridAutomaton::reset_function(DiscreteLocation location, DiscreteEvent event) const
{
	RealVectorFunction func;

	bool found = false;
    for(discrete_transition_const_iterator trans_iter=this->_transitions.begin();
    		trans_iter!=this->_transitions.end(); ++trans_iter) {
    	if (trans_iter->source() == location && trans_iter->event() == event) {
    		func = trans_iter->reset();
    		found = true;
    		break;
    	}
    }
    ARIADNE_ASSERT_MSG(found, "No transition with the given location and event is present into the automaton.");
    return func;
}


const std::list< DiscreteMode >&
HybridAutomaton::modes() const
{
    return this->_modes;
}



const DiscreteMode&
HybridAutomaton::mode(DiscreteLocation location) const
{
    // FIXME: This is a hack; we should use a logarithmic time real search to find a mode with the given discrete state.
    for(discrete_mode_const_iterator mode_iter=this->_modes.begin();
        mode_iter!=this->_modes.end(); ++mode_iter)
        {
            if(mode_iter->location()==location) {
                return *mode_iter;
            }
        }
    throw std::runtime_error("The hybrid automaton does not have a mode with the given id.");
}


const std::list< DiscreteTransition >&
HybridAutomaton::transitions() const
{
    return this->_transitions;
}


std::list< DiscreteTransition >
HybridAutomaton::transitions(DiscreteLocation source) const
{
    std::list< DiscreteTransition > result;
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->source()==source) {
                result.push_back(*transition_iter);
            }
        }
    return result;
}


std::map<DiscreteEvent,VectorFunction>
HybridAutomaton::blocking_guards(DiscreteLocation source) const
{
    std::map<DiscreteEvent,VectorFunction> result;
    const DiscreteMode& mode=this->mode(source);
    for(invariant_const_iterator invariant_iter=mode._invariants.begin();
        invariant_iter!=mode._invariants.end(); ++invariant_iter)
    {
        const DiscreteEvent event=invariant_iter->first;
        const VectorFunction invariant=invariant_iter->second;
        result[event]=invariant;
    }

    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->source()==source && transition_iter->kind() == URGENT) {
            const DiscreteEvent event=transition_iter->event();
            const VectorFunction guard=transition_iter->activation();
            result[event]=guard;
        }
    }
    return result;
}


std::map<DiscreteEvent,VectorFunction>
HybridAutomaton::permissive_guards(DiscreteLocation source) const
{
    std::map<DiscreteEvent,VectorFunction> result;

    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
    {
        if(transition_iter->source()==source && transition_iter->kind() == PERMISSIVE) {
            const DiscreteEvent event=transition_iter->event();
            const VectorFunction guard=transition_iter->activation();
            result[event]=guard;
        }
    }
    return result;
}



const DiscreteTransition&
HybridAutomaton::transition(DiscreteEvent event, DiscreteLocation source) const
{
    for(discrete_transition_const_iterator transition_iter=this->_transitions.begin();
        transition_iter!=this->_transitions.end(); ++transition_iter)
        {
            if(transition_iter->event()==event && transition_iter->source()==source) {
                return *transition_iter;
            }
        }
    ARIADNE_THROW(std::runtime_error, "HybridAutomaton::transition(event, source)",
        "The hybrid automaton does not have a transition with event \"" << event << "\" and source \"" << source << "\".");
}


const String&
HybridAutomaton::name() const
{
    return this->_name;
}

HybridAutomaton& 
HybridAutomaton::operator=(const std::pair< HybridAutomaton, RealSpace >& pair) 
{
    //std::cout << pair.first << std::endl;
    // Handle self-assignment
    if(this == &(pair.first)) return *this;
    this->_name = pair.first.name();
    this->_modes = pair.first.modes();
    this->_transitions = pair.first.transitions();
    return *this;
}


std::ostream&
operator<<(std::ostream& os, const HybridAutomaton& ha)
{
    return os << "HybridAutomaton( name=" << ha.name() <<
    							", parameters" << ha.parameters() <<
    							", modes=" << ha.modes() <<
    							", transitions=" << ha.transitions() <<
    							")";
}

void 
HybridAutomaton::substitute(const RealParameter& param)
{
	RealParameterSet parameters = this->parameters();

	bool found = false;
	for (RealParameterSet::const_iterator param_it = parameters.begin(); param_it != parameters.end(); ++param_it)
		if (param_it->name() == param.name()) {
			found = true;
			break;
		}

	ARIADNE_ASSERT_MSG(found, "The parameter to substitute is not present in the system.");

	for (std::list<DiscreteMode>::iterator modes_it=this->_modes.begin();modes_it!=this->_modes.end();modes_it++)
		modes_it->substitute(param);
	for (std::list<DiscreteTransition>::iterator trans_it=this->_transitions.begin();trans_it!=this->_transitions.end();trans_it++)
		trans_it->substitute(param);
}


void
HybridAutomaton::substitute(const RealParameterSet& params, bool use_midpoint)
{
	for (RealParameterSet::const_iterator param_it = params.begin(); param_it != params.end(); ++param_it) {
		if (use_midpoint)
			substitute(RealParameter(param_it->name(),param_it->value().midpoint()));
		else
			substitute(*param_it);
	}
}


RealParameterSet
HybridAutomaton::parameters() const
{
	RealParameterSet result;

	for (std::list<DiscreteMode>::const_iterator modes_it=this->_modes.begin();modes_it!=this->_modes.end();modes_it++) {
		RealParameterSet mode_result = modes_it->parameters();
		result.insert(mode_result.begin(),mode_result.end());
	}
	for (std::list<DiscreteTransition>::const_iterator trans_it=this->_transitions.begin();trans_it!=this->_transitions.end();trans_it++) {
		RealParameterSet trans_result = trans_it->parameters();
		result.insert(trans_result.begin(),trans_result.end());
	}
	return result;
}


RealParameterSet
nonsingleton_parameters(const RealParameterSet& parameters)
{
	RealParameterSet result;
	for (RealParameterSet::const_iterator param_it = parameters.begin(); param_it != parameters.end(); ++param_it) {
		if (!param_it->value().singleton())
			result.insert(*param_it);
	}

	return result;
}

Set<Identifier>
parameters_identifiers(const RealParameterSet& parameters)
{
	Set<Identifier> result;

	for (RealParameterSet::const_iterator param_it = parameters.begin(); param_it != parameters.end(); ++param_it) {
		result.insert(param_it->name());
	}

	return result;
}


}
