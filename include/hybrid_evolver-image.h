/***************************************************************************
 *            hybrid_evolver-image.h
 *
 *  Copyright  2007-10  Alberto Casagrande, Pieter Collins, Luca Geretti
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

/*! \file hybrid_evolver-image.h
 *  \brief Evolver for hybrid systems using image sets.
 */

#ifndef ARIADNE_HYBRID_EVOLVER_IMAGE_H
#define ARIADNE_HYBRID_EVOLVER_IMAGE_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "tuple.h"

#include "hybrid_set.h"

#include "hybrid_automaton.h"
#include "evolver_interface.h"
#include "evolver_base.h"

namespace Ariadne {

template<class Sys, class BS> class Evolver;

class VectorTaylorFunction;
class TaylorSet;
class HybridAutomatonInterface;
template<class ES> class Orbit;
class EvolutionSettings;
class TaylorModel;
template<class MDL> class CalculusInterface;
class HybridTime;
class DiscreteEvent;
class ImageSetHybridEvolverSettings;
class EvolutionData;
class ContinuousStepResult;

/*! \brief A class for computing the evolution of a hybrid system.
 *
 * The actual evolution steps are performed by the HybridEvolver class.
 */
class ImageSetHybridEvolver
    : public EvolverBase<HybridAutomatonInterface,LocalisedTaylorSet>
{
    friend class EvolutionData;
    friend class ContinuousStepResult;

    typedef ScalarFunction ScalarFunctionType;
    typedef VectorFunction VectorFunctionType;
    typedef Vector<Interval> BoxType;
    typedef VectorTaylorFunction FunctionModelType;
    typedef VectorTaylorFunction MapModelType;
    typedef VectorTaylorFunction FlowModelType;
    typedef ScalarTaylorFunction ConstraintModelType;
    typedef TaylorModel TimeModelType;
    typedef TaylorSet FlowSetModelType;
    typedef TaylorSet SetModelType;
  public:
    typedef ImageSetHybridEvolverSettings SettingsType;
    typedef HybridAutomatonInterface::TimeType TimeType;
    typedef int IntegerType;
    typedef uint AccuracyType;
    typedef Float RealType;
    typedef std::vector<DiscreteEvent> EventListType;
    typedef HybridAutomatonInterface SystemType;
    typedef TaylorSet ContinuousEnclosureType;
    typedef TaylorSet TimedSetModelType;
    typedef HybridBasicSet<TaylorSet> LocalisedEnclosureType;
    typedef LocalisedEnclosureType EnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<EnclosureType> EnclosureListType;
    typedef Float ContinuousTimeType;
    typedef Float StepType;
    typedef std::map< DiscreteEvent,tuple<TaylorModel,TaylorModel> > ActivationTimesType;
  public:

    //! \brief Construct from a system.
    ImageSetHybridEvolver(const SystemType& system);

    /*! \brief Make a dynamically-allocated copy. */
    ImageSetHybridEvolver* clone() const { return new ImageSetHybridEvolver(*this); }

    //@{
    //! \name Settings controlling the evolution.

    //! \brief A reference to the settings controlling the evolution.
    SettingsType& settings() { return *this->_settings; }
	//! \brief A constant reference to the settings controlling the evolution.
    const SettingsType& settings() const { return *this->_settings; }

	//! \brief Modifies the settings, given some metrics.
	virtual void tune_settings(
			const HybridBoxes& domain,
			const HybridGrid& grid,
			AccuracyType accuracy);

    //@}

  protected:
    virtual void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                            const EnclosureType& initial, const TimeType& time, bool ignore_activations,
                            ContinuousEvolutionDirection direction, Semantics semantics) const;

    virtual void _evolution_step(std::list<EvolutionData>& working_sets,
                                  EnclosureListType& reachable, EnclosureListType& intermediate,
                                  const EvolutionData& data, const TimeType& time,
                                  bool ignore_activations, ContinuousEvolutionDirection direction, Semantics semantics) const;

  protected:
    TimeModelType crossing_time(VectorFunction guard, const FlowSetModelType& flow_set, Semantics semantics) const;

    Interval normal_derivative(VectorFunction guard, const FlowSetModelType& flow_set, const TimeModelType& crossing_time) const;

    /*! \brief Computes the initially active events
     * \details Produces one entry for each possibly initially active event. Adds a convenience summary entry for the generic blocking_event,
     * which has activity FALSE iff no other entry exists, INDETERMINATE if other entries exist but no definitely active event exists, TRUE otherwise.
     */
    void compute_initially_active_events(
    		std::map<DiscreteEvent,tribool>&,
            const std::map<DiscreteEvent,RealScalarFunction>&,
            const ContinuousEnclosureType&,
            Semantics semantics) const;

    bool has_nonnegative_crossing(
    		const std::map<DiscreteEvent,RealScalarFunction>& urgent_guards,
			const RealVectorFunction dynamic,
			const Box& set_bounds,
			Semantics semantics) const;

    bool is_enclosure_to_be_discarded(
    		const ContinuousEnclosureType& enclosure,
            const std::map<DiscreteEvent,RealScalarFunction>& urgent_guards,
            const RealVectorFunction& dynamic,
            Semantics semantics) const;

    void compute_flow_model(
    		const DiscreteLocation&,
    		FlowSetModelType&,
    		BoxType&,
    		Float&,
    		RealVectorFunction,
            const SetModelType&,
            const TimeModelType&,
            Float) const;

    FlowSetModelType compute_flow_model_simplified(
            Float&,
            RealVectorFunction,
            Float&,
            const SetModelType&) const;

    std::pair<FlowSetModelType,Float>
    compute_flow_and_effective_step(
            Float& time_step,
            RealVectorFunction dynamic,
            Float& maximum_bounds_diameter,
            const SetModelType& starting_set_model) const;

    ContinuousStepResult
    compute_integration_step_result(const SetModelType& starting_set,
                   const DiscreteLocation& location,
                   ContinuousEvolutionDirection direction,
                   Float step) const;

    void compute_eventBlockingTimes_and_nonTransverseEvents(
    		std::map<DiscreteEvent,TimeModelType>&,
    		std::set<DiscreteEvent>&,
            const std::map<DiscreteEvent,RealScalarFunction>& urgent_guards,
            const FlowSetModelType& flow_set_model,
            Semantics semantics) const;

    void compute_blockingTime_and_relatedEvents(
    		std::set<DiscreteEvent>&, TimeModelType&,
            const std::map<DiscreteEvent,TimeModelType>&) const;

    void compute_activationTimes(
    		std::map<DiscreteEvent,tuple<TimeModelType,TimeModelType> >& activation_times,
            const std::map<DiscreteEvent,RealScalarFunction>& activations,
            const FlowSetModelType& flow_set_model,
            const TimeModelType& blocking_time_model,
            const Semantics semantics) const;

  private:

    std::map<uint,Vector<Float> > _indexed_set_models_widths(std::list<EvolutionData>& working_sets) const;

    void _box_on_contraction(TaylorSet& set_model,
        					 TaylorModel& time_model,
        					 const DiscreteLocation& loc,
        					 const TimeType& maximum_hybrid_time,
        					 ContinuousEvolutionDirection direction,
        					 Semantics semantics) const;

    ContinuousStepResult _adaptive_step_and_flow(const SetModelType& set_model,
                                                 const Float& previous_step,
                                                 const uint& refinement_width,
                                                 const DiscreteLocation& location,
                                                 ContinuousEvolutionDirection direction,
                                                 const Float& remaining_time,
                                                 const Float& maximum_step,
                                                 const SetModelType& maximum_flow_model,
                                                 const SetModelType& maximum_finishing_model
    				 	                         ) const;

    bool _is_enclosure_too_large(
    		const DiscreteLocation& loc,
    		const SetModelType& set_model,
    		const Vector<Float>& initial_set_model_widths) const;

    void _evolution_add_initialSet(std::list<EvolutionData>& working_sets,
    							   const EnclosureType& initial_set,
    							   Semantics semantics) const;

    void _add_models_subdivisions_autoselect(std::list<EvolutionData>& working_sets,
    										 const uint& set_index,
    		  	  	  	  	  	  	  		 const SetModelType& set_model,
    		  	  	  	  	  	  	  		 const TimeModelType& time_model,
    		  	  	  	  	  	  	  		 const DiscreteLocation& location,
    		  	                             const Float& previous_step,
    		  	  	  	  	  	  	  		 const EventListType& previous_events,
    		  	  	  	  	  	  	  		 Semantics semantics) const;

    void _add_models_subdivisions_time(std::list<EvolutionData>& working_sets,
			 	 	 	 	 	 	   const uint& set_index,
    		  	  	  	  	  	  	   const SetModelType& set_model,
    		  	  	  	  	  	  	   const TimeModelType& time_model,
    		  	  	  	  	  	  	   const DiscreteLocation& location,
    		  	                       const Float& previous_step,
    		  	  	  	  	  	  	   const EventListType& previous_events,
    		  	  	  	  	  	  	   Semantics semantics) const;

    void _add_subdivisions(std::list<EvolutionData>& working_sets,
    					   const array< TimedSetModelType >& subdivisions,
    					   const uint& set_index,
    					   const DiscreteLocation& location,
    					   const Float& previous_step,
    					   const EventListType& previous_events,
    					   const uint dimension) const;

    void _log_step_summary(const std::list<EvolutionData>& working_sets,
    					 const EnclosureListType& reach_sets,
    					 const EventListType& events,
    					 const TimeModelType& time_model,
    					 const SetModelType& set_model,
    					 const DiscreteLocation& location,
    					 const Float& step_size) const;

    void _computeEvolutionForEvents(std::list<EvolutionData>& working_sets,
			   	   	   	   	   	    EnclosureListType& intermediate_sets,
			   	   	   	   	   	    const uint& set_index,
			   	   	   	   	   	    const DiscreteLocation& location,
			   	   	   	   	   	    const std::set<DiscreteEvent>& blocking_events,
			   	   	   	   	   	    const EventListType& events,
			   	   	   	   	   	    const ActivationTimesType& activation_times,
			   	   	   	   	   	    const SetModelType& flow_set_model,
			   	   	   	   	   	    const TimeModelType& time_model,
			   	   	   	   	   	    const TimeModelType& blocking_time_model,
			   	   	   	   	   	    const Float& time_step,
			   	   	   	   	   	    const uint& next_refinement_width,
			   	   	   	   	   	    bool ignore_activations,
			   	   	   	   	   	    Semantics semantics) const;

    void _compute_blocking_info(std::set<DiscreteEvent>& non_transverse_events,
    				  	   std::set<DiscreteEvent>& blocking_events,
    				  	   TimeModelType& blocking_time_model,
    				  	   const TimeModelType& time_step_model,
    				  	   const SetModelType& flow_set_model,
    				  	   const std::map<DiscreteEvent,RealScalarFunction>& blocking_guards,
    				  	   double SMALL_RELATIVE_TIME,
    				  	   Semantics semantics) const;

    void _compute_activation_info(std::map<DiscreteEvent,RealScalarFunction>& activations,
    						 	  ActivationTimesType& activation_times,
    						 	  const std::set<DiscreteEvent>& non_transverse_events,
    						 	  const SetModelType& flow_set_model,
    						 	  const TimeModelType& blocking_time_model,
    						 	  const std::map<DiscreteEvent,RealScalarFunction>& blocking_guards,
    						 	  const std::map<DiscreteEvent,RealScalarFunction>& invariants,
    						 	  const Semantics semantics) const;

    void _compute_and_adjoin_reachableSet(EnclosureListType& reach_sets,
    									 SetModelType& reachable_set,
    									 const DiscreteLocation& location,
    									 const SetModelType& flow_set_model,
    									 const TimeModelType& zero_time_model,
    									 const TimeModelType& blocking_time_model,
    									 Semantics semantics) const;

    void _logEvolutionStepInitialState(const EventListType& previous_events,
    							  	   const TimeModelType& time_model,
    							  	   const DiscreteLocation& location,
    							  	   const SetModelType& set_model,
    							  	   const RealVectorFunction& dynamic,
    							  	   const std::map<DiscreteEvent,RealScalarFunction>& invariants,
    							  	   const std::map<DiscreteEvent,RealScalarFunction>& blocking_guards,
    							  	   const std::map<DiscreteEvent,RealScalarFunction>& permissive_guards) const;

  protected:
    // Special events
    static const DiscreteEvent starting_event;
    static const DiscreteEvent finishing_event;
    static const DiscreteEvent blocking_event;

  private:
    boost::shared_ptr< SettingsType > _settings;
    boost::shared_ptr< CalculusInterface<TaylorModel> > _toolbox;
};


//! \brief Settings for controlling the accuracy of evolution methods on enclosure sets.
class ImageSetHybridEvolverSettings {
    friend class ImageSetHybridEvolver;
  public:
    typedef uint UnsignedIntType;
    typedef HybridAutomatonInterface SystemType;

  protected:

    //! \brief Default constructor gives reasonable values.
    ImageSetHybridEvolverSettings(const SystemType& sys);

  private:

    const SystemType& _sys;

    //! \brief A fixed value for the maximum allowable step size for integration, different for each location.
    //! \details Decreasing the values increases the accuracy of the computation. The actual step size may be lower,
    //! if higher than a fixed threshold value related to the dynamics of the system.
    std::map<DiscreteLocation,Float> _maximum_step_size;

    //! \brief The reference enclosure widths that a discretised enclosure would have.
    //! \details If an enclosure starts evolution with widths strictly lesser than these, premature termination is performed
    //! when the widths are _maximum_enclosure_widths_ratio times the _reference_enclosure_widths. If not,
    //! premature termination is performed when the widths are _maximum_enclosure_widths_ratio times the initial enclosure widths.
    HybridFloatVector _reference_enclosure_widths;

    //! \brief The maximum ratio between enclosure widths and reference widths.
    //! \details Reference widths are the initial ones, if larger than _reference_enclosure_widths, or
    //! _reference_enclosure_widths themselves if not.
    Float _maximum_enclosure_widths_ratio;

    //! \brief Enforce the error rate on each step to be bound.
    //! \details Uses the _reference_enclosure_widths to identify a desired final error and, based on the evolution time, tries to
    //! keep the error-per-unit-time below the error/evolution_time ratio.
    bool _enable_error_rate_enforcement;

    //! \brief Boxes the enclosure if the dynamics is contractive (true by default).
    //! \details In order to reduce the error, it overapproximates the enclosure with a box. In such a way, the error (which is
    //! now affected by the dynamics) would be reduced if the set is shrinked. This operation is hence performed only if the local dynamics is contractive.
    bool _enable_boxing_on_contraction;

    //! \brief Recondition sets to handle the increase in the error terms.
    bool _enable_reconditioning;

    //! \brief Enable subdivision of basic sets (false by default).
    bool _enable_subdivisions;

    //! \brief Terminate evolution if basic sets became too large (true by default).
    //! \details In the case of upper semantics, if true and no subdivisions are present, the set is put into the final sets. In the case of lower semantics, the set is discarded.
    bool _enable_premature_termination_on_enclosure_size;

    //! \brief Terminate evolution if too many working sets are present (0 by default, hence disabled).
    unsigned int _maximum_number_of_working_sets;

  public:

    // Accessors

    const std::map<DiscreteLocation,Float>& maximum_step_size() const;
    void set_maximum_step_size(const Float&);
    void set_maximum_step_size(const std::map<DiscreteLocation,Float>&);

    const HybridFloatVector& reference_enclosure_widths() const;
    void set_reference_enclosure_widths(const Float&);
    void set_reference_enclosure_widths(const Vector<Float>&);
    void set_reference_enclosure_widths(const HybridFloatVector&);

    const Float& maximum_enclosure_widths_ratio() const;
    void set_maximum_enclosure_widths_ratio(const Float&);

    const bool& enable_error_rate_enforcement() const;
    void set_enable_error_rate_enforcement(const bool&);

    const bool& enable_boxing_on_contraction() const;
    void set_enable_boxing_on_contraction(const bool&);

    const bool& enable_reconditioning() const;
    void set_enable_reconditioning(const bool&);

    const bool& enable_subdivisions() const;
    void set_enable_subdivisions(const bool&);

    const bool& enable_premature_termination_on_enclosure_size() const;
    void set_enable_premature_termination_on_enclosure_size(const bool&);

    const unsigned int& maximum_number_of_working_sets() const;
    void set_maximum_number_of_working_sets(const unsigned int&);
};

struct EvolutionData {

    typedef ImageSetHybridEvolver::EventListType EventListType;
    typedef ImageSetHybridEvolver::SetModelType SetModelType;
    typedef ImageSetHybridEvolver::TimeModelType TimeModelType;
    typedef ImageSetHybridEvolver::StepType StepType;

    EvolutionData(uint set_index,
                  DiscreteLocation location,
                  StepType previous_step,
                  EventListType previous_events,
                  SetModelType set_model,
                  TimeModelType time_model,
                  uint refinement_width) :
                      _set_index(set_index),
                      _location(location),
                      _previous_step(previous_step),
                      _previous_events(previous_events),
                      _set_model(set_model),
                      _time_model(time_model),
                      _refinement_width(refinement_width){ }

    uint set_index() const { return _set_index; }
    DiscreteLocation location() const { return _location; }
    StepType previous_step() const { return _previous_step; }
    EventListType previous_events() const { return _previous_events; }
    SetModelType set_model() const { return _set_model; }
    TimeModelType time_model() const { return _time_model; }
    uint refinement_width() const { return _refinement_width; }

  private:

    uint _set_index;
    DiscreteLocation _location;
    StepType _previous_step;
    EventListType _previous_events;
    SetModelType _set_model;
    TimeModelType _time_model;
    uint _refinement_width;
};

struct ContinuousStepResult {

    typedef ImageSetHybridEvolver::SetModelType SetModelType;

    ContinuousStepResult(Float used_step, SetModelType flow_set_model, SetModelType finishing_set_model) :
        _used_step(used_step),
        _flow_set_model(flow_set_model),
        _finishing_set_model(finishing_set_model) { }

    Float used_step() const { return _used_step; }
    SetModelType flow_set_model() const { return _flow_set_model; }
    SetModelType finishing_set_model() const { return _finishing_set_model; }

  private:

    Float _used_step;
    SetModelType _flow_set_model;
    SetModelType _finishing_set_model;
};

std::ostream& operator<<(std::ostream& os, const ContinuousStepResult& r);

std::ostream& operator<<(std::ostream& os, const ImageSetHybridEvolverSettings& s);


/*! \brief Whether a box \a set_bounds under a given \a dynamic positively crosses an \a activation. */
tribool positively_crossing(const Box& set_bounds,
							const RealVectorFunction& dynamic,
							const RealScalarFunction& activation);


/*! \brief Set the minimum cell widths from the hybrid grid \a hgrid at \a maximum_grid_depth. */
HybridFloatVector getMinimumGridCellWidths(
		const HybridGrid& hgrid,
		int maximum_grid_depth);

/*! \brief Get the hybrid maximum integration step size, under the assumption that given the maximum derivatives \a hmad,
	all variables in a step must cover a length greater than a length determined by the \a hgrid with a given \a maximum_grid_depth.
	\details The actual result is scaled based on the \a semantics. */
std::map<DiscreteLocation,Float> getHybridMaximumStepSize(
		const HybridFloatVector& hmad,
		const HybridGrid& hgrid,
		int maximum_grid_depth);


} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_IMAGE_H
