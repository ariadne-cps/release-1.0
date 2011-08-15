/***************************************************************************
 *            parametric.cc
 *
 *  Copyright 2010  Luca Geretti
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

#include "hybrid_automaton_interface.h"
#include "parametric.h"
#include "box.h"
#include "textplot.h"
#include "graphics.h"

namespace Ariadne {

Real Parameterisable::parameter_value(String name) const
{
	RealParameterSet parameters = this->parameters();
	for (RealParameterSet::const_iterator parameter_it = parameters.begin();
												 parameter_it != parameters.end();
												 ++parameter_it) {
		if (parameter_it->name() == name)
			return parameter_it->value();
	}

	ARIADNE_FAIL_MSG("The parameter '" << name << "' was not found in the object.");
}

RealParameterSet
Parameterisable::nonsingleton_parameters() const
{
	RealParameterSet parameters = this->parameters();

	RealParameterSet result;
	for (RealParameterSet::const_iterator param_it = parameters.begin(); param_it != parameters.end(); ++param_it) {
		if (!param_it->value().singleton())
			result.insert(*param_it);
	}

	return result;
}


void
Parameterisable::substitute_all(const RealParameterSet& params, bool use_midpoints)
{
	for (RealParameterSet::const_iterator param_it = params.begin(); param_it != params.end(); ++param_it) {
		if (use_midpoints)
			substitute(RealParameter(param_it->name(),param_it->value().midpoint()));
		else
			substitute(*param_it);
	}
}


ParametricOutcome::ParametricOutcome(const RealParameterSet& params, tribool value)
{
	for (RealParameterSet::const_iterator param_it = params.begin(); param_it != params.end(); ++param_it)
		_params.insert(*param_it);

	_value = value;
}

ParametricOutcome::ParametricOutcome(const ParametricOutcome& other)
{
	for (RealParameterSet::const_iterator param_it = other.getParams().begin(); param_it != other.getParams().end(); ++param_it)
		_params.insert(*param_it);

	_value = other.getOutcome();
}

ParametricOutcome&
ParametricOutcome::operator=(const ParametricOutcome& other)
{
	for (RealParameterSet::const_iterator param_it = other.getParams().begin(); param_it != other.getParams().end(); ++param_it)
		_params.insert(*param_it);

	_value = other.getOutcome();

	return *this;
}

const RealParameterSet&
ParametricOutcome::getParams() const
{
	return _params;
}

const tribool&
ParametricOutcome::getOutcome() const
{
	return _value;
}

std::ostream&
ParametricOutcome::write(std::ostream& os) const
{
	os << "(" << _params << "->" << tribool_pretty_print(_value) << ")";
	return os;
}


void
draw(std::string basename, const std::list<ParametricOutcome>& outcomes)
{
	RealParameterSet _params = outcomes.back().getParams();

	ARIADNE_ASSERT_MSG(outcomes.size() > 0, "The outcomes list is empty.");
	ARIADNE_ASSERT_MSG(_params.size() > 1, "At least two parameters are required for drawing.");

	// Plots for each couple of parameters
	for (RealParameterSet::const_iterator xparam_it = _params.begin(); xparam_it != _params.end(); ++xparam_it) {

		RealParameterSet::const_iterator yparam_it = xparam_it;
		for (++yparam_it; yparam_it != _params.end(); ++yparam_it)
		{
			std::string currentname = basename + "[" + xparam_it->name() + ","
													 + yparam_it->name() + "]";

			TextPlot trueTxt((currentname + ".true.dump").c_str());
			TextPlot falseTxt((currentname + ".false.dump").c_str());
			TextPlot indeterminateTxt((currentname + ".indeterminate.dump").c_str());

			// Sets up the figure
			Figure fig;
			Box graphics_box(2);
			const RealParameterSet& params = outcomes.begin()->getParams();
			for (RealParameterSet::const_iterator param_it = params.begin(); param_it != params.end(); ++param_it) {
				if (param_it->name() == xparam_it->name())
					graphics_box[0] = param_it->value();
				if (param_it->name() == yparam_it->name())
					graphics_box[1] = param_it->value();
			}
			array<uint> xy(2,0,1);
			fig.set_projection_map(ProjectionFunction(xy,2));

			// Adds each outcome with a dedicated fill colour
			for (std::list<ParametricOutcome>::const_iterator outcome_it = outcomes.begin();
																		  outcome_it != outcomes.end();
																		  ++outcome_it) {
				Box outcome_box(2);
				const RealParameterSet& params = outcome_it->getParams();
				for (RealParameterSet::const_iterator param_it = params.begin(); param_it != params.end(); ++param_it) {
					if (param_it->name() == xparam_it->name())
						outcome_box[0] = param_it->value();
					if (param_it->name() == yparam_it->name())
						outcome_box[1] = param_it->value();
				}

				graphics_box[0].set_lower(min(graphics_box[0].lower(),outcome_box[0].lower()));
				graphics_box[0].set_upper(max(graphics_box[0].upper(),outcome_box[0].upper()));
				graphics_box[1].set_lower(min(graphics_box[1].lower(),outcome_box[1].lower()));
				graphics_box[1].set_upper(max(graphics_box[1].upper(),outcome_box[1].upper()));

				// Chooses the fill color and dumps the box
				tribool outcome = outcome_it->getOutcome();
				if (definitely(outcome)) {
					trueTxt.draw(outcome_box);
					fig.set_fill_colour(Colour(0.0,0.83,0.0));
				}
				else if (indeterminate(outcome)) {
					indeterminateTxt.draw(outcome_box);
					fig.set_fill_colour(Colour(1.0,1.0,0.0));
				}
				else {
					falseTxt.draw(outcome_box);
					fig.set_fill_colour(Colour(1.0,0.34,0.34));
				}

				fig.draw(outcome_box);
			}

			fig.set_bounding_box(graphics_box);
			fig.write(currentname.c_str());

			trueTxt.close();
			indeterminateTxt.close();
			falseTxt.close();
		}
	}
}


} //namespace Ariadne
