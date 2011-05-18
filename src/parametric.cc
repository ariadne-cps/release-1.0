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

#include "parametric.h"

namespace Ariadne {



ParametricOutcome::ParametricOutcome(const RealConstantSet params, const tribool value)
{
	for (RealConstantSet::const_iterator const_it = params.begin(); const_it != params.end(); ++const_it)
		_params.insert(*const_it);

	_value = value;
}

ParametricOutcome::ParametricOutcome(const ParametricOutcome& other)
{
	for (RealConstantSet::const_iterator const_it = other.getParams().begin(); const_it != other.getParams().end(); ++const_it)
		_params.insert(*const_it);

	_value = other.getOutcome();
}

ParametricOutcome&
ParametricOutcome::operator=(const ParametricOutcome& other)
{
	for (RealConstantSet::const_iterator const_it = other.getParams().begin(); const_it != other.getParams().end(); ++const_it)
		_params.insert(*const_it);

	_value = other.getOutcome();

	return *this;
}

const RealConstantSet&
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
	os << "(" << _params << "->" << pretty_print(_value) << ")";
	return os;
}


void
draw(std::string basename, const std::list<ParametricOutcome>& outcomes)
{
	RealConstantSet _params = outcomes.back().getParams();

	ARIADNE_ASSERT_MSG(outcomes.size() > 0, "The outcomes list is empty.");
	ARIADNE_ASSERT_MSG(_params.size() > 1, "At least two parameters are required for drawing.");

	// Plots for each couple of parameters
	for (RealConstantSet::const_iterator xparam_it = _params.begin(); xparam_it != _params.end(); ++xparam_it) {

		RealConstantSet::const_iterator yparam_it = xparam_it;
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
			graphics_box[0] = outcomes.begin()->getParams().find(*xparam_it)->value();
			graphics_box[1] = outcomes.begin()->getParams().find(*yparam_it)->value();
			array<uint> xy(2,0,1);
			fig.set_projection_map(ProjectionFunction(xy,2));

			// Adds each outcome with a dedicated fill colour
			for (std::list<ParametricOutcome>::const_iterator outcome_it = outcomes.begin();
																		  outcome_it != outcomes.end();
																		  ++outcome_it) {
				Box outcome_box(2);
				outcome_box[0] = outcome_it->getParams().find(*xparam_it)->value();
				outcome_box[1] = outcome_it->getParams().find(*yparam_it)->value();

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
