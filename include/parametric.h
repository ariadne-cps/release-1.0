/***************************************************************************
 *            parametric.h
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

/*! \file parametric.h
 *  \brief Data structures for handling parametric analysis.
 */

#include "tribool.h"
#include "variables.h"

#ifndef ARIADNE_PARAMETRIC_H
#define ARIADNE_PARAMETRIC_H

namespace Ariadne {

class HybridAutomatonInterface;

typedef Parameter<Real> RealParameter;
typedef std::set<RealParameter,ParameterSetComparator<Real> > RealParameterSet;


class Parameterisable
{
  public:
	//! \brief Extracts from the value of the parameter having identifier \a name.
	Real parameter_value(String name) const;

	//! \brief Extracts the subset of parameters whose values are not singleton.
	RealParameterSet nonsingleton_parameters() const;

	/*! \brief Substitute values from a set \a params. */
	void substitute_all(const RealParameterSet& params, bool use_midpoints = false);

	//! \brief Get the set of parameters.
    virtual RealParameterSet parameters() const = 0;

	//! \brief Substitute the parameter \a param where present.
	virtual void substitute(const RealParameter& param) = 0;
};

/**
 * \brief The data structure for the outcome over a configuration of parameters (i.e. constants of a system)
 */
struct ParametricOutcome
{
private:

	RealParameterSet _params;
	tribool _value;

public:

	ParametricOutcome(const RealParameterSet& params, tribool value);
	ParametricOutcome(const ParametricOutcome& other);

	ParametricOutcome& operator=(const ParametricOutcome& other);

	virtual std::ostream& write(std::ostream&) const;

	/** The parameters */
	const RealParameterSet& getParams() const;
	/** The outcome of the verification */
	const tribool& getOutcome() const;

};

inline std::ostream& operator<<(std::ostream& os, const ParametricOutcome& outcome) {
    return outcome.write(os); }

/*! \brief Draws the \a outcomeList in the current folder */
void draw(std::string basename, const std::list<ParametricOutcome>& outcomes);

} // namespace Ariadne

#endif // ARIADNE_PARAMETRIC_H
