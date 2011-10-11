/***************************************************************************
 *            verification_input.cc
 *
 *  Copyright 2011  Luca Geretti
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

#include <sys/stat.h>
#include <sys/types.h>
#include <string>

#include "verification_input.h"
#include "hybrid_set.h"

namespace Ariadne {


VerificationInput::VerificationInput(
		SystemType& system,
		HybridBoundedConstraintSet& initial_set,
		HybridBoxes& domain) :
		_system(system),
		_initial_set(initial_set),
		_domain(domain)
{
	_check_fields();
}


void VerificationInput::_check_fields() const
{
	HybridSpace hspace = _system.state_space();
	ARIADNE_ASSERT_MSG(hspace == _initial_set.space(), "The initial set space and the system space do not match.");

	for (HybridSpace::const_iterator hspace_it = hspace.begin(); hspace_it != hspace.end(); ++hspace_it) {
		HybridBoxes::const_iterator domain_it = _domain.find(hspace_it->first);
		ARIADNE_ASSERT_MSG(domain_it != _domain.end(),
						   "The location " << hspace_it->first.name() << "is not present into the domain.");
		ARIADNE_ASSERT_MSG(hspace_it->second == domain_it->second.dimension(),
						   "The dimension of the continuous space in the domain for location " << hspace_it->first.name() << " does not match the system space");
	}
}


std::ostream&
VerificationInput::write(std::ostream& os) const
{
	os << "(System: " << _system << "; Initial set: " << _initial_set << "; Domain: " << _domain << ")";
	return os;
}


SafetyVerificationInput::SafetyVerificationInput(
		SystemType& system,
		HybridBoundedConstraintSet& initial_set,
		HybridBoxes& domain,
		HybridConstraintSet& safety_constraint) :
		VerificationInput(system,initial_set,domain),
	    _safety_constraint(safety_constraint)
{
	_check_fields();
}


void SafetyVerificationInput::_check_fields() const
{
	VerificationInput::_check_fields();

	ARIADNE_ASSERT_MSG(_system.state_space() == _safety_constraint.space(),
			"The system space and the constraint space do not match.");
}


std::ostream&
SafetyVerificationInput::write(std::ostream& os) const
{
	os << "(System: " << getSystem() << "; Initial set: " << getInitialSet() << "; Domain: " <<
			getDomain() << "; Safety constraint: " << _safety_constraint << ")";
	return os;
}


DominanceVerificationInput::DominanceVerificationInput(
		SystemType& system,
		HybridBoundedConstraintSet& initial_set,
		HybridBoxes& domain,
		Vector<uint>& projection) :
		VerificationInput(system,initial_set,domain),
		_projection(projection)
{
}


std::ostream&
DominanceVerificationInput::write(std::ostream& os) const
{
	os << "(System: " << getSystem() << "; Initial set: " << getInitialSet() << "; Domain: " <<
			getDomain() << "; Projection: " << _projection << ")";
	return os;
}



}
