/***************************************************************************
 *            verification_input.h
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

/*! \file verification_input.h
 *  \brief Classes for holding input for verification routines.
 */

#ifndef ARIADNE_VERIFICATION_INPUT_H_
#define ARIADNE_VERIFICATION_INPUT_H_

#include "hybrid_automaton_interface.h"
#include "parametric.h"
#include "hybrid_set.h"

namespace Ariadne {

/** Provides a basic bundle for systems and their relative input for verification. */
class VerificationInput
{
  public:
	typedef HybridAutomatonInterface SystemType;
  protected:
	SystemType& _system;
  private:
	HybridBoundedConstraintSet& _initial_set;
	HybridBoxes& _domain;

  public:
	SystemType& getSystem() const { return _system; }
	HybridBoundedConstraintSet& getInitialSet() const { return _initial_set; }
	HybridBoxes& getDomain() const { return _domain; }

	VerificationInput(
			SystemType& system,
			HybridBoundedConstraintSet& initial_set,
			HybridBoxes& domain);

	virtual std::ostream& write(std::ostream&) const;

protected:

	/** \brief Checks fields consistency */
	void _check_fields() const;

};

inline std::ostream& operator<<(std::ostream& os, const VerificationInput& verInput) {
    return verInput.write(os); }

/** Provides a convenient bundle for systems and their relative input for safety verification. */
class SafetyVerificationInput : public VerificationInput
{
  private:
	HybridConstraintSet _safety_constraint;

  public:
	const HybridConstraintSet& getSafetyConstraint() const { return _safety_constraint; }

	SafetyVerificationInput(
			SystemType& system,
			HybridBoundedConstraintSet& initial_set,
			HybridBoxes& domain,
			HybridConstraintSet& safety_constraint);

	virtual std::ostream& write(std::ostream&) const;

protected:

	/** \brief Checks fields consistency */
	void _check_fields() const;

};

inline std::ostream& operator<<(
		std::ostream& os,
		const SafetyVerificationInput& verInput)
{
    return verInput.write(os);
}


/** Provides a convenient bundle for systems and their relative input for dominance verification. */
class DominanceVerificationInput : public VerificationInput
{
  private:
	Vector<uint> _projection;

  public:
	const Vector<uint>& getProjection() const { return _projection; }

	DominanceVerificationInput(
			SystemType& system,
			HybridBoundedConstraintSet& initial_set,
			HybridBoxes& domain,
			Vector<uint>& projection);

	virtual std::ostream& write(std::ostream&) const;

};

inline std::ostream& operator<<(std::ostream& os, const DominanceVerificationInput& verInput) {
    return verInput.write(os); }

}


#endif /* ARIADNE_VERIFICATION_INPUT_H_ */
