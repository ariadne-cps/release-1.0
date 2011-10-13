/***************************************************************************
 *            interruptible.h
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

/*! \file interruptible.h
 *  \brief Support for handling timeouts and throwing exceptions.
 */

#include <time.h>

#ifndef ARIADNE_INTERRUPTIBLE_H_
#define ARIADNE_INTERRUPTIBLE_H_

namespace Ariadne {


//! \brief Exception for handling timeouts.
class TimeoutException : public std::runtime_error {
  public:
	TimeoutException() : std::runtime_error("TimeoutException raised due to TTL being hit.") { }
};


//! \brief Class for holding data related to timeout checking.
class Interruptible {

  public:

	//! \brief Default constructor has a practically unlimited TTL
	Interruptible() : ttl(std::numeric_limits<uint>::max()) { }

  protected:

	//! \brief Reset the start time to the current time.
	void _reset_start_time() const { _start_time = time(NULL); }

	//! \brief Probe the remaining time.
	uint _remaining_time() const { return ttl - (time(NULL) - _start_time); }

	//! \brief Check for a timeout.
	void _check_timeout() const { if (_remaining_time() <= 0) throw TimeoutException(); }

  public:

	//! \brief The Time-To-Live, i.e. the time (in seconds) that any instrumented method is expected to return within,
	//! or raise a %TimeoutException.
	uint ttl;

  protected:

	//! \brief An internal variable to hold the start time reference.
    mutable uint _start_time;

};


}

#endif /* ARIADNE_INTERRUPTIBLE_H_ */
