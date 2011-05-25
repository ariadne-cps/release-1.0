/***************************************************************************
 *            discrete_location.h
 *
 *  Copyright  2004-9  Alberto Casagrande, Pieter Collins
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

/*! \file discrete_location.h
 *  \brief Class representing a discrete location.
 */

#ifndef ARIADNE_DISCRETE_LOCATION_H
#define ARIADNE_DISCRETE_LOCATION_H

#include "container.h"

namespace Ariadne {

class DiscreteLocation {
  public:
    DiscreteLocation() : _id("q?") { }
    DiscreteLocation(int n) : _id(std::string("q"+to_str(n))) { }
    DiscreteLocation(const std::string& s) : _id(s) { }
    std::string name() const { return this->_id; }
    bool operator==(const DiscreteLocation& q) const { return this->_id==q._id; }
    bool operator!=(const DiscreteLocation& q) const { return this->_id!=q._id; }
    bool operator<=(const DiscreteLocation& q) const { return this->_id<=q._id; }
    bool operator>=(const DiscreteLocation& q) const { return this->_id>=q._id; }
    bool operator< (const DiscreteLocation& q) const { return this->_id< q._id; }
    bool operator> (const DiscreteLocation& q) const { return this->_id> q._id; }
    friend std::ostream& operator<<(std::ostream& os, const DiscreteLocation& q);
  private:
    std::string _id;
};

inline std::ostream& operator<<(std::ostream& os, const DiscreteLocation& q) {
    return os << q._id; }

template<class A> inline void serialize(A& archive, DiscreteLocation& state, const uint version) {
    std::string& id=reinterpret_cast<std::string&>(state);
    archive & id;
}

} //namespace Ariadne

#endif /* ARIADNE_DISCRETE_LOCATION_H */
