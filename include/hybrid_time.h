/***************************************************************************
 *            hybrid_time.h
 *
 *  Copyright 2008  Pieter Collins
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
 
/*! \file hybrid_time.h
 *  \brief Hybrid times
 */

#ifndef ARIADNE_HYBRID_TIME_H
#define ARIADNE_HYBRID_TIME_H


namespace Ariadne {

//! \brief A value in a hybrid time domain, being a pair comprising a real \a continuous_time 
//! and an integer \a discrete_time. 
//!
//! When a %HybridTime is used to define a bound on a hybrid evolution, the evolution should
//! stop when <em>either</em> the continuous time or the discrete time reaches the bounding
//! value. This is to ensure that the evolution time is finite; in particular, that no
//! Zeno behaviour occurs.
struct HybridTime
{
    //! \brief The continuous (real, physical) time.
    double continuous_time;
    //! \brief The number of discrete steps taken.
    int discrete_time;
  public:
    HybridTime(double t, int n)
        : continuous_time(t), discrete_time(n) { } 
};

inline bool operator==(const HybridTime& ht1, const HybridTime& ht2) {
    return ht1.continuous_time==ht2.continuous_time &&
        ht1.discrete_time==ht2.discrete_time; 
}

inline bool operator!=(const HybridTime& ht1, const HybridTime& ht2) {
    return ht1.continuous_time!=ht2.continuous_time ||
        ht1.discrete_time!=ht2.discrete_time; 
}

inline bool operator<=(const HybridTime& ht1, const HybridTime& ht2) {
    return ht1.continuous_time<=ht2.continuous_time &&
        ht1.discrete_time<=ht2.discrete_time; 
}

inline bool operator<(const HybridTime& ht1, const HybridTime& ht2) {
    return (ht1<=ht2) && (ht1 != ht2);
}

inline std::ostream& operator<<(std::ostream& os, const HybridTime& ht) {
    return os << "("<<ht.continuous_time<<","<<ht.discrete_time<<")";
}

} // namespace Ariadne

#endif // ARIADNE_HYBRID_TIME_H
