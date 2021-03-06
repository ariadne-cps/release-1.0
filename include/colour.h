/***************************************************************************
 *            colour.h
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
 
/*! \file colour.h
 *  \brief Colours for graphics objects.
 */

#ifndef ARIADNE_COLOUR_H
#define ARIADNE_COLOUR_H

#include <iosfwd>
#include <string>

typedef unsigned int uint;

namespace Ariadne {


struct Colour {
    Colour();
    Colour(double rd, double gr, double bl, bool tr=true);
    Colour(const char* nm, double rd, double gr, double bl, bool tr=true);
    std::string name;
    double red, green, blue;
    bool transparant;
};

std::ostream& operator<<(std::ostream& os, const Colour& c);

extern const Colour transparant;

extern const Colour white;
extern const Colour black;
extern const Colour red;
extern const Colour green;
extern const Colour blue;
extern const Colour yellow;
extern const Colour cyan;
extern const Colour magenta;


} // namespace Ariadne

#endif // ARIADNE_COLOUR_H
