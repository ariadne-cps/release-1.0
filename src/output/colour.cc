/****************************************************************************
 *            colour.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include <iostream>

#include "output/colour.h"

namespace Ariadne {

std::ostream& 
Output::operator<<(std::ostream& os, const Colour& c) 
{
  return os << "Colour( name=" << c.name() << ", r=" << c.red() << ", g=" << c.green() << ", b=" << c.blue() << " )";
}


const Output::Colour Output::transparant=Colour();

const Output::Colour Output::white=Colour("white",255,255,255);
const Output::Colour Output::black=Colour("black",0,0,0);
const Output::Colour Output::red=Colour("red",255,0,0);
const Output::Colour Output::green=Colour("green",0,255,0);
const Output::Colour Output::blue=Colour("blue",0,0,255);
const Output::Colour Output::yellow=Colour("yellow",255,255,0);
const Output::Colour Output::cyan=Colour("cyan",0,255,255);
const Output::Colour Output::magenta=Colour("magenta",255,0,255);

}
