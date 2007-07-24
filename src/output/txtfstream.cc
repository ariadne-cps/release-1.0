/****************************************************************************
 *            txtfstream.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins, Davide Bresolin
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl, bresolin@sci.univr.it
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

#include "base/stlio.h"
#include "output/txtfstream.h"

namespace Ariadne { 

Output::txtfstream::txtfstream()
  : std::ofstream()
{
}



Output::txtfstream::~txtfstream() {
  this->close();
}

      
void
Output::txtfstream::open(const char* fn)
{
  this->std::ofstream::open(fn);
}


void 
Output::txtfstream::close() 
{
  this->std::ofstream::close();
}


void 
Output::txtfstream::writenl() 
{
  (*this) << std::endl;
}

}