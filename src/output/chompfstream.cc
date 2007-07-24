/****************************************************************************
 *            chompfstream.cc
 *
 *  Copyright  2007  Pieter Collins, Davide Bresolin
 *  Pieter.Collins@cwi.nl, bresolin@sci.univr.it
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

#include "output/chompfstream.h"

namespace Ariadne { 



Output::chompfstream::chompfstream()
 : _ofs()
{
}



Output::chompfstream::chompfstream(const char* fn)
 : _ofs(fn)
{
}



Output::chompfstream::~chompfstream() {
  this->close();
}

      
void
Output::chompfstream::open(const char* fn)
{
  this->_ofs.open(fn);
}


void 
Output::chompfstream::close() 
{
  this->_ofs.close();
}



Output::chompfstream& 
Output::operator<<(chompfstream& cfs, const char* str) 
{
  cfs._ofs << str; return cfs;
}

    
Output::chompfstream& 
Output::operator<<(chompfstream& cfs, const Combinatoric::LatticeCell& lc) 
{
  std::ofstream& ofs=cfs._ofs;
  if(lc.dimension()>0) {
    ofs << "(" << lc.lower_bound(0);
    for(dimension_type i=1; i!=lc.dimension(); ++i) {
      ofs << "," << lc.lower_bound(i);
    }
    ofs << ")";
  } else {
    ofs << "()";
  }
  return cfs;
}


Output::chompfstream& 
Output::operator<<(chompfstream& cfs, const Combinatoric::LatticeMaskSet& lms) 
{
  Combinatoric::LatticeCell lc;
  for(Combinatoric::LatticeMaskSet::const_iterator iter=lms.begin(); iter!=lms.end(); ++iter) {
    lc=*iter;
    cfs << lc;
    cfs._ofs << "\n";
  }
  return cfs;
}


Output::chompfstream& 
Output::operator<<(chompfstream& cfs, const Combinatoric::LatticeMultiMap& lmm) 
{
  Combinatoric::LatticeCell lc(lmm.argument_dimension());
  Combinatoric::LatticeCellListSet lcls(lmm.result_dimension());
  for(Combinatoric::LatticeMultiMap::const_iterator iter=lmm.begin(); iter!=lmm.end(); ++iter) {
    lc=iter->first;
    lcls=iter->second;
    cfs << lc << " -> { ";
    for(Combinatoric::LatticeCellListSet::const_iterator lcls_iter=lcls.begin(); 
        lcls_iter!=lcls.end(); ++lcls_iter) 
      {
        cfs << *lcls_iter << " ";
      }
    cfs << "}\n";
  }
  return cfs;
}


    
}