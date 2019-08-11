/***************************************************************************
 *            textplot.cc
 *
 *  Copyright 2009  Davide Bresolin
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
 
#include "config.h"

#include "macros.h"
#include "stlio.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "point.h"
#include "box.h"
#include "curve.h"
#include "polytope.h"
#include "textplot.h"
#include "denotable_set.h"
#include "list_set.h"

namespace Ariadne {


TextPlot::~TextPlot()
{
    this->_fstream.close();
}

 
TextPlot::TextPlot()
    : _fstream() 
{ 
}

TextPlot::TextPlot(const char* cfilename)
{ 
    std::string filename(cfilename);
    if(filename.rfind(".") != std::string::npos) {
    } else {
        filename=filename+".txt";
    }
    this->_fstream.open(filename.c_str(), ios::out | ios::trunc);
}


TextPlot::TextPlot(const char* cfilename, ios_base::openmode mode)
{ 
    std::string filename(cfilename);
    if(filename.rfind(".") != std::string::npos) {
    } else {
        filename=filename+".txt";
    }
    this->_fstream.open(filename.c_str(), mode);
}


void TextPlot::open(const char* cfilename) 
{
    std::string filename(cfilename);
    if(filename.rfind(".") != std::string::npos) {
    } else {
        filename=filename+".txt";
    }
    this->_fstream.open(filename.c_str(), ios::out | ios::trunc);
}


void TextPlot::open(const char* cfilename, ios_base::openmode mode) 
{
    std::string filename(cfilename);
    if(filename.rfind(".") != std::string::npos) {
    } else {
        filename=filename+".txt";
    }
    this->_fstream.open(filename.c_str(), mode);
}


void TextPlot::draw(const Point& pt) {
    for(uint i = 0; i < pt.dimension(); i++) {
        this->_fstream << double(pt[i]) << " ";
    }
    this->_fstream << std::endl;
}

void TextPlot::draw(const std::vector<Point>& pts) {
    for(std::vector<Point>::const_iterator iter = pts.begin() ; iter != pts.end() ; iter++) {
        this->draw(*iter);
    }
    this->draw(*pts.begin()); // Adds the initial point again, in order to provide a closed curve
    this->_fstream << std::endl;
}    


void TextPlot::draw(const Box& bx) {
    this->draw(bx.vertices());
}

void TextPlot::draw(const Polytope& p) {
    this->draw(p.vertices());
}

void TextPlot::draw(const InterpolatedCurve& c) {
    for(InterpolatedCurve::const_iterator iter = c.begin() ; iter != c.end() ; ++iter) {
        this->draw(iter->second);
    }
    this->_fstream << std::endl;
}

void TextPlot::draw(const DenotableSetType& ds) {
    for(DenotableSetType::const_iterator iter=ds.begin(); iter!=ds.end(); ++iter) {
        this->draw(iter->box());
    }
}

void TextPlot::close() {
    this->_fstream.close();
}


} // Namespace Ariadne

 
