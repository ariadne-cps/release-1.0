/***************************************************************************
 *            python/export_text_output.cc
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins, Davide Bresolin
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl, bresolin@sci.univr.it
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#include "python/python_float.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"

#include "output/textstream.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class S> inline void write(textfstream& txt, const S& s) { txt << s; }
template<class R> inline void write_rectangle(textfstream& txt, const Rectangle<R>& r) { txt << r; }
template<class R> inline void write_rectangular_set(textfstream& txt, const RectangularSet<R>& r) { txt << r; }
template<class R> inline void write_parallelotope(textfstream& txt, const Parallelotope<R>& p) { txt << p; }
template<class R0,class R1> inline void write_zonotope(textfstream& txt, const Zonotope<R0,R1>& z) { txt << z; }
template<class R> inline void write_polytope(textfstream& txt, const Polytope<R>& p) { txt << p; }
template<class R> inline void write_polyhedron(textfstream& txt, const Polyhedron<R>& p) { txt << p; }
template<class R> inline void write_polyhedral_set(textfstream& txt, const PolyhedralSet<R>& p) { txt << p; }
template<class BS> inline void write_list_set(textfstream& txt, const ListSet<BS>& ls) { txt << ls; }
template<class R> inline void write_polytope_list_set(textfstream& txt, const ListSet< Polytope<R> >& s) { txt << s; }
template<class R> inline void write_grid_cell(textfstream& txt, const GridCell<R>& r) { txt << Rectangle<R>(r); }
template<class R> inline void write_grid_block(textfstream& txt, const GridBlock<R>& r) { txt << Rectangle<R>(r); }
template<class R> inline void write_grid_cell_list_set(textfstream& txt, const GridCellListSet<R>& s) { txt << s; }
template<class R> inline void write_grid_mask_set(textfstream& txt, const GridMaskSet<R>& s) { txt << s; }
template<class R> inline void write_partition_tree_set(textfstream& txt, const PartitionTreeSet<R>& s) { txt << s; }
template<class R> inline void write_finite_grid(textfstream& txt, const FiniteGrid<R>& fg) { txt << fg; }
template<class R> inline void write_partition_tree(textfstream& txt, const PartitionTree<R>& s) { txt << s; }
template<class R> inline void textfstream_open(textfstream& txt) { txt.open("Ariadne"); }
inline void textfstream_close(textfstream& txt) { txt.close(); }

void export_text_output()
{
    
  class_<textfstream, boost::noncopyable>("TextFile",init<>())
    .def("open",(void(textfstream::*)(const char* fn))&textfstream::open)
    .def("close",&textfstream_close)
		.def("write",&write< Rectangle<Float> >)
    .def("write",&write< RectangularSet<Float> >)
    .def("write",&write< Parallelotope<Float> >)
    .def("write",&write< Zonotope<Float,Float> >)
    .def("write",&write< Zonotope<Interval<Float>,Float> >)
    .def("write",&write< Parallelotope<Float> >)
    .def("write",&write< Polytope<Float> >)
    .def("write",&write< Polyhedron<Float> >)
    .def("write",&write< PolyhedralSet<Float> >)
    .def("write",&write< ListSet< Rectangle<Float> > >)
    .def("write",&write< ListSet< Parallelotope<Float> > >)
    .def("write",&write< ListSet< Polytope<Float> > >)
    .def("write",&write< ListSet< Zonotope<Float,Float> > >)
    .def("write",&write< ListSet< Zonotope<Interval<Float>,Float> > >)
    .def("write",&write< GridCell<Float> >)
    .def("write",&write< GridBlock<Float> >)
    .def("write",&write< GridCellListSet<Float> >)
    .def("write",&write< GridMaskSet<Float> >)
    .def("write",&write< PartitionTreeSet<Float> >)
    .def("write",&write< FiniteGrid<Float> >)
    .def("write",&write< PartitionTree<Float> >)
  ;
  
}