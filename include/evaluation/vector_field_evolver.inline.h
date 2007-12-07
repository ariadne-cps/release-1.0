/***************************************************************************
 *            vector_field_evolver.inline.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 

#include "../evaluation/bounder.h"

namespace Ariadne {
    
template<class R> template<class BST> inline
Evaluation::VectorFieldEvolver<R>::VectorFieldEvolver(const EvolutionParameters<R>& parameters, const IntegratorInterface<BST>& plugin)
  : _parameters(new EvolutionParameters<R>(parameters)),
    _bounder(new Bounder<R>()),
    _integrator(plugin.clone())
{
}

}
