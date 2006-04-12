#!/usr/bin/python

##############################################################################
#            test_integration_step.py
#
#  Copyright 2006  Pieter Collins <Pieter.Collins@cwi.nl>
##############################################################################

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

from ariadne import Real
from ariadne.base import *
from ariadne.numeric import *
from ariadne.evaluation import *
from ariadne.geometry import *
from ariadne.linear_algebra import *
import sys

def plot(fn,bb,set):
  eps=EpsPlot(fn+".eps",bb)
  eps.set_fill_colour("green")
  for bs in set:
    eps.write(bs)
  eps.set_fill_colour("blue")
 # eps.write(set[0])
  eps.set_fill_colour("red")
#  eps.write(set[-1])
  eps.close()

bb=Rectangle("[-2,2]x[-2,2]")

A=Matrix("[-0.25,-1;+1,-0.25]")
b=Vector("[0,0]")
avf=AffineVectorField(A,b)

init_rect=Rectangle("[0.96,1.04]x[0.46,0.54]")
init_paral=Parallelotope(init_rect)

Real.__repr__=Real.__str__
reach_list=[init_paral]
step_size=Real(0.125)
step_sizes=[step_size]
print Rational(step_size)
for i in range(0,128):
  print "before:",Rational(step_size)
  reach_list.append(integration_step(avf,reach_list[-1],step_size))
  print Rational(step_size)
  step_sizes.append(step_size)
  step_size=min(step_size*2,step_sizes[0])
print Rational(step_size)
plot("integrationstep1",bb,reach_list)
print step_sizes, len(step_sizes)