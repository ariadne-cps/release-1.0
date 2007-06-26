#!/bin/python
#
# integration example
#
# Written by Davide Bresolin on May, 2007
#

from ariadne import *
from math import *

vector=extract_vector
matrix=extract_matrix

# y' = z & z' = 3z + 2y
dyn=AffineVectorField(matrix([[0,0,0],[0,0,1],[0,2,3]]),vector([1,0,0]))

grid=Grid(Vector("[0.1,0.1,0.1]"))
block=LatticeBlock("[-2,12]x[-12,22]x[-2,82]")
fgrid=FiniteGrid(grid,block)
print fgrid

print "Creating reference set"
reference_set=GridMaskSet(fgrid)
for i in range(11):
	x = 0.1*i+0.001
	y = exp(2*x)-2*exp(x)
	pstr = "["+str(x)+","+str(x)+"]x["+str(y)+","+str(y)+"]x[0.001,0.001]"
	point = Rectangle(pstr)
	print point
	reference_set.adjoin_outer_approximation(point)
	print i,":",pstr
	
print "reference_set.size(),capacity()=",reference_set.size(),reference_set.capacity()

print "Creating inital set"
initial_rectangle=Rectangle("[0.001,0.001]x[-0.999,-0.999]x[0.001,0.001]")
initial_set=GridMaskSet(fgrid)
print "Adjoining initial rectangle"
initial_set.adjoin_outer_approximation(initial_rectangle)
print "initial_set.size(),capacity()=",initial_set.size(),initial_set.capacity()

bounding_box=Rectangle("[-0.0,0.5]x[-1,2]x[0,8]")

print "Creating bounding set"
bounding_set=GridMaskSet(fgrid)
print "Adjoining bounding rectangle"
bounding_set.adjoin_over_approximation(bounding_box)
print "bounding_set.size(),capacity()=",bounding_set.size(),bounding_set.capacity()

maximum_step_size=0.005;
lock_to_grid_time=2;
maximum_set_radius=0.1;

integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius)

#set_integrator_verbosity(6)
#set_hybrid_evolver_verbosity(6)

print "Computing chainreach sets..."
chainreach_set=integrator.chainreach(dyn,initial_set,bounding_set)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to postscript output...",
eps=EpsPlot()
eps.open("pr36-2-0005.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."


  	
