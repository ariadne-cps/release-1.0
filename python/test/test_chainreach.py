#!/usr/bin/python

from ariadne import *
from evaluation import *
from geometry import *
from linear_algebra import *
import sys


h=HenonMap(Dyadic(1.5),Dyadic(0.75))
gbb=Rectangle("[-20,20]x[-20,20]")
g=FiniteGrid(gbb,2560);
r=Rectangle("[1.49,1.51]x[0.49,0.51]")
cb=Rectangle("[-4,4]x[-4,4]") # cutoff box
epsbb=Rectangle("[-4.1,4.1]x[-4.1,4.1]") # eps bounding box box
i=RectangleListSet(r)

cr=chainreach(h,i,g,cb)

eps=EpsPlot("cr.eps",epsbb)
eps.set_pen_colour("black")
eps.set_fill_colour("white")
eps.write(cb)
eps.set_line_style(0)
eps.set_fill_colour("green")
eps.write(cr)
eps.set_line_style(1)
eps.set_fill_colour("blue")
eps.write(r)
eps.close()
sys.exit()
