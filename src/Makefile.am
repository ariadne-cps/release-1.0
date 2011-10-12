#AM_CXXFLAGS = -O2 -Wall

phony: lib

lib: libariadne.la libcalculus.la

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = ariadne.pc

AM_CXXFLAGS = -Wall -fno-strict-aliasing

INCLUDES = -I$(top_srcdir)/include $(CAIRO_CFLAGS) $(GTK_CFLAGS)
LIBS +=  $(CAIRO_LIBS) $(GTK_LIBS)
LDADD =  $(LIBS) 

SRC = \
	numeric.cc \
	vector.cc \
	matrix.cc \
	linear_programming.cc \
	differential.cc \
	taylor_model.cc  \
	expression.cc \
	affine.cc \
	polynomial.cc \
	function.cc \
	taylor_function.cc  \
	simulation_toolbox.cc \
	taylor_calculus.cc \
	solver.cc \
	integrator.cc \
	point.cc \
	box.cc \
	zonotope.cc \
	polytope.cc \
	polyhedron.cc \
	curve.cc \
	function_set.cc \
	taylor_set.cc \
	grid.cc \
	denotable_set.cc \
	hybrid_set.cc \
	geometry.cc \
	hybrid_automaton.cc \
	hybrid_io_automaton.cc \
	orbit.cc \
	hybrid_evolver-image.cc \
	set_checker.cc \
	reachability_restriction.cc \
	reachability_analyser.cc \
	logging.cc \
	graphics.cc \
	textplot.cc \
	parametric.cc \
	verifier.cc \
	verification_input.cc

if USE_BUDDY
## BuDDy library is present, use BDDTreeSet
DENOTABLE_SRC = bdd_set.cc
else
## BuDDy library is not present (or disables), use GridTreeSet
DENOTABLE_SRC = grid_set.cc
endif

CALCULUS_SRC = \
	numeric.cc \
	vector.cc \
	matrix.cc \
	function.cc \
	polynomial.cc \
	differential.cc \
	taylor_model.cc  \
	taylor_set.cc  \
	taylor_function.cc  \
	taylor_calculus.cc

lib_LTLIBRARIES = libariadne.la libcalculus.la
libariadne_la_SOURCES = $(SRC) $(DENOTABLE_SRC)
libcalculus_la_SOURCES = $(CALCULUS_SRC)
nodist_EXTRA_libariadne_la_SOURCES = dummy.cc

EXTRA_DIST =
bin_PROGRAMS =
include_HEADERS = 

rounding_CXXFLAGS = -O0
