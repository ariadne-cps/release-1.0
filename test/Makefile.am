all: tests profiles

profile-lib: $(bin_PROGRAMS)
	for file in $(bin_PROGRAMS); do $(top_srcdir)/test/run_profile.sh $$file; done

#TESTS_ENVIRONMENT = (function run_tst { $1 > $1.log; }; run_tst )
#TESTS_ENVIRONMENT = $(top_srcdir)/test/run_test.sh
#PROFILES_ENVIRONMENT = $(top_srcdir)/test/run_profile.sh

#AM_CXXFLAGS = -Wall
AM_LDFLAGS = $(EXEC_LDFLAGS)

AM_CPPFLAGS = -I$(top_srcdir)/include $(CAIRO_CFLAGS) $(GTK_CFLAGS)
LIBS += $(CAIRO_LIBS) $(GTK_LIBS)
LDADD = ../src/libariadne.la $(LIBS)

if USE_BUDDY
## BuDDy library is present, use BDDTreeSet
TEST_DENOTABLE_SET_TYPE = test_bdd_set
else
## BuDDy library is not present (or disables), use GridTreeSet
TEST_DENOTABLE_SET_TYPE = test_grid_set
endif

noinst_PROGRAMS = \
	profile_arithmetic \
	profile_taylor_model \
	profile_taylor_calculus

check_PROGRAMS = \
	test_float \
	test_interval \
	test_vector \
	test_matrix \
	test_multi_index \
	test_expansion \
	test_polynomial \
	test_differential \
	test_taylor_model \
	test_taylor_function \
	test_taylor_set \
	test_function \
	test_function_set \
	$(TEST_DENOTABLE_SET_TYPE) \
	test_hybrid_set \
	test_hybrid_io_automaton \
	test_hybrid_evolution \
	test_discretised_evolution \
	test_reachability_restriction \
	test_reachability_analysis \
	test_serialization \
	test_graphics \
	test_textplot \
	test_misc

#	test_hybrid_simulator

tests: $(check_PROGRAMS)
profiles: $(bin_PROGRAMS)

#test_integer_LDADD = -lgmp $(LIBS)
#test_rational_LDADD = -lgmp $(LIBS)

# Explicitly list profile sources
profile_arithmetic_SOURCES=profile_arithmetic.cc
profile_taylor_model_SOURCES=profile_taylor_model.cc
profile_taylor_calculus_SOURCES=profile_taylor_calculus.cc

# Explicitly list test sources
test_float_SOURCES=test_float.cc
test_interval_SOURCES=test_interval.cc
test_vector_SOURCES=test_vector.cc
test_matrix_SOURCES=test_matrix.cc
test_multi_index_SOURCES=test_multi_index.cc
test_expansion_SOURCES=test_expansion.cc
test_polynomial_SOURCES=test_polynomial.cc
test_differential_SOURCES=test_differential.cc
test_taylor_model_SOURCES=test_taylor_model.cc
test_taylor_function_SOURCES=test_taylor_function.cc
test_taylor_set_SOURCES=test_taylor_set.cc
test_function_SOURCES=test_function.cc
test_function_set_SOURCES=test_function_set.cc
test_grid_set_SOURCES=test_grid_set.cc
test_bdd_set_SOURCES=test_bdd_set.cc
test_hybrid_io_automaton_SOURCES=test_hybrid_io_automaton.cc
test_hybrid_set_SOURCES=test_hybrid_set.cc
test_hybrid_evolution_SOURCES=test_hybrid_evolution.cc
test_discretised_evolution_SOURCES=test_discretised_evolution.cc
test_reachability_restriction_SOURCES=test_reachability_restriction.cc
test_reachability_analysis_SOURCES=test_reachability_analysis.cc
test_serialization_SOURCES=test_serialization.cc
test_graphics_SOURCES=test_graphics.cc
test_textplot_SOURCES=test_textplot.cc
test_misc_SOURCES=test_misc.cc

EXTRA_DIST =
SUBDIRS =

TESTS = $(check_PROGRAMS)
XFAIL_TESTS = 
testsubdir = ./
