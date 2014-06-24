AM_CPPFLAGS = -I$(top_srcdir)/include
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
AM_LDFLAGS = $(EXEC_LDFLAGS)
LDADD = ../src/libariadne.la $(LIBS)

tutorial: $(bin_PROGRAMS)

bin_PROGRAMS= monolithic-linear-forced \
              monolithic-linear-unforced \
              monolithic-nonlinear-unforced \
              compositional-linear-forced \
              compositional-linear-unforced \
              compositional-nonlinear-unforced
		
monolithic_linear_forced_SOURCES= monolithic-linear-forced.cc
monolithic_linear_unforced_SOURCES= monolithic-linear-unforced.cc
monolithic_nonlinear_unforced_SOURCES= monolithic-nonlinear-unforced.cc
compositional_linear_forced_SOURCES= compositional-linear-forced.cc
compositional_linear_unforced_SOURCES= compositional-linear-unforced.cc
compositional_nonlinear_unforced_SOURCES= compositional-nonlinear-unforced.cc


