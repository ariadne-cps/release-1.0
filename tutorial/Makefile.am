AM_CPPFLAGS = -I$(top_srcdir)/include
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
AM_LDFLAGS = $(EXEC_LDFLAGS)
LDADD = ../src/libariadne.la $(LIBS)

tutorial: $(bin_PROGRAMS)

bin_PROGRAMS= monolithic-forced \
              monolithic-unforced \
              compositional
		
monolithic_forced_SOURCES= monolithic-forced.cc
monolithic_unforced_SOURCES= monolithic-unforced.cc
compositional_SOURCES= compositional.cc

