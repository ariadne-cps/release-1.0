AM_CXXFLAGS = -Wall

INCLUDES = -I$(top_srcdir)/include
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
AM_LDFLAGS = $(EXEC_LDFLAGS)
LDADD = ../src/libariadne.la $(LIBS)

include Makefile-bin.am

examples: $(bin_PROGRAMS)

