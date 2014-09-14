AM_CPPFLAGS = -I$(top_srcdir)/include
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
AM_LDFLAGS = $(EXEC_LDFLAGS)
LDADD = ../src/libariadne.la $(LIBS)

tutorial: $(bin_PROGRAMS)

noinst_PROGRAMS = monolithic compositional
		
monolithic_SOURCES = monolithic.cc
compositional_SOURCES= compositional.cc

