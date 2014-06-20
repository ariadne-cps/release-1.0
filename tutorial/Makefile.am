AM_CPPFLAGS = -I$(top_srcdir)/include
LIBS +=  $(BOOST_LIBS) $(CAIRO_LIBS) $(GTK_LIBS)
AM_LDFLAGS = $(EXEC_LDFLAGS)
LDADD = ../src/libariadne.la $(LIBS)

tutorial: $(bin_PROGRAMS)

bin_PROGRAMS=  monolithic-linear-nonurgent-parametric \
        	monolithic-nonlinear-nonurgent-parametric \
		compositional-linear-nonurgent-parametric
		
monolithic_linear_nonurgent_parametric_SOURCES= monolithic-linear-nonurgent-parametric.cc
monolithic_nonlinear_nonurgent_parametric_SOURCES= monolithic-nonlinear-nonurgent-parametric.cc
compositional_linear_nonurgent_parametric_SOURCES= compositional-linear-nonurgent-parametric.cc


