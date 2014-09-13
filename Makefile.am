SUBDIRS = include src test examples tutorial

EXTRA_DIST = Doxyfile

AM_CXXFLAGS = -Wall -fno-strict-aliasing

ACLOCAL_AMFLAGS = -I m4

DOCS_DESTDIR = $(prefix)/share/doc/libariadne-doc
EXAMPLES_DESTDIR = $(DOCS_DESTDIR)

.PHONY: shiny doc all lib examples profiles tutorial check-lib 

all: lib tests examples profiles doc tutorial

# Show the values of the environment variables
evn:
	@echo 'INCLUDES =' $(INCLUDES)
	@echo 'LIBS =' $(LIBS)

shiny: maintainer-clean
	-rm -rf *~ */*~
	-mv doc/evolve.eps doc/evolve.eps.tmp 
	-rm -f *.log */*.log */*.eps
	-mv doc/evolve.eps.tmp doc/evolve.eps
	-rm -rf docs 
	-rm -rf test/*.png test/*.txt test/*.log  
	-rm -rf */.deps */.libs
	-rm -f aclocal.m4 config.* configure depcomp install-sh libtool ltmain.sh missing mkinstalldirs stamp-h* .tm_project.cache
	-rm -rf autom4te.cache
	-rm -f Makefile.in */Makefile.in
	-rm -f */Makefile
	-rm -f Makefile

if HAVE_DOXYGEN
doc: docs docs/refman.pdf docs/html

docs: 
	doxygen 

docs/html: docs

docs/refman.pdf: docs
	(cd docs/latex; make)
	cp docs/latex/refman.pdf docs/refman.pdf 

install-doc: doc
	install -d $(DOCS_DESTDIR)
	-cp -rp docs/refman.pdf $(DOCS_DESTDIR)/ariadne-manual.pdf
	-cp -rp docs/html $(DOCS_DESTDIR)
else
doc:

install-doc:	
endif

lib:
	(cd src; make lib)

tests: lib
	(cd test; make tests)

profiles: lib
	(cd test; make profiles)

check-lib: lib
	(cd test; make check)

profile-lib: lib profiles
	(cd test; make profile-lib)

install-lib: lib
	(cd src; make install)

install-examples: install-cxx-examples 

install-cxx-examples:
	install -d $(EXAMPLES_DESTDIR)
	-mv examples/Makefile examples/Makefile.orig
	-cp examples/Makefile-examples examples/Makefile
	tar -zcf $(EXAMPLES_DESTDIR)/c++-examples.tar.gz examples/*.h \
		examples/*.cc examples/*.plt examples/Makefile-bin.am \
		examples/Makefile 
	-mv examples/Makefile.orig examples/Makefile

install-am: install-doc
