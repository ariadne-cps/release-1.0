/***************************************************************************
 *            main_page.h
 *
 *  Copyright  2004-7  Pieter Collins
 *  Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! 

\file main_page.h
\brief Main page of Doxygen documentation



\mainpage

\section Introduction

%Ariadne is a C++ package for set-based analysis of dynamical and control systems, including reachability analysis and verification.

\section Download

The homepage of %Ariadne is <a href="http://fsv.dimi.uniud.it/ariadne/">http://fsv.dimi.uniud.it/ariadne/</a>.

You can check out the latest version on the Subversion repository by typing:

  \code svn checkout https://fsv.dimi.uniud.it/svn/ariadne ariadne \endcode

To make the code documentation, change to the ariadne/trunk/ directory and type:

  \code doxygen \endcode

\section Requirements

To compile %Ariadne, you will need:
  - A C++ compiler (we recommend g++ version 4.0.2 or higher, which can be downloaded from <a href="http://gcc.gnu.org/">http://gcc.gnu.org</a>) 
  - A version of the 'make' utility (such as GNU make <a href="http://www.gnu.org/software/make/">http://www.gnu.org/software/make/</a>)

You will also need the following libraries:
  - The GNU Multiple-Precision Library (version 4.1.2 or higher)  <a href="http://www.swox.com/gmp/">http://www.swox.com/gmp/</a>.
  - MPFR Library (version 2.2.1 or higher) <a href="http://www.mpfr.org/">http://www.mpfr.org/</a>.
  - The Boost C++ Libraries (version 1.33.1) <a href="http://www.boost.org/">http://www.boost.org/</a>.
  - The Parma Polyhedra Library (version 0.9 or higher) <a href="http://www.cs.unipr.it/ppl/">http://www.cs.unipr.it/ppl/</a>.
  - TBLAS (version 0.4.1 or higher) <a href="http://homepages.cwi.nl/~collins/software/">http://homepages.cwi.nl/~collins/software/</a>.
  - TLAPACK (version 0.4.0 or higher) <a href="http://homepages.cwi.nl/~collins/software/">http://homepages.cwi.nl/~collins/software/</a>.

For the Python interface, you will also need:
  - Python (version 2.4) <a href="http://www.python.org/">http://www.python.org/</a>.

To make the source code documentation, you will need:
  - Doxygen (version 1.4.6 or higher is recommended) 
       <a href="http://www.stack.nl/~dimitri/doxygen/">http://www.stack.nl/~dimitri/doxygen/</a>

To build the source code from the Subversion repository, you will also need
  - GNU autotools
      - autoconf <a href="http://www.gnu.org/software/autoconf/">http://www.gnu.org/software/autoconf/</a>.
      - automake <a href="http://www.gnu.org/software/automake/">http://www.gnu.org/software/automake/</a>.
      - libtool <a href="http://www.gnu.org/software/libtool/">http://www.gnu.org/software/libtool/</a>.
  - GNU m4 <a href="http://www.gnu.org/software/m4/">http://www.gnu.org/software/m4/</a>.
  - Perl <a href="http://www.perl.org/">http://www.perl.org/</a>.

\section Installation

From main source directory, type
\code 
  ./configure
  make
  make install
\endcode
The default installation directory for the library is $/usr/local/lib/, and for the Python interface is /usr/local/python/.
These defaults can be changed by using the --prefix flag to ./configure, as in
\code 
  ./configure --prefix=$HOME
\endcode

If installing from the Subversion repository, change to the ariadne/trunk/ directory and type 
\code 
  ./bootstrap
\endcode
Then follow the directions above.

*/