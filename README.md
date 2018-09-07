# Ariadne 
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Build Status](https://travis-ci.org/ariadne-cps/release-1.0.svg?branch=master)](https://travis-ci.org/ariadne-cps/release-1.0) [![codecov](https://codecov.io/gh/ariadne-cps/release-1.0/branch/master/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/release-1.0)

This is the *stable* release of the Ariadne C++ library.

### Installation ###

These installation instructions have been tested on Ubuntu 16.04 and macOS 10.12.

For the Ubuntu installation, we will refer to packages available on Aptitude. The macOS installation instead will assume you are using the Brew package manager.

The build system is CMake. The compiler we tested the library under Ubuntu is g++, while for macOS is clang. To build the library in a clean way, it is preferable that you set up a build subdirectory:

```
$ mkdir build
$ cd build
```

#### Dependencies

The library dependencies of Ariadne are the following:

##### Ubuntu
Aptitude packages required: `git cmake libboost-system-dev libboost-serialization-dev libboost-thread-dev libgtk2.0-dev libcairo2-dev libbdd-dev`

##### macOS
1. Install the Command Line Developer Tools (will also be asked when installing Homebrew) from the Apple Store

2. Install Homebrew from http://brew.sh/

    Homebrew packages required: `cmake boost gtk cairo`

3. No Buddy package is offered, you need to compile the library from https://sourceforge.net/projects/buddy/

    Download and extract the Buddy package, then from the extracted directory:
    
    ```
$ ./configure
$ make
$ make install
    ```

Optionally, if you want to build the documentation, you need Doxygen and a working Latex distribution (including the Math packages).

#### Building

Then you can prepare the build environment:

```
$ cmake ..
```

At this point, if no error arises, you can build the library itself:

```
$ make
```

Optionally, you can also run the test suite for the library:

```
$ make test
```

where no error should appear.

If you want to build the documentation, you have to issue the following:

```
$ make doc
```

### Installing globally

To install the library globally, you must do
```
$ make install
```

To find the installed library under Ubuntu, you may need to set the `LD_LIBRARY_PATH` variable in the .bashrc file:
```
export LD_LIBRARY_PATH=/usr/local/lib
```

To test the correct installation, you may build the `tutorial/` example, which relies on Ariadne being present within system library and include directories.
