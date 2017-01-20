This is the *stable* release of the Ariadne C++ library.

We are in the process of producing a first official release. As such, documentation is not completed yet. In the end, a tutorial on the definition of hybrid automata and their analysis/verification will be available.

### Installation ###

These installation instructions have been tested on Ubuntu 16.04 and OSX 10.11.

For the Ubuntu installation, we will refer to packages available on Aptitude. The OSX installation instead will assume you are using the Brew package manager.

The build system is CMake. The compiler we tested the library under Ubuntu is g++, while for OSX is clang. To build the library in a clean way, it is preferable that you set up a build subdirectory:

```
$ mkdir build
$ cd build
```

#### Dependencies

The library dependencies of ARIADNE are the following:

##### Ubuntu
Aptitude packages required: `git cmake libboost-system-dev libboost-serialization-dev libboost-thread-dev libgtk2.0-dev libcairo2-dev libbdd-dev`

##### OSX
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

To find the installed library under Ubuntu, you may need to set the LD_LIBRARY_PATH in the .bashrc file:
```
export LD_LIBRARY_PATH=/usr/local/lib
```