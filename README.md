This is the *stable* release of the Ariadne C++ library.

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
Install the Command Line Developer Tools (will also be asked when installing Homebrew) from the Apple Store
Install Homebrew from http://brew.sh/
Homebrew packages required: `cmake boost gtk cairo`
No Buddy package is offered, you need to compile the library from https://sourceforge.net/projects/buddy/
Download and extract the Buddy package, then from the extracted directory:
```
$ ./configure
$ make
$ make install
```