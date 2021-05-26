![Image](test.gif "header")


SWESolver
=========
SWESolver is a software written in C++ for the numerical solution of the 2D shallow water
equations. 

The software solves the equations using a finite volume method (FVM) centered on nodes of a 
2D structured quadrilateral mesh. 

### Prerequisites
The only prerequisites for compiling SWESolver are:
- [CMake](https://cmake.org/)
- [Doxygen](https://www.doxygen.nl/index.html) (optional, used for the documentation)

For Debian/Ubuntu they can be downloaded using the following commands on terminal:
```
$ sudo apt-get install cmake
$ sudo apt-get install doxygen
```

### Build
To build SWESolver, move into folder `build` and generate the makefile using CMake:
```
$ cmake ..
```
Then the software can be compiled using the command `make`:
```
$ make -j [jobs]
```
where the option `-j [jobs]` can be used for speeding up the compilation specifying the 
number of jobs.

### Documentation
To generate the documentation (optional), just use the command `doxygen` in the main folder:
```
$ doxygen
```

### Path setup
After installation, the `PATH` to SWESolver can be added to `.bashrc` telling the OS where
to look for the binary. Edit `.bashrc` adding the following line at the bottom:
```
export PATH=/your/path/to/SWESolver/bin:$PATH
```

### Running a simulation
Copy the configuration file `config.cfg` into the simulation folder. Open the terminal and
run the command:
```
$ SWE config.cfg
```

### Test cases
Several test cases for SWESolver are available. They can be downloaded directly from GitHub 
moving into the [TestCases](https://github.com/dellanoce/TestCases) repository or cloned running on terminal:
```
$ git clone git@github.com:dellanoce/TestCases.git
```

### License 
This software is distributed under the GNU GPLv3 license ([COPYING](COPYING)).
