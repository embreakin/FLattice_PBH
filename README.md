# FLattice_PBH
## Features

This is a modified version of [FLattice](https://github.com/Axion243/FLattice) (though it appears that FLattice itself has been updated and improved in its independent path since I have started using it for my own use.), written predominantly in C++. This program is focused on solving and simulating the behavior of fields of double inflation type models. For the two inflationary periods of double inflation, the program solves first-order perturbation equations using the Quality-Controlled Runge-Kutta method in the momentum space. For the oscillatory period between the two inflationary periods, the program conducts lattice simulation in real space so that it captures non-linear effects of scalar fields.

Some notable features are...

- Variables are implemented in reduced Planck mass.  
- Calculation of first-order metric perturbation are also implemented in lattice simulation.
- Output of power spectrum, with which one can analyze PBH abundance of the model.
- OpenMP parallelization is implemented to speed up the computation process.
- JSON file is added just in case so that one can utilize it to manage certain data (such as names of output files) if needed. However, it is certainly optional to use this feature and can be ignored. I have attempt to use it to manage parameters of models, but I have found out that the program actually runs faster when one defines the parameters using macros than when setting them using a json file.
- Doxyfile is added so that one can create documents for the program.
- **Quantum vacuum fluctuations** are used for inintial fluctuations of the scalar fields and metric perturbations when the oscillation period starts. The initial conditions are set in momentum space and then Fourier transformed to real space to give the initial values of the fields and their derivatives at each grid point on the lattice. This is a feature that is implemented in the well-known lattice simulation code [LatticeEasy](http://www.felderbooks.com/latticeeasy/). More info about the feature can be found there. 

## Requirements

The following items are necessary and must be installed in your computer to run the program.

- Cmake  
   Cmake is used to control the software compilation process. 
- C++ Compiler  
   The program can be compiled in one of the following three compiliers: 
   - Intel C++ Compiler
   - GCC (g++) 
   - Clang (also AppleClang).  
   
   Note that in order to compile in Clang, one must also install the `libomp` library in order to use OpenMP.
- FFTW  
   The `fftw` library is used for computing discrete Fourier transforms (DFTs). `fftw` files are assumed to be installed in `/opt/fftw/`. If you have installed `fftw` in a different directory, change `include_directories` and `link_directories` in the cmake file.
- Boost  
   The `boost` library is used to manage output files and direcotries.
- Doxygen   
   This must be installed to generate the document for the program. To generate the document,
   
   ```bash
   cd [PATH to the downloaded directory]
   doxygen
   open html/index.html
   ```
- Graphviz  
   If one wants to see diagrams in the document created using Doxygen, one must also install this tool.

## How to Use the Program

Here is a simple explanation on how to use the program.

1. Create `Makefile` by `cmake`

   Go to the directory and create `Makefile` in a `build` directory.

   ```bash
   cd [PATH to the downloaded directory]
   mkdir build
   cd build
   cmake ..
   ```

   `cmake` creates the `Makefile` automatically.

2. Run

   Typing

   ```bash
   make
   ```

   in the `build` directory creates the executable file as `FLattice`. Then, run the code by

   ```bash
   ./FLattice
   ```
3. Selecting Another C++ compiler (Optional)

   If necessary, one can switch to a different C++ compiler. CMake stores the path of the selected compiler inside a variable called
   `CMAKE_CXX_COMPILER`. This can be set for example by using an environment variable:
   ```bash
   export CXX=<compiler name>
   ```
   Compiler names for Intel C++ Compiler, GCC and Clang are `icpc`, `g++` and `clang++`, respectively.  
   Then, you can rerun the program by
   
   ```bash
   cd [PATH to the build directory]
   rm -rf *
   cmake ..
   make
   ```
  
3. Results

   Data such as time, field averages, field variances, etc. will be written in `status.txt`. The field values in each time step are stored as `.vtk` files in the `data` directory when `dim` is 2 or 3. I recommend you use `Paraview.app` to visualize these files.

## Contents

- `main.cpp`

 `main.cpp` is designed to run only the main calculation, which can be devided into three phases. First we have the smooth hybrid inflation phase (1), then the oscillatory phase (2) and lastly the new inflation phase (3). (1) and (3) are both calculated by evolving the linear perturbation equations in wave number space. (2) can be calculated either by evolving the linear perturbation equations just as in (1) and (3), or by switching to real space and running lattice simulations. 
 
- `src`

 Classes and functions used in `main.cpp` are defined in other CPP files in the `src` directory, and the header files for each one of them are stored in the `include` directory. We also have the `lattice_src` directory in the `src` directory, which contains files necessary to perform lattice simulations.

   - `calculation.cpp`

     This file contains functions necessary to calculate linear perturbation equations in wave number space. We have the `Zeromode` class that calculates first only the zeormode of fields. This is necessary to find out the value of the potential minimum and make sure that it is set to zero. After that we use the `Perturbation` class to calculate both the zeromode and the linear perturbation for each wave number.

   - `equations.cpp`

     This file contains functions necessary to actually evolve the linear perturbation equations in wave number space, including the actual potential and its derivatives and so on. 

   - `nr.cpp`

     Functions in this file describes numerical methods used for evolving the linear perturbation equations in wave number space. They are implemented based on `Numerical Recipes`. Hence the file name `nr.cpp`.

   - `parameter.cpp`

     Parameters necessary for the linear perturbation equations are stored in this file. Parameters necessary for lattice simulations are also stored in here.

   - `uc.cpp`

      Three different units are used for wave modes depending on the situation. They are reduced Planck mass, Mpc and the number of the mode counting from the lowest one. This file contains all the functions necessary to convert from one unit to another. They are the unit conversion functions, hence the file name  `uc.cpp`.

   - `utilities.cpp`

      This file contains miscellaneous functions that don't fit into other files for both linear perturbation calculations and lattice simulations. These include functions to output power spectums, set the gaussian distributed perturbations for the lattice, Fourier transformations and output all sorts of detailed data into files.

   - `lattice_src`

      This directory stores the files necessary for running lattice simulations.

      - `lattice.cpp`

        This file corresponds to the main function for lattice simulations, and describes the main flow of the calculation.
      - `lattice_initialize.cpp`
      
         This is where fields are initialized using quantum vacuum fluctuations.

      - `lattice_field.cpp`

        This is where functions related to fields such as scalar potentials, their derivatives, field averages and field variances are implemented. Fields        are also finalized here.

      - `lattice_evol.cpp`

        Functions necessary to calculate time evolution are stored here. 2nd-order and 4th-order simplectic integration scheme ( leap-frog method ) are implemented here, but only the 2nd-order should be trusted for now.

      - `lattice_calc.cpp`

        This file is in charge of calculating values regarding energy density. 

