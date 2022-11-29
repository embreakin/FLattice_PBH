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
(I also assume that you are using Mac and Homebrew is already installed).

- Cmake  
   Cmake is used to control the software compilation process. Install it by
   
   ```bash
   brew install cmake
   ```
- C++ Compiler  
   The program can be compiled in one of the following three compiliers: 
   - Intel C++ Compiler
   - GCC (g++) 
   - Clang (also AppleClang).  
   
   Clang is installed in MacOS by default. Intel and GCC are not and must be installed manually.
 
 - OpenMP  
   OpenMP parallelization is implemented to speed up the computation process. For Intel and GCC, OpenMP can be used by default.
   However for Clang, you must install the `libomp` library in order to use OpenMP.
  
   ```bash
   brew install libomp
   ```
- FFTW  
   The `fftw` library is used for computing discrete Fourier transforms (DFTs). Install it by
    ```bash
   brew install fftw
   ```
   `fftw` files are assumed to be installed in `/opt/fftw/`. If you have already installed `fftw` in a different directory, change `include_directories` and `link_directories` in the cmake file.
- Boost  
   The `boost` library is used to manage output files and direcotries.  
   If you are using Intel C++ Compiler or Clang, then 
   ```bash
   brew install boost
   ```
   will install the library. However for GCC, this doesn't work.  
   For GCC, do the following:  
1. First, make sure that you don't have boost installed. If you have installed boost using Homebrew, then make sure you uninstall it by
   ```bash
   brew uninstall --ignore-dependencies boost
   ```
1. Download boost from [Boost Download page](https://www.boost.org/users/download/) and unzip it.
1. Go to where it was downloaded by
   ```bash
   cd ~/Downloads/boost_1_80_0/
   ```
   The number of boost depends on the version you have downloaded (I have downloaded `boost_1_80_0`).
1. Start by bootstrapping the installer:
    ```bash
   ./bootstrap.sh
   ```
1. Once the above is finished, open project-config.jam in a text editor and comment out these lines:
   ```bash
    #if ! clang in [ feature.values <toolset> ]
    #{
    #    using clang ;
    #}
   
   #project : default-build <toolset>clang ;
   ```
1. Now, specify the version of GCC you want to use by adding an extra line to project-config.jam:
   ```bash
   #project : default-build <toolset>clang ;
   
   using gcc : 12 : g++-12 ;
   ```
   Again, the number depends on the gcc version you have installed (Mine is version 12).
1. At this point, you are ready to build and install Boost by GCC:
    ```bash
   sudo ./b2 --toolset=gcc install
   ```
1. You might also need to add the Boost libraries to the dynamic libraries path with:
   ```bash
   export DYLD_LIBRARY_PATH=/usr/local/boost-1.80.0/lib:$DYLD_LIBRARY_PATH
   ```
1. You are now ready to run the program using GCC (detail of this is explained below).
- Doxygen   
   This must be installed to generate the document for the program. Install it by
   ```bash
   brew install doxygen
   ```
   To generate the document,
   
   ```bash
   cd [PATH to the downloaded directory]
   doxygen
   open html/index.html
   ```
   
- Graphviz  
   If one wants to see diagrams in the document created using Doxygen, one must also install this tool. Install it by
   ```bash
   brew install graphviz
   ```
   
## How to Run the Program

Here is a simple explanation on how to run the program.

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

   in the `build` directory creates the executable file as `FLattice_PBH`. Then, run the code by

   ```bash
   ./FLattice_PBH
   ```
3. Selecting Another C++ compiler (Optional)

   If necessary, one can switch to a different C++ compiler (assuming of course that multiple compilers are installed in your computer). CMake stores the path of the selected compiler inside a variable called
   `CMAKE_CXX_COMPILER`. This can be set for example by using an environment variable:
   ```bash
   export CXX=<compiler name>
   ```
   Compiler names for Intel C++ Compiler, GCC and Clang are `icpc`, `g++` and `clang++`, respectively.  
   Check to see if the compiler has been successfully switched by doing
    ```bash
   echo $CXX
   ```
   Then, you can rerun the program by
   
   ```bash
   cd [PATH to the build directory]
   rm -rf *
   cmake ..
   make
   ./FLattice_PBH
   ```
   
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

   - `parameters.cpp`

     Parameters necessary for the linear perturbation equations and lattice simulations are stored in this file. 

   - `uc.cpp`

      Three different units are used for wave modes depending on the situation. They are reduced Planck mass, Mpc and the number of the mode counting from the lowest one. This file contains all the functions necessary to convert from one unit to another. They are the unit conversion functions, hence the file name  `uc.cpp`.

   - `utilities.cpp`

      This file primarily contains functions to output data.

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
 
## Results
   - `<parameter set name>`

      All the files will be stored in a single directory. It is recommended that the directory is named after the parameter set of the model that one wants to analyze.

      - `<parameter set name>_unpWMAP5.txt`

        Data regarding the evolution of zero mode of fields are stored here.
        
      - `<parameter set name>_spectrum_bfosc.txt`
      
         Power spectrum of curvature purerbation right before oscillation period starts is stored here.

      - `<parameter set name>_spectrum.txt` 

         The final power spectrum of curvature purerbation is stored here.

      - `<parameter set name>_kAnalyze_nonlattice`

         This directory contains data regarding the evolution of zero mode + perturbations when only the linear perturbation is calculated and lattice simulation is turned off.

      - `<parameter set name>_kAnalyze_lattice`

         This directory contains data regarding the evolution of zero mode + perturbations when lattice simulation is conducted. 
        
      - `<parameter set name>_status.txt`

         Data during lattice simulation are written in here.  
        
      - `<parameter set name>_energy`

         The energy density values for each time step during lattice simulation are stored as `.vtk` files in this directory when `dim` is 2 or 3. `Paraview.app` is necessary to visualize these files. 
        
      - `<parameter set name>_field`

         The field values for each time step during lattice simulation are stored as `.vtk` files in this directory when `dim` is 2 or 3. `Paraview.app` is necessary to visualize these files.
