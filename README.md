# FLattice_PBH
## Features

This is a modified version of [FLattice](https://github.com/Axion243/FLattice) (though it seems that FLattice itself has been updated and improved in its independent path since I have started using it for my own use.) This program is focused on solving and simulating double inflation type models. For the two inflationary periods of double inflation, the program solves first-order perturbation equations using the Runge-Kutta method. For the oscillatory period between the two inflationary periods, the program conducts lattice simulation so that it captures non-linear effects of scalar fields.

Some notable features are...

- Calculation of first-order gravitational perturbation also for lattice simulation.
- Output of power spectrum, which leads to the ability to analyze PBH abundance.
- OpenMP parallelization is implemented.
- Variables are implemented in reduced Planck mass.  
- **Quantum vacuum fluctuations** are used for inintial fluctuations of the scalar fields when the oscilatory phase starts. The initial conditions are set in momentum space and then Fourier transformed to give the initial values of the fields and their derivatives at each grid point. This is a feature that is included in the well-known lattice simulation code [LatticeEasy](http://www.felderbooks.com/latticeeasy/) and more info about the feature can be found there. 

## Requirements

The following items are required to run the code.
- cmake
- Intel compiler
- fftw

fftw files are assumed to be installed in  `/opt/fftw/`. If you have installed fftw in a different directory, change `include_directories` and `link_directories` in the cmake file.

## How to Use

Here is a simple explanation on how to use the code.

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

         Three different units are used for wave modes depending on the situation. They are reduced Planck mass, Mpc and the number of the mode counting from the lowest one. This file contains all the functions necessary to convert from one unit to another different unit. They are the unit conversion functions, hence the file name  `uc.cpp`.

      - `utilities.cpp`

         This file contains miscellaneous functions that doesn't really fit into other files for both linear perturbation calculations and lattice simulations. These include functions to output powerspectum, set the gaussian distributed perturbations for the lattice, Fourier transformations and output all sorts of detailed data into files.

      - `lattice_src`

         This directory stores the files necessary for running lattice simulations.

         - `lattice.cpp`

           This file corresponds to the main function for lattice simulations, and describes the main flow of the calculation.

         - `lattice_field.cpp`

           This is where the fields are initialized using quantum vacuum fluctuations, and also finalized. Functions related to fields such as scalar potentials, their derivatives, field averages and field variances are also implemented here.

         - `lattice_evol.cpp`

           Functions necessary to calculate time evolution are stored here. 2nd-order and 4th-order simplectic integration scheme ( leap-frog method ) are implemented here, but only the 2nd-order should be trusted for now.

         - `lattice_calc.cpp`

           This file is in charge of calculating values regarding energy density. 

