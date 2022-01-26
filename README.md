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

 `main.cpp` is designed to run only the main calculation. Classes and functions used in `main.cpp` are defined in the other five CPP files in the `lib` directory, and the header files for each one of them are stored in the `include` directory.
The five CPP files are the following.

- `parameter.cpp`

  Parameters necessary for the simulation are stored in this file.

- `field.cpp`

  This is where the fields are initialized using quantum vacuum fluctuations. Other configurations such as the scalar potential, effective mass, etc. are set in `field.hpp`.

- `evolution.cpp`

  All classes necessary to calculate the time evolution are in `evolution.cpp`. 2nd-order and 4th-order simplectic integration scheme ( leap-frog method ) are implemented here, but only the 2nd-order should be trusted.
  
- `calculation.cpp`
  
  All ouput values such as field averages, field variances, energy densities etc. are calculated in this file.   
  
- `utilities.cpp`
  
  Miscellaneous functions necessary to initialize the fields are implemented here. Functions for the final output of data are also written here.  
