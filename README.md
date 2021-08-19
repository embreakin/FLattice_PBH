# FLattice

FLattice is the Fast Lattice code to caluclate the scalar feild evolution ( still in the experimental stage ).

The following items are necessary to run the code ( or fix the relevant parts ).

- cmake
- Intel compiler
- fftw

I assume that FFTW files are installed in  `/opt/fftw/`. If you installed FFTW in other directory, change `include_directories` and `link_directories` in the cmake file.

## How to use

Here, I introduce a simple procedure to use FLattice.

1. Create make file by `cmake`

   Move to the download directory and create `Makefile` in a `build` directory.

   ```bash
   cd /PATH/to/the/FLattice-master
   mkdir build
   cd build
   cmake ..
   ```

   `cmake` create the `Makefile` automatically.

2. Run

   Typing

   ```bash
   make
   ```

   in the `build` directory creates the executable file as `FLattice`. Then, run the code by

   ```bash
   ./FLattice
   ```

3. Check the result

   The data such as time, field averages, field variances, etc.  are written in `status.txt`. The field values in each time step are stored as `.vtk` files in `data` directory when `dim` is 2 or 3. I recommend you to use `Paraview.app` to visualize these files.

## Basic Concept

 `main.cpp` is designed to run only the main calculation. Classes or functions used in `main.cpp` are defined in the other five fiies in the `lib` directory.

- `parameter.cpp`

  Simulation parameters are stored in this file.

- `field.cpp`

  We arrange the field array in `field.cpp`. The class `Field` includes all functions involving the field arrangement, such as seting initial condition, potential, etc.

- `evolutioin.cpp`

  All Classes to calculate the time evolutioin are included in `evolution.cpp`. I still only implement the 2nd-order or 4th-order simplectic integration scheme ( leap-frog method ), but I will add other integration scheme such as 4th-order Runge-Kutta method.