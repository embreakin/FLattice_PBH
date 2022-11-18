cmake_minimum_required ( VERSION 3.14 )

project( FLattice CXX )

if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
# using Intel C++

# -parallel is for automatic parallelization
# -qopenmp is for OpenMP parallelization
set ( CMAKE_CXX_FLAGS "-ipo -xhost -no-prec-div -qopenmp -parallel" )

elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # using GCC
set ( CMAKE_CXX_FLAGS " -Ofast   -fopenmp" ) #-Wall

elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")


else ()

endif()

set ( CMAKE_CXX_STANDARD 17 )

# Add include files (-I option)
include_directories ( ${PROJECT_SOURCE_DIR}/include )
include_directories ( ${PROJECT_SOURCE_DIR}/include/lattice_include )
include_directories ( /opt/fftw/include )

# Specify the Library directory (-L option)
link_directories ( /opt/fftw/lib )


# Add executing files
file ( GLOB src_codes ${PROJECT_SOURCE_DIR}/src/*.cpp )
file ( GLOB lattice_src_codes ${PROJECT_SOURCE_DIR}/src/lattice_src/*.cpp )
#add_library(${Boost_LIBRARIES} SHARED main.cpp ${lib_codes})
add_executable ( ${PROJECT_NAME} main.cpp ${src_codes} ${lattice_src_codes} )

# Boost
set(boost_min_ver  1.69.0)
set(boost_libs filesystem system)
find_package(Boost ${boost_min_ver})

if(Boost_FOUND)
find_package(Boost ${boost_min_ver} COMPONENTS ${boost_libs})
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
endif()

# Link external libraries
target_link_libraries ( ${PROJECT_NAME} fftw3 )


# Link project "library" when compile FLattice
# target_link_libraries ( FLattice library_code )

# Add sub-directory
# add_subdirectory ( lib )

