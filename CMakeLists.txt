cmake_minimum_required ( VERSION 3.14 )

project( FLattice_PBH CXX )

# using Intel C++
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
# -parallel is for automatic parallelization
# -qopenmp is for OpenMP parallelization
set(CMAKE_CXX_FLAGS "-ipo -xhost -no-prec-div -qopenmp -parallel")


# using GCC
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
set ( CMAKE_CXX_FLAGS " -Ofast -Wall  -fopenmp" ) #


# using Clang
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

#Get the path where libomp is installed by executring the following:
execute_process(COMMAND brew --prefix libomp
                OUTPUT_VARIABLE OpenMP_HOME
                OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "OpenMP Root : ${OpenMP_HOME}")

# Check if OpenMP is installed
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY "${OpenMP_HOME}/lib/")

# The FORCE option is appended since it is a cash variable
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -Wno-unused-command-line-argument -I${OpenMP_HOME}/include -lomp -L${OpenMP_omp_LIBRARY}" CACHE STRING "" FORCE)


# Look for the libomp package.
find_package(OpenMP)
message(STATUS "OpenMP_CXX_Found : ${OpenMP_CXX_FOUND}")

set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -Ofast ")

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

