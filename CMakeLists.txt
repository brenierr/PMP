cmake_minimum_required(VERSION 2.8)
project(CPF C CXX)

# Prefix
if (NOT PREFIX)
    set(PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif (NOT PREFIX)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE STRING "Installation prefix." FORCE)

if (CXX)
    set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()


# Detect MPI
set(MPI_LIBS)
find_package(MPI REQUIRED)
if (MPI_C_FOUND AND MPI_CXX_FOUND)
    include_directories(${MPI_C_INCLUDE_PATH})
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set (MPI_LIBS ${MPI_LIBS} ${MPI_C_LIBRARIES})
    set (MPI_LIBS ${MPI_LIBS} ${MPI_CXX_LIBRARIES})
else (MPI_C_FOUND AND MPI_CXX_FOUND)
    message(STATUS "MPI not found.")
endif (MPI_C_FOUND AND MPI_CXX_FOUND)

include_directories(${PROJECT_SOURCE_DIR})
# Detect Boost
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS program_options REQUIRED)
message(STATUS "Boost program_options component: ${Boost_PROGRAM_OPTIONS_LIBRARY}")
FIND_PACKAGE(Boost COMPONENTS serialization REQUIRED)
MESSAGE(STATUS "Boost serialization : ${Boost_SERIALIZATION_LIBRARY}")
include_directories(${Boost_INCLUDE_DIR})
set(LIBS ${LIBS} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_SERIALIZATION_LIBRARY})

# Detect PNL
find_package(Pnl REQUIRED)
set(LIBS ${LIBS} ${PNL_LIBRARIES})
include_directories(${PNL_INCLUDE_DIRS})
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
endif()

# Detect OPM
find_package(OpenMP QUIET)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS ${CMAKE_CXX_CFLAGS} ${OpenMP_CXX_FLAGS})
    message(STATUS "OpenMP found.")
else (OPENMP_FOUND)
    message(STATUS "OpenMP not found.")
endif (OPENMP_FOUND)

# Add math lib
set(LIBS ${LIBS} m)

# Add math lib
set(MPI_LIBS ${MPI_LIBS} m)

	add_executable(mc-pricer src/main.cpp src/Asian.cpp  src/Basket.cpp  src/Performance.hpp src/Performance.cpp src/BlackScholesModel.cpp  src/MonteCarlo.cpp  src/Option.cpp  src/parser.cpp src/Asian.hpp  src/Basket.hpp  src/BlackScholesModel.hpp src/MonteCarlo.hpp  src/Option.hpp  src/parser.hpp)
        target_link_libraries(mc-pricer ${MPI_LIBS} ${LIBS})

# Print compiler flags
get_directory_property(DirDefs COMPILE_DEFINITIONS)
message(STATUS "COMPILE_DEFINITIONS = ${DirDefs}")
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}" )
    message(STATUS "C Flags: ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}" )
else()
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}" )
    message(STATUS "C Flags: ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}" )
endif()
