cmake_minimum_required(VERSION 3.5)
project(md_simulation_confined_ions)

set(CMAKE_CXX_STANDARD 14)

add_executable(md_simulation_confined_ions
        NanoconfinementMd.h
        databin.h
        forces.h
        interface.h
        thermostat.h
        vector3d.h
        control.h
        energies.h
        functions.h
        particle.h
        utility.h
        vertex.h
        NanoconfinementMd.cpp
        functions.cpp
        interface.cpp
        main.cpp
        md.cpp
        mdenergies.cpp
        mdforces.cpp)

find_package(GSL REQUIRED)
target_link_libraries(md_simulation_confined_ions GSL::gsl GSL::gslcblas)

find_package(Boost COMPONENTS system filesystem program_options mpi serialization REQUIRED)
target_link_libraries(md_simulation_confined_ions ${Boost_LIBRARIES})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(md_simulation_confined_ions ${MPI_CXX_LIBRARIES})
find_package(OpenMP)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

