cmake_minimum_required(VERSION 3.25)
project(swapMC)

set(CMAKE_CXX_STANDARD 17)

#add_compile_options(-O3 -flto)
include_directories(.)

add_executable(swapMC
        energy.cpp
        energy.h
        main.cpp
        MonteCarlo.cpp
        MonteCarlo.h
        random.cpp
        random.h
        readSaveFile.cpp
        readSaveFile.h
        types.h
        unittests.cpp
        unittests.h
        util.cpp
        util.h
        pressure.cpp
        pressure.h
        INPUT/Parameter.h POTENTIALS/PairPotentials.h POTENTIALS/BondPotentials.h POTENTIALS/PairPotentials.cpp POTENTIALS/BondPotentials.cpp NEIGHBORS/Neighbors.cpp PARTICLES/Particles.h PARTICLES/Particles.cpp)