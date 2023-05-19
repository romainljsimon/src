cmake_minimum_required(VERSION 3.25)
project(swapMC)

set(CMAKE_CXX_STANDARD 17)

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
        util.h)