cmake_minimum_required(VERSION 3.25)
project(swapMC)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

#add_compile_options(-lcnpy -lz)

#list(APPEND CMAKE_PREFIX_PATH "/home/rsimon/.local/include/")
#find_package(n REQUIRED)
#include_directories(~/.local/include/arrow/)
#add_subdirectory(~/.local/include/arrow/)
#include(~/.local/include/arrow/)

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
        input/Parameter.h)

#target_include_directories(swapMC PRIVATE "/home/rsimon/.local/include/")
#target_link_directories(swapMC PRIVATE "/home/rsimon/.local/lib/")
#target_link_libraries(swapMC cnpy z)

target_include_directories(swapMC PRIVATE "/home/rsimon/.local/include/")
target_link_directories(swapMC PRIVATE "/home/rsimon/.local/lib/")
target_link_libraries(swapMC cnpy z)