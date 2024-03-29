cmake_minimum_required(VERSION 2.8)

####
# project settings
####

set(sample_VERSION_MAJOR 1)
set(sample_VERSION_MINOR 0)

set (CMAKE_CXX_COMPILER "/home/users/markert/bin/g++")
set (CMAKE_C_COMPILER "/home/users/markert/bin/gcc")


project(lensingCalculator)

include_directories(include)

####
# sources
####

#set(SOURCES lensingCalculator.cpp myUtilities.cpp inputFunctions.cpp lensFitter.cpp pixelmapFunctions.cpp)
file(GLOB SOURCES "src/lensingCalculator.cpp" "src/inputFunctions.cpp" "src/pixelmapFunctions.cpp" "src/lensFitter.cpp" "src/myUtilities.cpp" "src/lensingClasses.cpp" "src/outputFunctions.cpp")

####
# dependencies
####

find_package(GLAMER NO_MODULE REQUIRED)
include_directories(${GLAMER_INCLUDE_DIRS})

####
# omp
####
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

####
# build
####
message(${GLAMER_LIBRARIES})
message(${GLAMER_INCLUDE_DIRS})
add_executable(lens_fitter ${SOURCES})

target_link_libraries(lens_fitter ${GLAMER_LIBRARIES} "-lcfitsio" "-lCCfits")