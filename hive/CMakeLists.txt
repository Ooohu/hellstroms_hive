cmake_minimum_required (VERSION 3.3)
project (Hive VERSION 2.0.0 )
include(xgboost/cmake/Utils.cmake)

# The version number.
set (Hive_VERSION_MAJOR 2)
set (Hive_VERSION_MINOR 0)


if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
          set (CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS}")
    endif ()
    else ()
	  set (CMAKE_CXX_STANDARD 14)
endif ()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#This works if your root as compiled with cmake, 
find_package(ROOT REQUIRED COMPONENTS TMVA TreePlayer)# COMPONENTS Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread)
include(${ROOT_USE_FILE})

set (ROOTXTRAS "-lTMVA -lTreePlayer")
execute_process(COMMAND "root-config" "--cflags" "--libs" OUTPUT_VARIABLE rootcflag)
execute_process(COMMAND "root-config" "--incdir" OUTPUT_VARIABLE ROOT_INCDIR)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${ROOTXTRAS} ")

#set (CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
message(STATUS "OUTPUT of CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

include_directories("${PROJECT_BINARY_DIR}")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc/)
include_directories(${ROOT_INCDIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tinyxml/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/xgboost/include/)

add_subdirectory(tinyxml)
add_subdirectory(src)
add_subdirectory(xgboost)

#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/xml/default.xml DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/build/src/)

#CMAKE Error on macOS Mojave (10.14.4), saying  Could NOT find OpenMP_C (missing: OpenMP_C_FLAGS OpenMP_C_LIB_NAMES)? 
# IT IS Clang, mac default compile's fault.
# to make cmake works: https://www.iditect.com/how-to/53889826.html
# install gcc before make: brew install gcc; and use it instead of Clang
# export CC=<your gcc>; see `ls /usr/local/bin/gc*`, dont use Clang.. (check by looking at CC --version)
# similarly CXX=<your g++>;

# if export does not work, us `ln -s`
# see linking command: https://apple.stackexchange.com/questions/99077/how-to-set-gcc-4-8-as-default-gcc-compiler


#CMake Error at /usr/local/Cellar/cmake/3.11.4/share/cmake/Modules/TestBigEndian.cmake:49 (message):
#  no suitable type found
#
# (NO) might need to update xcode? see https://stackoverflow.com/questions/34617452/how-to-update-xcode-from-command-line
# down load xcode, then libraries will be available.
