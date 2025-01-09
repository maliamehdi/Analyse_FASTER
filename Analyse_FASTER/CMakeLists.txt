# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

############ Make your changes here ############

set(name    Analyse)     # equivalent to "name := exampleN01" in G4 GNUmakefiles

file(COPY PID DESTINATION .)
file(COPY Parameter_Files DESTINATION .)

############ What's below should not need to be changed ############

# http://www.cmake.org/cmake/help/cmake_tutorial.html
# http://www.cmake.org/cmake/help/cmake2.6docs.html
# http://www.cmake.org/Wiki/CMake_FAQ

set(mymain ${name}_Main.cxx )
set(mymainlib ${name}_Main.h)
set(myexe  my${name}    )
set(myproj ${name}_   )  # any name as long as it is different from myexe

project(Analyse)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()


#Set C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wno-unused-parameter -Wno-sign-compare -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -lm -ldl -lSpectrum" )



set(CMAKE_POSITION_INDEPENDENT_MODE ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

#----------------------------------------------------------------------------
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# Nota Bene: headers are included so they will show up in IDEs
#
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${ROOTSYS}/include ${INCLUDE_DIR})
add_definitions(${ROOT_CXX_FLAGS})

ROOT_GENERATE_DICTIONARY(RootLibDict *.h LINKDEF ${INCLUDE_DIR}/DictLinkDef.h)
ROOT_LINKER_LIBRARY(RootLib *.cxx RootLibDict.cxx LIBRARIES Core)

#——————————————————————————————————————
#
add_executable(${myexe} ${mymain} ${SOURCE_LIST} ${HEADER_LIST})
target_include_directories(${myexe} PUBLIC include)
target_link_libraries(${myexe} RootLib Spectrum ${ROOT_LIBRARIES} rt)


