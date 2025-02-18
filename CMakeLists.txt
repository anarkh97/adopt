cmake_minimum_required(VERSION 3.24)

# use gcc and g++ to compile
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)

# set the project name
project(ADOPT VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install")

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
  CACHE STRING "Set the RPATH in Dakota executables and libraries")

# compiler flags (turn on most of the warnings, except unused variables and parameters
list(APPEND MY_CXX_FLAGS "-Wall -Wextra -Wno-reorder")
list(APPEND MY_CXX_FLAGS "-Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter")
list(APPEND MY_CXX_FLAGS "-Wno-unused-but-set-variable -Wno-dev")
list(APPEND MY_CXX_FLAGS "-Wdeprecated-declarations")

set(CMAKE_CXX_FLAGS:STRING ${MY_CXX_FLAGS})

set(CMAKE_CXX_FLAGS_DEBUG "-g3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# -----------------------------
# for version control
# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
add_custom_target(VersionHeader
  ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/version.cmake
)

# Variable used by dakota subdirectories to export libraries
set(ExportTarget MyLibraryTargets CACHE STRING 
  "Export set name for ${CMAKE_PROJECT_NAME}")

#---------------------------------------------------------------------------------------------------------
# External packages required.
#find_package(Boost 1.71 COMPONENTS filesystem REQUIRED)

set(Dakota_CONFIG_DIR 
    "$ENV{DAK_BUILD}/install/lib/cmake/Dakota/" 
    CACHE PATH "Path to the Dakota configuration directory.")

# build JEGA libraries
# could not get this to work -- hard coding Targets.cmake
##add_subdirectory(JEGA)
##install(EXPORT ${ExportTarget}
##  DESTINATION lib/cmake/MyLibrary)
##include(${CMAKE_SOURCE_DIR}/install/lib/cmake/MyLibrary/MyLibraryTargets.cmake)

#file(GLOB DAK_LIBS "$ENV{DAK_BUILD}/install/lib/lib*.so")
#file(GLOB DAK_TPL_LIBS "$ENV{DAK_BUILD}/install/bin/lib*.so")
#include(${CMAKE_SOURCE_DIR}/MyLibraryTargets.cmake)

#---------------------------------------------------------------------------------------------------------

ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_SOURCE_DIR}/version.h
           ${CMAKE_SOURCE_DIR}/_version.h
    COMMAND cmake -P
            ${CMAKE_SOURCE_DIR}/version.h)

# link to dakota libraries
link_directories($ENV{DAK_BUILD}/install/lib)

# -----------------------------
# add include and link directories 
include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${Dakota_INCLUDE_DIRS}")
#include_directories("$ENV{DAK_BUILD}/install/include")
#include_directories("$ENV{DAK_SRC}/packages/external/JEGA/include")
#include_directories("$ENV{DAK_SRC}/packages/external/JEGA/eddy")

add_executable(adopt
Main.cpp
AdaptiveJegaOptimizer.cpp
JegaDriver.cpp
JegaEvaluatorCreator.cpp
JegaEvaluator.cpp
)

# link to libraries
target_link_libraries(adopt ${Dakota_LIBRARIES})
#target_link_libraries(adopt ${DAK_LIBS} ${DAK_TPL_LIBS})
#target_link_libraries(adopt jega utilities jega_fe moga soga eutils)
#target_link_libraries(adopt dakota_src)
#target_link_libraries(adopt teuchoscomm teuchoscore teuchosnumerics teuchosparameterlist teuchosparser teuchosremainder)
#target_link_libraries(adopt pebbl pecos_src pecos_util nomad)
#target_link_libraries(adopt ${DAK_TPL_LIBS})
## unneccessary but required
#target_link_libraries(adopt dakota_dll_api dakota_src_fortran dakota_util ddace dfftpack dream interfaces ncsuopt nidr optpp
#tinyxml utilib) 

add_dependencies(adopt VersionHeader)
