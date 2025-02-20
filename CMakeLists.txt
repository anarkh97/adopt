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

# -----------------------------
# External packages required.

# Temporary workaround: Dakota's build tree 
# option is currently not functioning as expected.
set(Dakota_BUILD_TREE
    ON CACHE BOOL "Enable Dakota build tree model")

set(Dakota_DIR 
    "$ENV{DAK_BUILD}/install/lib/cmake/Dakota/" 
    CACHE PATH "Path to the Dakota configuration directory.")
find_package(Dakota REQUIRED PATHS ${Dakota_DIR} NO_DEFAULT_PATH)

if(Dakota_FOUND)
	message(STATUS "Found Dakota: ${Dakota_CONFIG}")
else()
	message(FATAL_ERROR "ERROR: Could not find Dakota ...")
endif()

# -----------------------------
# Version header
ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_SOURCE_DIR}/version.h
           ${CMAKE_SOURCE_DIR}/_version.h
    COMMAND cmake -P
            ${CMAKE_SOURCE_DIR}/version.h)

# -----------------------------
# add include and link directories 
include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${Dakota_INCLUDE_DIRS}")
include_directories("${Dakota_TPL_INCLUDE_DIRS}")

foreach(def ${Dakota_DEFINES})
  add_definitions(-D${def})
endforeach()

# -----------------------------
# add executable
add_executable(adopt
Main.cpp
DakotaEnvironmentWrapper.cpp
AdaptiveJegaOptimizer.cpp
JegaDriver.cpp
JegaEvaluatorCreator.cpp
JegaEvaluator.cpp
)

# -----------------------------
# link to libraries
target_link_libraries(adopt ${Dakota_LIBRARIES})
target_link_libraries(adopt ${Dakota_TPL_LIBRARIES})
add_dependencies(adopt VersionHeader)
