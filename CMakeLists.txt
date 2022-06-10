# Check CMake Version 
cmake_minimum_required(VERSION 3.19.0 FATAL_ERROR)

#
#---------------------------------------------------------------------
# Project Setup
#---------------------------------------------------------------------
#
project("CGIL" VERSION 1.0.0.0 LANGUAGES CXX)
#set(PROJECT_VERSION_MAJOR 0)
#set(PROJECT_VERSION_MINOR 0)
#set(PROJECT_VERSION_PATCH 0)
#set(PROJECT_VERSION_TWEAK 1)

#---------------------------------------------------------------------
# Set location of *.cmake modules
#---------------------------------------------------------------------
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#---------------------------------------------------------------------
# User Configure Build Process
#---------------------------------------------------------------------
include(GNUInstallDirs)

set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_STATUS_MAKEFILE TRUE)
set(CMAKE_COLOR_MAKEFILE TRUE)

option(CMAKE_BUILD_TYPE "Release")
option(CGIL_BUILD_TESTING "Build Testing" OFF)

#
#---------------------------------------------------------------------
# Force build in seperate directory
#---------------------------------------------------------------------
#
include(InSourceBuild)

#
#---------------------------------------------------------------------
# Directory Report
#---------------------------------------------------------------------
#
message(STATUS "")
message(STATUS "--------------------- Directory Report -----------------------")
message(STATUS " ") # CMAKE directories are from top level CMakeLists.txt
message(STATUS "Top Level Directories:")
message(STATUS "CMAKE_SOURCE_DIR  = ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_BINARY_DIR  = ${CMAKE_BINARY_DIR}")
message(STATUS " ") # PROJECT directories are for recent project call
message(STATUS "Project Level Directories:")
message(STATUS "PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")


#---------------------------------------------------------------------
# Locations of Installation & Report
#---------------------------------------------------------------------

# Source locations in source tree
set(MY_PROJECT_PREFIX ${PROJECT_SOURCE_DIR})
set(MY_PROJECT_INCDIR ${MY_PROJECT_PREFIX}/include)
set(MY_PROJECT_SRCDIR ${MY_PROJECT_PREFIX}/src)
set(MY_PROJECT_TSTDIR ${MY_PROJECT_PREFIX}/test)
set(MY_PROJECT_DOCDIR ${MY_PROJECT_PREFIX}/docs)

# Where to place all libraries
set(MY_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
set(MY_INSTALL_BINDIR ${MY_INSTALL_PREFIX}/bin)
set(MY_INSTALL_INCDIR ${MY_INSTALL_PREFIX}/include)
set(MY_INSTALL_SRCDIR ${MY_INSTALL_PREFIX}/src)
set(MY_INSTALL_TSTDIR ${MY_INSTALL_PREFIX}/test)
set(MY_INSTALL_DOCDIR ${MY_INSTALL_PREFIX}/docs)

message(STATUS " ")
message(STATUS "Using Directories:")
message(STATUS "--- Original Locations ---")
message(STATUS "MY_PROJECT_PREFIX = ${MY_PROJECT_PREFIX}")
message(STATUS "MY_PROJECT_INCDIR = ${MY_PROJECT_INCDIR}")
message(STATUS "MY_PROJECT_SRCDIR = ${MY_PROJECT_SRCDIR}")
message(STATUS "MY_PROJECT_TSTDIR = ${MY_PROJECT_TSTDIR}")
message(STATUS "MY_PROJECT_DOCDIR = ${MY_PROJECT_DOCDIR}")
message(STATUS " ")
message(STATUS "--- Installation Locations ---")
message(STATUS "MY_INSTALL_PREFIX = ${MY_INSTALL_PREFIX}")
message(STATUS "MY_INSTALL_BINDIR = ${MY_INSTALL_BINDIR}")
message(STATUS "MY_INSTALL_INCDIR = ${MY_INSTALL_INCDIR}")
message(STATUS "MY_INSTALL_TSTDIR = ${MY_INSTALL_TSTDIR}")
message(STATUS "MY_INSTALL_DOCDIR = ${MY_INSTALL_DOCDIR}")

#---------------------------------------------------------------------
# Detect Library includes, flags, etc.
#---------------------------------------------------------------------

message("")
message("================================================================")
message("                    Searching for Libraries                     ")
message("================================================================")

find_package(_Boost REQUIRED)

# find_package(_CGAL)
# if(NOT CGAL_FOUND)
# message(STATUS "CGAL will be Downloaded & Installed")
add_subdirectory(third_party)
# endif()

#
#---------------------------------------------------------------------
# Compiler Feature Detection
#---------------------------------------------------------------------
# Hasn't proven to be very usefull. I find testing for the exact 
# feature of interest in the following section much better.
# Partialy because each compiler names features differently or 
# some not at all.
#
message(STATUS "")
message(STATUS "--------------- Compiler Feature Detection -----------------")
message(STATUS "")
message(STATUS "CMAKE_CXX_COMPILER         = ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER_ID      = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")


#---------------------------------------------------------------------
# Set Project Level Compiler Options
#---------------------------------------------------------------------
# ONLY those options which impact level of compiler errors
add_compile_options("-Werror")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  add_compile_options("-stdlib=libc++")
  add_link_options("-stdlib=libc++")
endif()

#---------------------------------------------------------------------
# 
#---------------------------------------------------------------------

# Build Libraries & Applications
add_subdirectory(src)

# Build Tests for Libraries
include(CTest)
if(CGIL_BUILD_TESTING)
	add_subdirectory(test)
endif()
