PROJECT( paladin CXX )
cmake_minimum_required( VERSION 3.1.0 )

set( MAJOR_VERSION "1" )
set( MINOR_VERSION "0" )
set( PATCH_VERSION "0" )


# warn about building in source
# ------------------------------------------------
if( CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR )
  message( WARNING "\nIt appears that you are building directly in the source directory.\n"
                   "This is strongly discouraged.\n"
                   "Suggestion:\n\tmkdir build; cd build\n"
                   "Then run cmake from the 'build' directory" ) 
endif( CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR )
# ------------------------------------------------



# default to release builds
# ------------------------------------------------
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE "Release" CACHE STRING
       "Choose the type of build, options are: Debug Release" FORCE )
endif( NOT CMAKE_BUILD_TYPE )
# ------------------------------------------------



# Set the location for the installation.
# ------------------------------------------------
if( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  set( CMAKE_INSTALL_PREFIX
       ${CMAKE_CURRENT_BINARY_DIR}
       CACHE PATH "Installation directory" FORCE )
endif()
# ------------------------------------------------
  


# git & version information
# ------------------------------------------------
set( PALADIN_REPO_DATE "\"DATE NOT FOUND\"" )
set( PALADIN_REPO_HASH "\"HASH NOT FOUND\"" )
find_package( Git )
if( GIT_FOUND )
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 "--pretty=format:\"%H\""
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    RESULT_VARIABLE RESULT
    OUTPUT_VARIABLE PALADIN_REPO_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if( NOT ${RESULT} EQUAL 0 )
    set( PALADIN_REPO_HASH "\"HASH NOT FOUND\"" )
  endif()
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log -1 "--pretty=format:\"%cd\""
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE PALADIN_REPO_DATE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
   if( NOT ${RESULT} EQUAL 0 )
     set( PALADIN_REPO_DATE "\"DATE NOT FOUND\"" )
   endif()
endif( GIT_FOUND )

configure_file(
  ${PROJECT_SOURCE_DIR}/src/paladin-configure.h.in
  ${PROJECT_BINARY_DIR}/config/paladin-configure.h
  @ONLY
  )
include_directories( ${PROJECT_BINARY_DIR}/config )
# ------------------------------------------------



# MPI
# ------------------------------------------------
FIND_PACKAGE( MPI REQUIRED )
# ------------------------------------------------



# LAPACK & BLAS
# ------------------------------------------------
IF( APPLE )
	MESSAGE( "Searching for framework Accelerate..." )
	FIND_LIBRARY( BLASLAPACK Accelerate )
	IF( NOT BLASLAPACK )
		MESSAGE( FATAL_ERROR "Mac OS X: framework Accelerate not found!" )
    ELSE()
        MESSAGE( "-- Found framework Accelerate!" )
	ENDIF()
ELSE()
    MESSAGE( "Searching for BLAS and LAPACK..." )
    FIND_PACKAGE( BLAS REQUIRED )
    FIND_PACKAGE( LAPACK REQUIRED )
    SET( BLASLAPACK ${LAPACK_LIBRARIES} )
    IF( NOT BLASLAPACK )
        MESSAGE( FATAL_ERROR "BLAS and LAPACK not found!" )
    ELSE()
        MESSAGE( "-- Found BLAS and LAPACK!" )
    ENDIF()
ENDIF()
# ------------------------------------------------



# C++11
# ------------------------------------------------
SET( CMAKE_CXX_STANDARD 11 )
SET( CMAKE_CXX_STANDARD_REQUIRED ON )
# ------------------------------------------------



# build flags
# ------------------------------------------------
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic" )
# ------------------------------------------------



# includes
# ------------------------------------------------
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/src
                     ${MPI_INCLUDE_PATH} )
# ------------------------------------------------



# source code
# ------------------------------------------------
ADD_SUBDIRECTORY( src )
# ------------------------------------------------



# testing
# ------------------------------------------------
option( ENABLE_TESTS "Enable tests" ON )
if( ENABLE_TESTS )
  enable_testing()
  add_subdirectory( test )
endif( ENABLE_TESTS )
# ------------------------------------------------



message( STATUS "Files will be installed to: " ${CMAKE_INSTALL_PREFIX} )


