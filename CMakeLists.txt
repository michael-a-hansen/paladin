PROJECT( paladin )
cmake_minimum_required( VERSION 3.1.0 )

SET( CMAKE_CXX_STANDARD 11 )

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic" )

ADD_SUBDIRECTORY( src )