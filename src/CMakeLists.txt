set( HEADERS
     command-line-parser.hpp
     lapack-wrapper.hpp
     mpi-util.hpp
     paladin.hpp
     spectrum-util.hpp
     square-matrix.hpp
     types.hpp
     header.hpp
     make-unique.hpp
     )

# create dummy/empty file to circumvent restriction on add_library
set( paladin_sources compute-write-with-timings.cpp )
add_library( paladin STATIC ${paladin_sources} )
target_link_libraries( paladin ${BLASLAPACK} ${MPI_LIBRARIES} )

install( FILES
         ${HEADERS}
         DESTINATION include/paladin
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
         )
         
install( TARGETS
         paladin
         EXPORT Paladin
         DESTINATION lib/paladin/
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
         )
install( EXPORT Paladin DESTINATION lib/paladin )


SET( SOURCES
     compute-write-with-timings.cpp ) 

set( libs paladin ${BLASLAPACK} ${MPI_LIBRARIES} )

ADD_EXECUTABLE( exec-paladin compute-write-with-timings.cpp )
TARGET_LINK_LIBRARIES( exec-paladin ${libs} )