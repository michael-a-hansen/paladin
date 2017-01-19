
SET( SOURCES
     compute-and-write-eigenvalues.cpp ) 

set( libs paladin ${BLASLAPACK} ${MPI_LIBRARIES} )

ADD_EXECUTABLE( test-paladin compute-and-write-eigenvalues.cpp )
TARGET_LINK_LIBRARIES( test-paladin ${libs} )

include( compare_spectra.cmake )

set( MPIRANKS 1 2 4 )

foreach( np ${MPIRANKS} )
	compare_spectra( identity_matrix_x1  identity-matrix 10x10.listing1  ${np} exact )
	compare_spectra( identity_matrix_x2  identity-matrix 10x10.listing2  ${np} exact )
	compare_spectra( identity_matrix_x10 identity-matrix 10x10.listing10 ${np} exact )
	compare_spectra( canonical_2x2_x1    canonical-2x2   listing1        ${np} exact )
	compare_spectra( canonical_2x2_x2    canonical-2x2   listing2        ${np} exact )
	compare_spectra( combustion_x1       combustion      listing         ${np} gold  )
	compare_spectra( combustion_x2       combustion      listing2        ${np} gold  )
endforeach( np ${MPIRANKS} )
