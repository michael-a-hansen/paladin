function( compare_spectra testName testDir matrixListingName numProc referenceType )
	
	file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${testDir} )
	file( GLOB testfiles ${CMAKE_CURRENT_SOURCE_DIR}/${testDir}/* )
	file( REMOVE_RECURSE ${CMAKE_CURRENT_SOURCE_DIR}/${testDir}/*matrix.* )
	file( COPY ${testfiles} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${testDir}/ )

	add_test( NAME ${testName}_compute_np${numProc}
	          COMMAND mpirun -np ${numProc} test-paladin
	          --listing=${CMAKE_CURRENT_BINARY_DIR}/${testDir}/${matrixListingName}
	          --rootdir=${CMAKE_CURRENT_BINARY_DIR}/${testDir}
	          --left
	          --right )

  	FILE( READ ${CMAKE_CURRENT_BINARY_DIR}/${testDir}/${matrixListingName} matrixFiles )
  	STRING( REGEX REPLACE ";" "\\\\;" matrixFiles "${matrixFiles}" )
  	STRING( REGEX REPLACE "\n" ";" matrixFiles "${matrixFiles}" )
 
  	set( INC 0 )
	FOREACH( matrixFile ${matrixFiles} )
		MATH(EXPR INC "${INC}+1")

		set( computed_spectrum ${CMAKE_CURRENT_BINARY_DIR}/${testDir}/${matrixFile}.spectrum )
		set( reference_spectrum ${CMAKE_CURRENT_BINARY_DIR}/${testDir}/${matrixFile}.spectrum.${referenceType} )
		set( computed_leftvecs ${CMAKE_CURRENT_BINARY_DIR}/${testDir}/${matrixFile}.leftvecs )
		set( reference_leftvecs ${CMAKE_CURRENT_BINARY_DIR}/${testDir}/${matrixFile}.leftvecs.${referenceType} )		
		set( computed_rightvecs ${CMAKE_CURRENT_BINARY_DIR}/${testDir}/${matrixFile}.rightvecs )
		set( reference_rightvecs ${CMAKE_CURRENT_BINARY_DIR}/${testDir}/${matrixFile}.rightvecs.${referenceType} )
		
    	add_test( NAME ${testName}_vs_${referenceType}_spectrum_np${numProc}_matrix${INC}
    	          COMMAND ${CMAKE_COMMAND} -E compare_files ${computed_spectrum} ${reference_spectrum} )
    	add_test( NAME ${testName}_vs_${referenceType}_leftvecs_np${numProc}_matrix${INC}
    	          COMMAND ${CMAKE_COMMAND} -E compare_files ${computed_leftvecs} ${reference_leftvecs} )
    	add_test( NAME ${testName}_vs_${referenceType}_rightvecs_np${numProc}_matrix${INC}
    	          COMMAND ${CMAKE_COMMAND} -E compare_files ${computed_rightvecs} ${reference_rightvecs} )
    	          
		set_property( TEST ${testName}_vs_${referenceType}_spectrum_np${numProc}_matrix${INC} APPEND PROPERTY DEPENDS ${testName}_compute_np${numProc} )
		set_property( TEST ${testName}_vs_${referenceType}_leftvecs_np${numProc}_matrix${INC} APPEND PROPERTY DEPENDS ${testName}_compute_np${numProc} )
		set_property( TEST ${testName}_vs_${referenceType}_rightvecs_np${numProc}_matrix${INC} APPEND PROPERTY DEPENDS ${testName}_compute_np${numProc} )
	ENDFOREACH( matrixFile ${matrixFiles} )

endfunction()