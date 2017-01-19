function( compare_spectra testName testDir matrixListingName numProc referenceType )
	
	file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${testDir} )
	file( GLOB testfiles ${CMAKE_CURRENT_SOURCE_DIR}/${testDir}/* )
	file( REMOVE_RECURSE ${CMAKE_CURRENT_SOURCE_DIR}/${testDir}/*matrix.* )
	file( COPY ${testfiles} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${testDir}/ )

	add_test( NAME ${testName}_compute_np${numProc}
	          COMMAND mpirun -np ${numProc} cawe
	          --listing=${CMAKE_CURRENT_BINARY_DIR}/${testDir}/${matrixListingName}
	          --rootdir=${CMAKE_CURRENT_BINARY_DIR}/${testDir} )

  	FILE( READ ${CMAKE_CURRENT_BINARY_DIR}/${testDir}/${matrixListingName} matrixFiles )
  	STRING( REGEX REPLACE ";" "\\\\;" matrixFiles "${matrixFiles}" )
  	STRING( REGEX REPLACE "\n" ";" matrixFiles "${matrixFiles}" )
 
  	set( INC 0 )
	FOREACH( matrixFile ${matrixFiles} )
		MATH(EXPR INC "${INC}+1")
		
		set( computed ${CMAKE_CURRENT_BINARY_DIR}/${testDir}/${matrixFile}.spectrum )
		set( reference ${CMAKE_CURRENT_BINARY_DIR}/${testDir}/${matrixFile}.${referenceType} )
		
    	add_test( NAME ${testName}_vs_${referenceType}_np${numProc}_matrix${INC}
    	          COMMAND ${CMAKE_COMMAND} -E compare_files ${computed} ${reference} )
    	          
		set_property( TEST ${testName}_vs_${referenceType}_np${numProc}_matrix${INC} APPEND PROPERTY DEPENDS ${testName}_compute_np${numProc} )
	ENDFOREACH( matrixFile ${matrixFiles} )

endfunction()