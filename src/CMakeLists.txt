SET( SOURCES
     paladin.cpp )

FIND_PACKAGE( MPI REQUIRED )

INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )

IF( APPLE )
	MESSAGE( "Searching for framework Accelerate..." )
	FIND_LIBRARY( BLASLAPACK Accelerate )
	IF( NOT BLASLAPACK )
		MESSAGE( FATAL_ERROR "Mac OS X: framework Accelerate not found!" )
    ELSE()
        MESSAGE( "Found framework Accelerate!" )
	ENDIF()
ELSE()
    MESSAGE( "Searching for BLAS and LAPACK..." )
    FIND_PACKAGE( BLAS REQUIRED )
    FIND_PACKAGE( LAPACK REQUIRED )
    SET( BLASLAPACK ${LAPACK_LIBRARIES} )
    IF( NOT BLASLAPACK )
        MESSAGE( FATAL_ERROR "BLAS and LAPACK not found!" )
    ELSE()
        MESSAGE( "Found BLAS and LAPACK!" )
    ENDIF()
ENDIF()


ADD_EXECUTABLE( paladin paladin.cpp )
TARGET_LINK_LIBRARIES( paladin ${BLASLAPACK} ${MPI_LIBRARIES} )