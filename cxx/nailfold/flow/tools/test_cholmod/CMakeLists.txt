PROJECT( test_cholmod )

cmake_minimum_required(VERSION 2.6)

#INCLUDE_DIRECTORIES( ${SPARSELM_INCLUDE_DIR} )
#LINK_DIRECTORIES( ${LAPACK_LIBRARY_DIR} ${SuiteSparse_LIBS_DIR} ${SPARSELM_LIBRARY_DIR} )
#ADD_EXECUTABLE( test_sparselm test_sparselm.cxx )
#TARGET_LINK_LIBRARIES( test_sparselm lapack blas libf2c suitesparse metis splm )

LINK_DIRECTORIES( ${LAPACK_LIBRARY_DIR} ${LAPACK_BLAS_LIBRARY_DIR} ${LAPACK_F2C_LIBRARY_DIR} ${SPARSELM_LIBRARY_DIR} )
ADD_EXECUTABLE( test_cholmod test_cholmod.cxx )
TARGET_LINK_LIBRARIES( test_cholmod lapack blas libf2c suitesparse metis )


