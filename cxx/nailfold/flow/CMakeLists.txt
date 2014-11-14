# nailfold/CMakeLists.txt
#
#   Phil Tresadern
#

FIND_PACKAGE( SuiteSparse REQUIRED )

#SET( SuiteSparse_INCLUDE_DIR "U:/uman-src/suitesparse-metis-for-windows-1.0.0/SuiteSparse"
#     CACHE PATH "Location of the SuiteSparse header files" )
SET( SPARSELM_INCLUDE_DIR "U:/uman-src/sparselm-1.3"
     CACHE PATH "Location of the sparselm header files" )
SET( NAG_INCLUDE_DIR "C:/Program Files/NAG/CL23/clw3223dal/include"
     CACHE PATH "Location of the NAG header files" )

SET( SPARSELM_LIBRARY_DIR C:/uman-bin/sparselm
     CACHE PATH "Location of the sparselm libraries" )
SET( LAPACK_LIBRARY_DIR {$LAPACK_LIBRARY_DIR }
     CACHE PATH "Location of the lapack libraries" )
SET( LAPACK_BLAS_LIBRARY_DIR {$LAPACK_BLAS_LIBRARY_DIR}
     CACHE PATH "Location of the lapack libraries" )
SET( LAPACK_F2C_LIBRARY_DIR {$LAPACK_F2C_LIBRARY_DIR}
     CACHE PATH "Location of the lapack libraries" )
SET( NAG_LIBRARY_DIR "C:/Program Files/NAG/CL23/clw3223dal/lib"
     CACHE PATH "Location of the NAG libraries" )

LINK_DIRECTORIES( ${LAPACK_LIBRARY_DIR} ${LAPACK_BLAS_LIBRARY_DIR} ${LAPACK_F2C_LIBRARY_DIR} ${SuiteSparse_LIBS_DIR} )

INCLUDE_DIRECTORIES( ${SPARSELM_INCLUDE_DIR} )
LINK_DIRECTORIES( ${SPARSELM_LIBRARY_DIR} )

#INCLUDE_DIRECTORIES( ${NAG_INCLUDE_DIR} )
#LINK_DIRECTORIES( ${NAG_LIBRARY_DIR} )


FOREACH( class 
	cost_function
#		cost_function_sparselm
#		cost_function_vxl
#		cost_function_nag
		cost_function_cholmod
		brightness_cost
		smoothness_cost
	robustifier
		lorentzian_robustifier
		charbonnier_robustifier
		geman-mcclure_robustifier
		null_robustifier
	minimizer
	differentiator
	interpolator
	field
	estimator
#		estimator_sparselm
#		estimator_vxl
#		estimator_nag
		estimator_cholmod
)
  SET ( ncm_flow_sources 
    ${ncm_flow_sources} ncm_flow_${class}.h ncm_flow_${class}.cxx
  )
ENDFOREACH()

ADD_LIBRARY ( ncm_flow ${ncm_flow_sources} )
TARGET_LINK_LIBRARIES( ncm_flow 
	vbl vil vil_algo vnl vsl vul ncm 
	lapack blas libf2c suitesparse metis 
#	splm
#	nagc_nag_MT
)

ADD_SUBDIRECTORY( tools )
