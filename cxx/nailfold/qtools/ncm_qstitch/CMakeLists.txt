# nailfold/qtools/ncm_qstitch/CMakeLists.txt
#
#   Phil Tresadern
#

FIND_PACKAGE( SuiteSparse REQUIRED )

# Add folder where .h files are created (from the .ui files) to the search path
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

# Convert these .ui files to .h
QT4_WRAP_UI( ncm_qstitch_ui_h_files
	ncm_qstitch_main.ui
)	

# Convert these .h files to .cxx
QT4_WRAP_CPP( ncm_qstitch_moc_files 
	ncm_qstitch_main.h
	ncm_qstitch_graphicsview.h
)

SET(ncm_qstitch_sources
	ncm_qstitch.cpp
	ncm_qstitch_main.h ncm_qstitch_main.cpp 
	ncm_qstitch_graphicsview.h ncm_qstitch_graphicsview.cxx
)

LINK_DIRECTORIES( ${LAPACK_LIBRARY_DIR} ${LAPACK_BLAS_LIBRARY_DIR} ${LAPACK_F2C_LIBRARY_DIR} ${SPARSELM_LIBRARY_DIR} )
ADD_EXECUTABLE( ncm_qstitch 
	${ncm_qstitch_sources}
	${ncm_qstitch_ui_h_files}
	${ncm_qstitch_moc_files} 
)

TARGET_LINK_LIBRARIES( ncm_qstitch 
	vbl vil vil_algo vnl 
  lapack libf2c blas suitesparse metis
	ncm ncm_flow ncm_qt qcore qvcr ${QT_LIBRARIES} ${QT_QTSQL_LIBRARY})
