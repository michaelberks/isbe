# nailfold/qtools/ncm_qcapture/CMakeLists.txt
#
#   Phil Tresadern
#

# Add folder where .h files are created (from the .ui files) to the search path
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

# Convert these .ui files to .h
QT4_WRAP_UI( ncm_qsandpit_ui_h_files
)	

# Convert these .h files to .cxx
QT4_WRAP_CPP( ncm_qsandpit_moc_files 
)

SET(ncm_qsandpit_sources
	ncm_qsandpit.cxx
)

ADD_EXECUTABLE( ncm_qsandpit
	${ncm_qsandpit_sources}
	${ncm_qsandpit_ui_h_files}
	${ncm_qsandpit_moc_files} 
)

TARGET_LINK_LIBRARIES( ncm_qsandpit 
	vil vil_algo vidl vimt 
	ncm ncm_qt ncm_cameras 
	qcore qvcr ${QT_LIBRARIES}
)

