# nailfold/qtools/ncm_qtool/CMakeLists.txt
#
#   Phil Tresadern
#

# Add folder where .h files are created (from the .ui files) to the search path
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

# Convert these .ui files to .h
QT4_WRAP_UI( ncm_qtool_ui_h_files
	ncm_qtool_tmp.ui
)	

# Convert these .h files to .cxx
QT4_WRAP_CPP( ncm_qtool_moc_files 
	ncm_qtool_tmp.h
)

SET(ncm_qtool_sources
	ncm_qtool.cpp
	ncm_qtool_tmp.h ncm_qtool_tmp.cpp
)
ADD_EXECUTABLE( ncm_qtool 
	${ncm_qtool_sources}
	${ncm_qtool_ui_h_files}
	${ncm_qtool_moc_files} 
)
TARGET_LINK_LIBRARIES( ncm_qtool vil vidl vimt ncm qcore qvcr ncm_qt ${QT_LIBRARIES})

IF(DIRECTSHOW_FOUND)
#  TARGET_LINK_LIBRARIES( ncm_qtool ${DIRECTSHOW_LIBRARIES} )
ENDIF(DIRECTSHOW_FOUND)

