# nailfold/qtools/ncm_qmarkup/CMakeLists.txt
#
#   Phil Tresadern
#

# Add folder where .h files are created (from the .ui files) to the search path
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

# Convert src/<myname>.ui files to bin/ui_<myname>.h equivalents
QT4_WRAP_UI( ncm_qmarkup_ui_h_files
	ncm_qmarkup_gui.ui
	ncm_qmarkup_options.ui
)	

# Convert src/<myname>.h files to bin/moc_<myname>.cxx equivalents
QT4_WRAP_CPP( ncm_qmarkup_moc_files 
	ncm_qmarkup_gui.h
	ncm_qmarkup_options.h
)

# Convert these resource files
QT4_ADD_RESOURCES( ncm_qmarkup_res_files
	ncm_qmarkup.qrc	
)

SET(ncm_qmarkup_sources
	ncm_qmarkup_preferences.h ncm_qmarkup_preferences.cxx

	ncm_qmarkup.cpp
	ncm_qmarkup_gui.h ncm_qmarkup_gui.cpp
	ncm_qmarkup_options.h ncm_qmarkup_options.cpp
	ncm_qmarkup_server.h ncm_qmarkup_server.cxx
)
ADD_EXECUTABLE( ncm_qmarkup 
	${ncm_qmarkup_sources}
	${ncm_qmarkup_ui_h_files}
	${ncm_qmarkup_moc_files}
	${ncm_qmarkup_res_files} 
)
TARGET_LINK_LIBRARIES( ncm_qmarkup vil vidl vimt ncm qcore qvcr ncm_qt ${QT_LIBRARIES})

IF(DIRECTSHOW_FOUND)
#  TARGET_LINK_LIBRARIES( ncm_qmarkup ${DIRECTSHOW_LIBRARIES} )
ENDIF(DIRECTSHOW_FOUND)

