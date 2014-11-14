# nailfold/qtools/ncm_qseries/CMakeLists.txt
#
#   Mike Berks
#

# Add folder where .h files are created (from the .ui files) to the search path
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

# Convert src/<myname>.ui files to bin/ui_<myname>.h equivalents
QT4_WRAP_UI( ncm_qseries_ui_h_files
	ncm_qseries_gui.ui
	ncm_qseries_options.ui
)	

# Convert src/<myname>.h files to bin/moc_<myname>.cxx equivalents
QT4_WRAP_CPP( ncm_qseries_moc_files 
	ncm_qseries_gui.h
	ncm_qseries_options.h
)

# Convert these resource files
QT4_ADD_RESOURCES( ncm_qseries_res_files
	ncm_qseries.qrc	
)

SET(ncm_qseries_sources
	ncm_qseries_preferences.h ncm_qseries_preferences.cxx

	ncm_qseries.cpp
	ncm_qseries_gui.h ncm_qseries_gui.cpp
	ncm_qseries_options.h ncm_qseries_options.cpp
	ncm_qseries_image_manager.h ncm_qseries_image_manager.cxx
	ncm_qseries_grade.h ncm_qseries_grade.cxx
)
ADD_EXECUTABLE( ncm_qseries 
	${ncm_qseries_sources}
	${ncm_qseries_ui_h_files}
	${ncm_qseries_moc_files}
	${ncm_qseries_res_files} 
)
TARGET_LINK_LIBRARIES( ncm_qseries vil vidl vimt ncm qcore qvcr ncm_qt ${QT_LIBRARIES})

IF(DIRECTSHOW_FOUND)
#  TARGET_LINK_LIBRARIES( ncm_qseries ${DIRECTSHOW_LIBRARIES} )
ENDIF(DIRECTSHOW_FOUND)

