# nailfold/qtools/ncm_qfinal/CMakeLists.txt
#
#   Mike Berks
#

# Add folder where .h files are created (from the .ui files) to the search path
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

# Convert src/<myname>.ui files to bin/ui_<myname>.h equivalents
QT4_WRAP_UI( ncm_qfinal_ui_h_files
	ncm_qfinal_gui.ui
	ncm_qfinal_options.ui
	ncm_qfinal_search.ui
)	

# Convert src/<myname>.h files to bin/moc_<myname>.cxx equivalents
QT4_WRAP_CPP( ncm_qfinal_moc_files 
	ncm_qfinal_gui.h
	ncm_qfinal_options.h
	ncm_qfinal_search.h
)

# Convert these resource files
QT4_ADD_RESOURCES( ncm_qfinal_res_files
	ncm_qfinal.qrc	
)

SET(ncm_qfinal_sources
	ncm_qfinal_preferences.h ncm_qfinal_preferences.cxx

	ncm_qfinal.cpp
	ncm_qfinal_gui.h ncm_qfinal_gui.cpp
	ncm_qfinal_options.h ncm_qfinal_options.cpp
	ncm_qfinal_search.h ncm_qfinal_search.cpp
)
ADD_EXECUTABLE( ncm_qfinal 
	${ncm_qfinal_sources}
	${ncm_qfinal_ui_h_files}
	${ncm_qfinal_moc_files}
	${ncm_qfinal_res_files} 
)
TARGET_LINK_LIBRARIES( ncm_qfinal vil vidl vimt ncm qcore qvcr ncm_qt ${QT_LIBRARIES} ${QT_QTSQL_LIBRARY})
MESSAGE( STATUS "NCM final QT libraries: ${QT_LIBRARIES}" )
MESSAGE( STATUS "NCM final QT SQL library: ${QT_QTSQL_LIBRARY}" )
IF(DIRECTSHOW_FOUND)
#  TARGET_LINK_LIBRARIES( ncm_qfinal ${DIRECTSHOW_LIBRARIES} )
ENDIF(DIRECTSHOW_FOUND)

