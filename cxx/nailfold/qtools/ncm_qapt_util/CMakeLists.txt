# nailfold/qtools/CMakeLists.txt
#
#   Phil Tresadern
#

# Add folder where .h files are created (from the .ui files) to the search path
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

# Convert src/<myname>.ui files to bin/ui_<myname>.h equivalents
QT4_WRAP_UI( ncm_qapt_util_ui_h_files
	ncm_qapt_util_gui.ui
)	

# Convert src/<myname>.h files to bin/moc_<myname>.cxx equivalents
QT4_WRAP_CPP( ncm_qapt_util_moc_files 
	ncm_qapt_gui_graphicsview.h
	ncm_qapt_util_gui.h
)

# Convert these resource files
QT4_ADD_RESOURCES( ncm_qapt_util_res_files
#	ncm_qapt_util_gui.qrc
)

SET(ncm_qapt_util_sources
	ncm_qapt_gui_graphicsview.h ncm_qapt_gui_graphicsview.cxx
	ncm_qapt_util_gui.h ncm_qapt_util_gui.cpp
	ncm_qapt_util.cpp
)

ADD_EXECUTABLE( ncm_qapt_util
	${ncm_qapt_util_sources}
	${ncm_qapt_util_ui_h_files}
	${ncm_qapt_util_moc_files}
	${ncm_qapt_util_res_files} 
)
TARGET_LINK_LIBRARIES( ncm_qapt_util vil vidl vimt ncm qcore qvcr ncm_qt ${QT_LIBRARIES} )

IF(DIRECTSHOW_FOUND)
#  TARGET_LINK_LIBRARIES( ncm_qapt_util_gui ${DIRECTSHOW_LIBRARIES} )
ENDIF(DIRECTSHOW_FOUND)

