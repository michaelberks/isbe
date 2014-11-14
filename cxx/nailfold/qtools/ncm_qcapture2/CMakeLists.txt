# nailfold/qtools/ncm_qcapture/CMakeLists.txt
#
#   Phil Tresadern
#

# Add folder where .h files are created (from the .ui files) to the search path
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

IF ( USE_DMK_CAMERAS )
  LINK_DIRECTORIES( ${DMK_LIB_DIR} )

#  Note: This directory is not needed so long as DMK_LIB_DIR points to the root directory
#        i.e. without the 'Debug' or 'Release' suffix. This is added automatically at build
#        time.
#  LINK_DIRECTORIES( ${DMK_LIB_DIR_R} )

  MESSAGE( STATUS "DMK link lib: ${DMK_LIB_DIR}" )
ENDIF ( USE_DMK_CAMERAS )

# Convert these .ui files to .h
QT4_WRAP_UI( ncm_qcapture2_ui_h_files
	ncm_qcapture2_gui.ui
)	

# Convert these .h files to .cxx
QT4_WRAP_CPP( ncm_qcapture2_moc_files 
	ncm_qcapture2_gui.h
	ncm_qcapture2_saver.h
	ncm_qcapture2_processor.h
	ncm_qcapture2_graphicsview.h
)

# Convert these resource files
QT4_ADD_RESOURCES( ncm_qcapture2_res_files
	ncm_qcapture2.qrc	
)

SET(ncm_qcapture2_sources
	ncm_qcapture2.cpp
	ncm_qcapture2_gui.h ncm_qcapture2_gui.cpp 
	
	ncm_qcapture2_data_manager.h ncm_qcapture2_data_manager.cxx
	ncm_qcapture2_saver.h ncm_qcapture2_saver.cpp
	ncm_qcapture2_processor.h ncm_qcapture2_processor.cxx
	ncm_qcapture2_graphicsview.h ncm_qcapture2_graphicsview.cxx
)

ADD_EXECUTABLE( ncm_qcapture2
	${ncm_qcapture2_sources}
	${ncm_qcapture2_ui_h_files}
	${ncm_qcapture2_moc_files} 
	${ncm_qcapture2_res_files}
)

TARGET_LINK_LIBRARIES( ncm_qcapture2 vil vil_algo vidl vimt ncm qcore qvcr ncm_qt ncm_cameras ${QT_LIBRARIES})

IF(DIRECTSHOW_FOUND)
#  TARGET_LINK_LIBRARIES( ncm_qcapture2 ${DIRECTSHOW_LIBRARIES} )
ENDIF(DIRECTSHOW_FOUND)

