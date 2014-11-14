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
QT4_WRAP_UI( ncm_qcapture_ui_h_files
	ncm_qcapture_gui.ui
	ncm_qcapture_search.ui
	ncm_qcapture_new_subject.ui
	ncm_qcapture_end_sequence.ui
	ncm_qcapture_end_session.ui
	ncm_qcapture_options.ui
)	

# Convert these .h files to .cxx
QT4_WRAP_CPP( ncm_qcapture_moc_files 
	ncm_qcapture_gui.h
	ncm_qcapture_saver.h
	ncm_qcapture_processor.h
	ncm_qcapture_aligner.h
	ncm_qcapture_mosaic_maker.h
	ncm_qcapture_graphicsview.h
	ncm_qcapture_scene.h
	ncm_qcapture_search.h
	ncm_qcapture_new_subject.h
	ncm_qcapture_end_sequence.h
	ncm_qcapture_end_session.h
	ncm_qcapture_options.h
)

# Convert these resource files
QT4_ADD_RESOURCES( ncm_qcapture_res_files
	ncm_qcapture.qrc	
)
	
SET(ncm_qcapture_sources
	ncm_qcapture.cpp
)

FOREACH( class
	data_manager data_stream
	processor aligner mosaic_maker
	scene graphicsview
)
	SET( ncm_qcapture_sources
		${ncm_qcapture_sources} ncm_qcapture_${class}.h ncm_qcapture_${class}.cxx
	)
ENDFOREACH( class )

# Would be easier to rename files but I won't for now
FOREACH( class
	gui saver search new_subject end_sequence end_session options preferences
)
	SET( ncm_qcapture_sources
		${ncm_qcapture_sources} ncm_qcapture_${class}.h ncm_qcapture_${class}.cpp
	)
ENDFOREACH( class )

ADD_EXECUTABLE( ncm_qcapture 
	${ncm_qcapture_sources}
	${ncm_qcapture_ui_h_files}
	${ncm_qcapture_moc_files} 
	${ncm_qcapture_res_files}
)
TARGET_LINK_LIBRARIES( ncm_qcapture vil vil_algo vidl vimt ncm qcore qvcr ncm_qt ncm_cameras ${QT_LIBRARIES} ${QT_QTSQL_LIBRARY})

# Disable the annoying warnings caused by Thorlabs's API that uses a deprecated 
# command (std::_Copy_opt).
SET_TARGET_PROPERTIES( ncm_qcapture PROPERTIES COMPILE_FLAGS /wd4996 )

IF(DIRECTSHOW_FOUND)
#  TARGET_LINK_LIBRARIES( ncm_qcapture ${DIRECTSHOW_LIBRARIES} )
ENDIF(DIRECTSHOW_FOUND)

