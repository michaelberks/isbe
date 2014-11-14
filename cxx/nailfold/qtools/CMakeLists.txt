# nailfold/qtools/CMakeLists.txt
#
#   Phil Tresadern
#

# include Tim's Qt components
INCLUDE_DIRECTORIES( ${ISBE_CORE_QT4_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${ISBE_MARK_QT4_INCLUDE_DIR} )

# Add support for Qt4
INCLUDE(${CMAKE_ROOT}/Modules/FindQt4.cmake)
INCLUDE(${QT_USE_FILE})

# Add folder where .h files are created (from the .ui files) to the search path
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

# add option to link in DMK camera libraries
OPTION(USE_DMK_CAMERAS "Should we build the DMK cameras module?" YES)
IF ( USE_DMK_CAMERAS )
  SET( DMK_INCLUDE_DIR ${DMK_INCLUDE_DIR} CACHE PATH "Where is the include dir for the DMK camera SDK" )
  SET( DMK_LIB_DIR ${DMK_LIB_DIR} CACHE PATH "Library directory for DMK Camera SDK (debug)" )
  SET( DMK_LIB_DIR_R ${DMK_LIB_DIR_R} CACHE PATH "Library directory for DMK Camera SDK (release)" )
  MESSAGE( STATUS "DMK include dir: ${DMK_INCLUDE_DIR}" )
  INCLUDE_DIRECTORIES( ${DMK_INCLUDE_DIR} )
ENDIF ( USE_DMK_CAMERAS )

# add option to link in Lumenera camera libraries
#OPTION(LUMENERA_INCLUDE_DIR "Include directory for Lumenera Camera SDK" "")
#OPTION(LUMENERA_LIB "Library directory for Lumenera Camera SDK" "")
# include Lumenera files if defined
#IF ( (DEFINED ${LUMENERA_INCLUDE_DIR}) AND
#     (DEFINED ${LUMENERA_LIB_DIR}) )
#INCLUDE_DIRECTORIES( ${LUMENERA_INCLUDE_DIR} )

# add files to library of Qt widgets for NCM applications

# Convert src/<myname>.ui files to bin/ui_<myname>.h equivalents
QT4_WRAP_UI( ncm_qt_ui_h_files
	ncm_qsite_selector.ui
)	

# Convert these .h files to .cxx (these classes all have the capacity to use signals and slots, unlike pure graphicsitems)
QT4_WRAP_CPP( ncm_qt_moc_files
	nailfold_graphics_view.h

	ncm_qimagehandler.h
	ncm_qframe_aligner.h
	ncm_qframe_aligner_file.h
	ncm_qmosaic_maker.h

	ncm_qscene.h
	ncm_qsceneview.h
	ncm_qthumbnail.h
	ncm_qvessel_properties.h
	ncm_exposure_spinbox.h
	
#	ncm_camera.h
	ncm_video_frame.h
	ncm_video_frame_queue.h
#	ncm_dmk_camera.h
#	ncm_dmk_grab_frame_filter.h

	ncm_qdata_manager.h
	ncm_qdata_saver.h
	ncm_qdata_processor.h
	ncm_qdata_stream.h

	ncm_qapt_server.h

	ncm_qapexitem.h

	ncm_qsite_selector.h
)

SET( ncm_qt_sources
	ncm_qglobal.h

	qshow_image.h

	ncm_qvesselitem_appearance.h ncm_qvesselitem_appearance.cxx
	ncm_qapexitem_appearance.h ncm_qapexitem_appearance.cxx
	ncm_qhaemorrhageitem_appearance.h ncm_qhaemorrhageitem_appearance.cxx

# Classes based on QGraphicsItem
	ncm_qvesselitem.h ncm_qvesselitem.cxx 
	ncm_qapexitem.h ncm_qapexitem.cxx
	ncm_qhaemorrhageitem.h ncm_qhaemorrhageitem.cxx
	ncm_qscenegrid.h ncm_qscenegrid.cxx 

# QObjects
	nailfold_graphics_view.h nailfold_graphics_view.cxx

	ncm_qimagehandler.h ncm_qimagehandler.cxx
	ncm_qframe_aligner.h ncm_qframe_aligner.cxx
	ncm_qframe_aligner_file.h ncm_qframe_aligner_file.cxx
	ncm_qmosaic_maker.h ncm_qmosaic_maker.cxx

	ncm_qscene.h ncm_qscene.cxx
	ncm_qsceneview.h ncm_qsceneview.cxx
	ncm_qthumbnail.h ncm_qthumbnail.cxx
	ncm_qvessel_properties.h ncm_qvessel_properties.cxx
	ncm_exposure_spinbox.h ncm_exposure_spinbox.cpp

	ncm_qapt_server.h ncm_qapt_server.cxx
	
# Classes to implement cameras
	ncm_video_frame.h ncm_video_frame.cxx
	ncm_video_frame_queue.h ncm_video_frame_queue.cxx

#	ncm_camera.h ncm_camera.cxx
#  	ncm_lumenera_camera.h ncm_lumenera_camera.cxx
#  	ncm_dmk_camera.h ncm_dmk_camera.cxx
#	ncm_dmk_grab_frame_filter.h ncm_dmk_grab_frame_filter.cxx	

	ncm_qdata_manager.h ncm_qdata_manager.cxx
	ncm_qdata_saver.h ncm_qdata_saver.cxx
	ncm_qdata_processor.h ncm_qdata_processor.cxx
	ncm_qdata_stream.h ncm_qdata_stream.cxx

#Self-contained GUI to select host site
	ncm_qsite_selector.h ncm_qsite_selector.cpp

	convert_ncm_video_frame_to_qimage.h

	#NCM classes that use QMySQL libraries
	ncm_user.h ncm_user.cxx ncm_user_io.cxx
	ncm_subject.h ncm_subject.cxx ncm_subject_io.cxx
	ncm_image_session.h ncm_image_session.cxx ncm_image_session_io.cxx
	ncm_image_sequence.h ncm_image_sequence.cxx ncm_image_sequence_io.cxx
)

ADD_LIBRARY( ncm_qt
	${ncm_qt_sources}
	${ncm_qt_moc_files}
	${ncm_qt_ui_h_files}
)

# Qt interface to Thorlabs motorized platforms
#SUBDIRS( ncm_qapt_util )

# Tools that use the camera classes
SUBDIRS( ncm_cameras )
SUBDIRS( ncm_qcapture ) # standalone image capture tool (for testing camera interfaces)
SUBDIRS( ncm_qcapture2 ) # standalone image capture tool for multiple camera
#SUBDIRS( ncm_qcapture_mdi ) # standalone image capture tool (for testing camera interfaces)

#SUBDIRS( ncm_qsandpit )

#SUBDIRS( dmk_qcapture ) # standalone image capture tool for DMK cameras

# Tool that stitches images into a mosaic and displays the result to the user
SUBDIRS( ncm_qstitch )

# Standalone markup tool used in RSA study
SUBDIRS( ncm_qmarkup ) 

# Standalone markup tool for grading temporal series of images
SUBDIRS( ncm_qseries )

# Integrated final system
SUBDIRS( ncm_qfinal ) 

# Integrated test system Phil was developing
SUBDIRS( ncm_qtool ) 


