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
  SET( DMK_LIB_DIR ${DMK_LIB_DIR} CACHE PATH "Root directory for DMK Camera SDK libraries" )
#  SET( DMK_LIB_DIR_R ${DMK_LIB_DIR_R} CACHE PATH "Library directory for DMK Camera SDK (release)" )
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

# Convert these .h files to .cxx (these classes all have the capacity to use signals and slots)
QT4_WRAP_CPP( ncm_cameras_moc_files
	ncm_camera.h
	ncm_dmk_camera.h
	ncm_dmk_grab_frame_filter.h
)

SET( ncm_cameras_sources
# Classes to implement cameras
	ncm_camera.h ncm_camera.cxx
  	ncm_lumenera_camera.h ncm_lumenera_camera.cxx
  	ncm_dmk_camera.h ncm_dmk_camera.cxx
	ncm_dmk_grab_frame_filter.h ncm_dmk_grab_frame_filter.cxx
	ncm_qcamera_gui.h ncm_qcamera_gui.cxx
)

ADD_LIBRARY( ncm_cameras
	${ncm_cameras_sources}
	${ncm_cameras_moc_files}
)
