# nailfold/CMakeLists.txt
#
#   Phil Tresadern
#

doxygen_add_book(isbe_apps/nailfold 
  "Nailfold Capillaroscopy Documentation")
doxygen_add_package(isbe_apps/nailfold "Nailfold Capillaroscopy")

#INCLUDE_DIRECTORIES( ${ISBE_SEARCH_INCLUDE_DIR} )

OPTION(BUILD_ISBE_APPS_NAILFOLD_QTOOLS "Should the tools in isbe_apps/nailfold be built?" NO)

# Add DirectShow includes if needed
INCLUDE( ${MODULE_PATH}/NewCMake/FindDirectShow.cmake )
IF(DIRECTSHOW_FOUND)
  INCLUDE_DIRECTORIES( ${DIRECTSHOW_INCLUDE_DIRS} )
ENDIF(DIRECTSHOW_FOUND)

# add option to link in Lumenera camera libraries
#OPTION(LUMENERA_INCLUDE_DIR "Include directory for Lumenera Camera SDK" "")
#OPTION(LUMENERA_LIB "Library directory for Lumenera Camera SDK" "")
# include Lumenera files if defined
#IF ( (DEFINED ${LUMENERA_INCLUDE_DIR}) AND
#     (DEFINED ${LUMENERA_LIB_DIR}) )
#INCLUDE_DIRECTORIES( ${LUMENERA_INCLUDE_DIR} )

SET (ncm_sources
  vil_gaussian_derivatives.h vil_gaussian_derivatives.txx
  vil_suppress_non_max_dir.h vil_suppress_non_max_dir.cxx
  vil_connected_components.h vil_connected_components.cxx
  vil_find_connected_components.h vil_find_connected_components.cxx
  vil_dtcwt.h vil_dtcwt.cxx
	vil_multiscale_gaussian_2nd_derivative.h vil_multiscale_gaussian_2nd_derivative.cxx

  ncm_sharpness_evaluator.h ncm_sharpness_evaluator.cxx

  ncm_frame_aligner.h ncm_frame_aligner.cxx
  ncm_mosaic_maker.h ncm_mosaic_maker.cxx

  ncm_encrypt.h ncm_encrypt.cxx
  ncm_image_server.h ncm_image_server.cxx
  ncm_registry_accessor.h ncm_registry_accessor.cxx

  ncm_server_handler.h ncm_server_handler.cxx
  ncm_scp_handler.h ncm_scp_handler.cxx
  ncm_sftp_handler.h ncm_sftp_handler.cxx

  ncm_vessel.h ncm_vessel.cxx ncm_vessel_io.cxx
  ncm_vessel_properties.h ncm_vessel_properties.cxx ncm_vessel_properties_io.cxx
  ncm_apex.h ncm_apex.cxx ncm_apex_io.cxx
  ncm_haemorrhage.h ncm_haemorrhage.cxx ncm_haemorrhage_io.cxx
  ncm_annotation.h ncm_annotation.cxx ncm_annotation_io.cxx
  ncm_image_grade.h ncm_image_grade.cxx

  ncm_apt_server.h ncm_apt_server.cxx
  ncm_apt_controller_base.h ncm_apt_controller_base.cxx
  ncm_apt_controller.h ncm_apt_controller.cxx
  ncm_apt_null_controller.h ncm_apt_null_controller.cxx

	ncm_hand.h

#  ncm_camera.h ncm_camera.cxx
#  ncm_lumenera_camera.h ncm_lumenera_camera.cxx

  ncm_timebar.h
)

AUX_SOURCE_DIRECTORY(Templates ncm_sources)

# Find headers and libraries for Thorlabs motorized platforms
#SET( THORLABS_INCLUDE_DIR "C:/Program Files/Thorlabs/APT/APT API/APTAPI" )
SET( THORLABS_INCLUDE_DIR ${THORLABS_INCLUDE_DIR} CACHE PATH "Where is the include dir for the Thorlabs APT API" )
INCLUDE_DIRECTORIES( ${THORLABS_INCLUDE_DIR} )
FIND_LIBRARY(THORLABS_APTAPI_LIB ${THORLABS_INCLUDE_DIR})
MESSAGE( STATUS "Thorlabs include dir: ${THORLABS_INCLUDE_DIR}" )
MESSAGE( STATUS "Thorlabs APT API lib: ${THORLABS_APTAPI_LIB}" )

ADD_LIBRARY ( ncm ${ncm_sources} )
TARGET_LINK_LIBRARIES( ncm vil vil_algo ${THORLABS_APTAPI_LIB})

# Disable the annoying warnings caused by Thorlabs's API that uses a deprecated 
# command (std::_Copy_opt).
SET_TARGET_PROPERTIES( ncm PROPERTIES COMPILE_FLAGS /wd4996 )

# Build optical flow library
ADD_SUBDIRECTORY( flow )

# include tools in build
ADD_SUBDIRECTORY( tools )

IF( BUILD_ISBE_APPS_NAILFOLD_QTOOLS )
  ADD_SUBDIRECTORY( qtools )
ENDIF( BUILD_ISBE_APPS_NAILFOLD_QTOOLS )

#IF( BUILD_TESTING )
  ADD_SUBDIRECTORY(tests)
#ENDIF( BUILD_TESTING )
