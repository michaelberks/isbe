PROJECT( ncm_sandpit )

ADD_EXECUTABLE( ncm_sandpit ncm_sandpit.cxx )
#TARGET_LINK_LIBRARIES( ncm_sandpit vil )
TARGET_LINK_LIBRARIES( ncm_sandpit vil vil_algo vnl vsl vul vbl vgl ncm )

ADD_EXECUTABLE( ncm_encrypt_string ncm_encrypt_string.cxx )
TARGET_LINK_LIBRARIES( ncm_encrypt_string ncm )

ADD_EXECUTABLE( ncm_thorlabs ncm_thorlabs.cxx )
TARGET_LINK_LIBRARIES( ncm_thorlabs vsl vgl ncm ${THORLABS_APTAPI_LIB} )

ADD_SUBDIRECTORY( ncm_stitch )
ADD_SUBDIRECTORY( ncm_sandpit_mb )

