ADD_EXECUTABLE( ncm_test_gaussian_2nd_derivative ncm_test_gaussian_2nd_derivative.cxx )
TARGET_LINK_LIBRARIES( ncm_test_gaussian_2nd_derivative vil vil_algo vnl vsl ncm )

ADD_EXECUTABLE( ncm_test_find_connected_components ncm_test_find_connected_components.cxx )
TARGET_LINK_LIBRARIES( ncm_test_find_connected_components vil vil_algo vnl vsl ncm )

ADD_EXECUTABLE( ncm_test_dtcwt ncm_test_dtcwt.cxx )
TARGET_LINK_LIBRARIES( ncm_test_dtcwt vil vil_algo vnl vsl ncm )

ADD_EXECUTABLE( ncm_test_registry_accessor ncm_test_registry_accessor.cxx )
TARGET_LINK_LIBRARIES( ncm_test_registry_accessor ncm )

ADD_EXECUTABLE( ncm_test ncm_test.cxx )
TARGET_LINK_LIBRARIES( ncm_test vsl vgl mbl ncm )

