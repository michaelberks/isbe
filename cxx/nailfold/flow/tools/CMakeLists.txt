# nailfold/CMakeLists.txt
#
#   Phil Tresadern
#

ADD_EXECUTABLE ( ncm_estimate_flow ncm_estimate_flow.cxx )
TARGET_LINK_LIBRARIES( ncm_estimate_flow vil_algo ncm_flow )

#ADD_SUBDIRECTORY( test_nag )
#ADD_SUBDIRECTORY( test_sparselm )
ADD_SUBDIRECTORY( test_cholmod )
