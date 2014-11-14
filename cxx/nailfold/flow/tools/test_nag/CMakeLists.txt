PROJECT( test_nag )

cmake_minimum_required(VERSION 2.6)

#INCLUDE_DIRECTORIES( ${NAG_INCLUDE_DIR} )
#LINK_DIRECTORIES( ${NAG_LIBRARY_DIR} )

ADD_EXECUTABLE( test_nag test_nag.cxx )
TARGET_LINK_LIBRARIES( test_nag nagc_nag_MT )


