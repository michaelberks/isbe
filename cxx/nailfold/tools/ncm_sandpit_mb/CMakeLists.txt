# nailfold/tools/ncm_sandpit_mb/CMakeLists.txt
#
#   Mike Berks
#

# Add folder where .h files are created (from the .ui files) to the search path
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

SET(ncm_sandpit_mb_sources
	ncm_sandpit_mb.cxx
	ncm_test_apex_aam.cxx ncm_test_apex_aam.h
	ncm_base_iih_patch_trainer.cxx ncm_base_iih_patch_trainer.h
	ncm_base_iih_patch_extractor.cxx ncm_base_iih_patch_extractor.h
	ncm_base_iih_patch_predictor.cxx ncm_base_iih_patch_predictor.h
	ncm_apex_regression_sampler.cxx ncm_apex_regression_sampler.h
	ncm_apex_classification_sampler.cxx ncm_apex_classification_sampler.h
	ncm_vessel_patch_trainer.cxx ncm_vessel_patch_trainer.h
	ncm_vessel_detection_sampler.cxx ncm_vessel_detection_sampler.h
	ncm_gaussian_patch_extractor.cxx ncm_gaussian_patch_extractor.h
	
)
ADD_EXECUTABLE( ncm_sandpit_mb 
	${ncm_sandpit_mb_sources} 
)
TARGET_LINK_LIBRARIES( ncm_sandpit_mb vil vil_algo vnl vsl ncm vimt3d sm3d vaam vapm mdpm mrfr vsml vpdfl )

