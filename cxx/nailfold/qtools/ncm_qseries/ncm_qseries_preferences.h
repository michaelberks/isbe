#ifndef __ncm_qseries_preferences_h__
#define __ncm_qseries_preferences_h__

//:
// \file
// \brief Structure to store the preferences of the system
//
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_string.h>

#include <nailfold/qtools/ncm_qglobal.h>
#include <nailfold/qtools/ncm_qvesselitem_appearance.h>

// forward declaration
class ncm_qseries_options;

class ncm_qseries_preferences
{
  friend ncm_qseries_options;

//  INTERFACE

public:
  //  No member variables here, please

  enum ImageZoom { ZoomNoResize = 0,
                   ZoomFit,
                   ZoomFitWidth,
                   ZoomFitHeight,
  
                   ZoomFirst = ZoomNoResize,
                   ZoomLast = ZoomFitHeight };

  //: Default constructor
  ncm_qseries_preferences();

  //: Destructor
  ~ncm_qseries_preferences();

  //: Write preferences to Windows registry
  void write_to_registry();

  //: Read preferences from Windows registry
  void read_from_registry();

  double version() const;

  void set_username(const vcl_string& username);
  vcl_string username() const;
  void set_advanced_user(bool is_advanced_user = true);
  bool is_advanced_user() const;

  void set_grid_pixels_per_mm(double grid_pixels_per_mm);
  double grid_pixels_per_mm() const;
  void set_grid_spacing_mm(double grid_spacing_mm);
  double grid_spacing_mm() const;

  //: Access functions
  ImageZoom label_image_zoom() const;
  ImageZoom add_vessels_zoom() const;
  ImageZoom set_vessel_props_zoom() const;
  ImageZoom add_haemorrhages_zoom() const;

  int normal_vessel_zoom() const;
  int enlarged_zoom_relative() const;
  int giant_zoom_relative() const;
  double vessel_path_yshift() const;

  void set_use_remote_server(bool use_remote_server = true);
  bool use_remote_server();

  //: Enable/disable warnings
  bool warn_if_ungraded;
  bool warn_if_vessels_missing;
  bool warn_if_sizes_missing;
  bool warn_if_shapes_missing;
  bool warn_if_apices_missing;
  bool warn_if_paths_missing;
  bool warn_if_haemorrhages_missing;
	bool warn_if_no_reasons;

  
//  IMPLEMENTATION

protected:

private:

  // Put preferences here

  vcl_string system_username() const;

  //: Apex labelling parameters
  //: By how much to zoom in when centring on a vessel (-1 = don't rezoom)
  //int rezoom_factor_;

  //: Current software version
  static const double version_;

  ImageZoom label_image_zoom_;
  ImageZoom add_vessels_zoom_;
  ImageZoom set_vessel_props_zoom_;
  ImageZoom add_haemorrhages_zoom_;

  double grid_pixels_per_mm_;
  double grid_spacing_mm_;

  int normal_vessel_zoom_;
  int enlarged_zoom_relative_;
  int giant_zoom_relative_;

  double vessel_path_yshift_;

  vcl_string username_;
  bool is_advanced_user_;

  bool use_remote_server_;
};

#endif // __ncm_qseries_preferences_h__