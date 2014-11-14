#include "ncm_qfinal_preferences.h"

#include <vul/vul_string.h>

#include <nailfold/ncm_vessel.h>
#include <nailfold/ncm_registry_accessor.h>

#include <QString>

//
// Static member variables
//

const double ncm_qfinal_preferences::version_ = 1.5;

//
//  Define class member functions

//: Constructor
ncm_qfinal_preferences::ncm_qfinal_preferences()
: label_image_zoom_(ZoomFit),
  add_vessels_zoom_(ZoomFitHeight),
  set_vessel_props_zoom_(ZoomFit),
  add_haemorrhages_zoom_(ZoomFit),
  normal_vessel_zoom_(200.0),
  enlarged_zoom_relative_(1),
  giant_zoom_relative_(1),
  vessel_path_yshift_(0.6),
  grid_pixels_per_mm_(803.85),
  grid_spacing_mm_(1.0),
  username_(system_username()),
  is_advanced_user_(false),
  warn_if_ungraded(true),
  warn_if_vessels_missing(true),
  warn_if_sizes_missing(true),
  warn_if_shapes_missing(true),
  warn_if_apices_missing(true),
  warn_if_paths_missing(false),
  warn_if_haemorrhages_missing(false),
  use_remote_server_(false)
{
  // Read values from the registry if they exist; 
  // (Write default values to the registry if they don't.)
  read_from_registry();
}

//: Destructor
ncm_qfinal_preferences::~ncm_qfinal_preferences()
{
  // When the preferences object dies, save its values in the registry.
  write_to_registry();
}

//  Access functions


double ncm_qfinal_preferences::version() const
{
  return version_;
}

ncm_qfinal_preferences::ImageZoom ncm_qfinal_preferences::label_image_zoom() const
{
  return label_image_zoom_;
}
ncm_qfinal_preferences::ImageZoom ncm_qfinal_preferences::add_vessels_zoom() const
{
  return add_vessels_zoom_;
}
ncm_qfinal_preferences::ImageZoom ncm_qfinal_preferences::set_vessel_props_zoom() const
{
  return set_vessel_props_zoom_;
}
ncm_qfinal_preferences::ImageZoom ncm_qfinal_preferences::add_haemorrhages_zoom() const
{
  return add_haemorrhages_zoom_;
}

int ncm_qfinal_preferences::normal_vessel_zoom() const
{
  return normal_vessel_zoom_;
}
int ncm_qfinal_preferences::enlarged_zoom_relative() const
{
  return enlarged_zoom_relative_;
}
int ncm_qfinal_preferences::giant_zoom_relative() const
{
  return giant_zoom_relative_;
}
double ncm_qfinal_preferences::vessel_path_yshift() const
{
  return vessel_path_yshift_;
}

void ncm_qfinal_preferences::set_use_remote_server(
    bool use_remote_server /* = true */)
{
  use_remote_server_ = use_remote_server;
}

bool ncm_qfinal_preferences::use_remote_server()
{
  return use_remote_server_;
}

//
//: Write preferences to Windows registry
void ncm_qfinal_preferences::write_to_registry()
{
  ncm_registry_accessor ra;
  ra.get_current_user_key();

  long result = 
      ra.get_subkey("Software\\University of Manchester\\NCM QMarkup");

  if (result != ERROR_SUCCESS)
    ra.create_subkey("Software\\University of Manchester\\NCM QMarkup");

  // Write name/value pairs to registry
  QString s;

  // Always save the version number with the preferences
  s = QString::number(version());
  ra.write_string("Version", s.toStdString());

  s = QString::number(normal_vessel_zoom_);
  ra.write_string("NormalVesselZoom", s.toStdString());
  s = QString::number(vessel_path_yshift_);
  ra.write_string("VesselPathYShift", s.toStdString());
  s = QString::number(grid_pixels_per_mm_);
  ra.write_string("GridPixelsPerMm", s.toStdString());
  s = QString::number(grid_spacing_mm_);
  ra.write_string("GridSpacingMm", s.toStdString());

  ra.write_numeric("LabelImageZoom", label_image_zoom_);
  ra.write_numeric("AddVesselsZoom", add_vessels_zoom_);
  ra.write_numeric("SetVesselPropsZoom", set_vessel_props_zoom_);
  ra.write_numeric("AddHaemorrhagesZoom", add_haemorrhages_zoom_);

  ra.write_numeric("EnlargedZoomRelative", enlarged_zoom_relative_);
  ra.write_numeric("GiantZoomRelative", giant_zoom_relative_);

  ra.write_boolean("IsAdvancedUser", is_advanced_user_);
  ra.write_boolean("WarnIfUngraded", warn_if_ungraded);
  ra.write_boolean("WarnIfVesselsMissing", warn_if_vessels_missing);
  ra.write_boolean("WarnIfSizesMissing", warn_if_sizes_missing);
  ra.write_boolean("WarnIfShapesMissing", warn_if_shapes_missing);
  ra.write_boolean("WarnIfApicesMissing", warn_if_apices_missing);
  ra.write_boolean("WarnIfPathsMissing", warn_if_paths_missing);
  ra.write_boolean("WarnIfHaemorrhagesMissing", warn_if_haemorrhages_missing);
}

//
//: Read preferences from Windows registry
void ncm_qfinal_preferences::read_from_registry()
{
  ncm_registry_accessor ra;
  ra.get_current_user_key();

  long result = 
      ra.get_subkey("Software\\University of Manchester\\NCM QMarkup");

  if (result != ERROR_SUCCESS)
  {
    // Write defaults to the registry and return if preferences have not been 
    // stored previously. 
    write_to_registry();
    return;
  }

  // Find out which software version was used to write the preferences and read
  // them back in an appropriate way.
  vcl_string s = "";
  ra.read_string("Version", s);
  double prefs_version = vul_string_atof(s);

  if (prefs_version == 1.3)
  {
    // Read name/value pairs to registry
    ra.read_string("NormalVesselZoom", s);
    normal_vessel_zoom_ = vul_string_atof(s);
    ra.read_string("VesselPathYShift", s);
    vessel_path_yshift_ = vul_string_atof(s);
    ra.read_string("GridPixelsPerMm", s);
    grid_pixels_per_mm_ = vul_string_atof(s);
    ra.read_string("GridSpacingMm", s);
    grid_spacing_mm_ = vul_string_atof(s);

    // Need to cast from int to ncm_qfinal_preferences::ImageZoom explicitly
    int zoom_type;
    ra.read_numeric("LabelImageZoom", zoom_type);
    label_image_zoom_ = 
        static_cast<ncm_qfinal_preferences::ImageZoom>(zoom_type);
    ra.read_numeric("AddVesselsZoom", zoom_type);
    add_vessels_zoom_ = 
        static_cast<ncm_qfinal_preferences::ImageZoom>(zoom_type);
    ra.read_numeric("SetVesselPropsZoom", zoom_type);
    set_vessel_props_zoom_ = 
        static_cast<ncm_qfinal_preferences::ImageZoom>(zoom_type);
    ra.read_numeric("AddHaemorrhagesZoom", zoom_type);
    add_haemorrhages_zoom_ = 
        static_cast<ncm_qfinal_preferences::ImageZoom>(zoom_type);

    ra.read_numeric("EnlargedZoomRelative", enlarged_zoom_relative_);
    ra.read_numeric("GiantZoomRelative", giant_zoom_relative_);

    ra.read_boolean("IsAdvancedUser", is_advanced_user_);
    ra.read_boolean("WarnIfUngraded", warn_if_ungraded);
    ra.read_boolean("WarnIfVesselsMissing", warn_if_vessels_missing);
    ra.read_boolean("WarnIfSizesMissing", warn_if_sizes_missing);
    ra.read_boolean("WarnIfShapesMissing", warn_if_shapes_missing);
    ra.read_boolean("WarnIfApicesMissing", warn_if_apices_missing);
    ra.read_boolean("WarnIfPathsMissing", warn_if_paths_missing);
    ra.read_boolean("WarnIfHaemorrhagesMissing", warn_if_haemorrhages_missing);
  }
}

//
//: User properties
void ncm_qfinal_preferences::set_username(const vcl_string& username)
{
  username_ = username;
}
vcl_string ncm_qfinal_preferences::username() const
{
  return username_;
}
void ncm_qfinal_preferences::set_advanced_user(
    bool is_advanced_user /* = true */)
{
  is_advanced_user_ = is_advanced_user;
}
bool ncm_qfinal_preferences::is_advanced_user() const
{
  return is_advanced_user_;
}

vcl_string ncm_qfinal_preferences::system_username() const
{
#if defined(VCL_WIN32)
  if (getenv("USERNAME") != NULL)
    return vcl_string(getenv("USERNAME"));
  else
    return "";
#else
  return "";
#endif
}

void ncm_qfinal_preferences::set_grid_pixels_per_mm(double grid_pixels_per_mm)
{
  grid_pixels_per_mm_ = grid_pixels_per_mm;
}
double ncm_qfinal_preferences::grid_pixels_per_mm() const
{
  return grid_pixels_per_mm_;
}

void ncm_qfinal_preferences::set_grid_spacing_mm(double grid_spacing_mm)
{
  grid_spacing_mm_ = grid_spacing_mm;
}
double ncm_qfinal_preferences::grid_spacing_mm() const
{
  return grid_spacing_mm_;
}