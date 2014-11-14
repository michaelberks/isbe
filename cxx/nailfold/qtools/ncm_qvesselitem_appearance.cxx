#include "ncm_qvesselitem_appearance.h"
#include "ncm_qglobal.h" // for colour definitions

//
//  Define class member functions

//: Constructor
ncm_qvesselitem_appearance::ncm_qvesselitem_appearance()
: show_placeholder_(false),
  show_path_(true),
  placeholder_radius_(0.02),
  line_width_(0.002),
  enlarged_scale_(2),
  giant_scale_(2),
  colour_selected_(NcmQt::yellow), // Selected
  colour_undefined_(NcmQt::white), // Undefined
  colour_normal_(NcmQt::red), // Normal
  colour_enlarged_(NcmQt::yellow), // Size: Enlarged
  colour_giant_(NcmQt::cyan), // Size: Giant
  colour_tortuous_(NcmQt::yellow), // Shape: Tortuous
  colour_ramified_(NcmQt::cyan), // Shape: Ramified
	colour_auto_distal_(NcmQt::red), // Auto detected vessel
	colour_auto_nondistal_(NcmQt::cyan), // Auto detected vessel
  opacity_(255),
  opacity_selected_(64),
  outline_when_selected_(true),
  warn_on_overwrite_size(true),
  warn_on_overwrite_shape(true)
{
}

//
//: Choose which elements of the vessel to show
void ncm_qvesselitem_appearance::set_show_placeholder(bool show_placeholder)
{
  show_placeholder_ = show_placeholder;
}
bool ncm_qvesselitem_appearance::placeholder_is_visible() const
{
  return show_placeholder_;
}
void ncm_qvesselitem_appearance::set_show_path(bool show_path)
{
  show_path_ = show_path;
}
bool ncm_qvesselitem_appearance::path_is_visible() const
{
  return show_path_;
}

//
//: Placeholder radius
void ncm_qvesselitem_appearance::set_placeholder_radius(double placeholder_radius)
{
  placeholder_radius_ = placeholder_radius;
}
double ncm_qvesselitem_appearance::placeholder_radius() const
{
  return placeholder_radius_;
}

//
//: Line widths
void ncm_qvesselitem_appearance::set_line_width(double line_width)
{
  line_width_ = line_width;
}
double ncm_qvesselitem_appearance::line_width() const
{
  return line_width_;
}
void ncm_qvesselitem_appearance::set_scale_enlarged(double enlarged_scale)
{
  enlarged_scale_ = enlarged_scale;
}
double ncm_qvesselitem_appearance::scale_enlarged() const
{
  return enlarged_scale_;
}
void ncm_qvesselitem_appearance::set_scale_giant(double giant_scale)
{
  giant_scale_ = giant_scale;
}
double ncm_qvesselitem_appearance::scale_giant() const
{
  return giant_scale_;
}

//
//: Colours
void ncm_qvesselitem_appearance::set_colour_selected(QRgb vessel_colour)
{
  colour_selected_ = vessel_colour;
}
QRgb ncm_qvesselitem_appearance::colour_selected() const
{
  return colour_selected_;
}
void ncm_qvesselitem_appearance::set_colour_undefined(QRgb vessel_colour)
{
  colour_undefined_ = vessel_colour;
}
QRgb ncm_qvesselitem_appearance::colour_undefined() const
{
  return colour_undefined_;
}


void ncm_qvesselitem_appearance::set_colour_normal(QRgb vessel_colour)
{
  colour_normal_ = vessel_colour;
}
QRgb ncm_qvesselitem_appearance::colour_normal() const
{
  return colour_normal_;
}
void ncm_qvesselitem_appearance::set_colour_enlarged(QRgb vessel_colour)
{
  colour_enlarged_ = vessel_colour;
}
QRgb ncm_qvesselitem_appearance::colour_enlarged() const
{
  return colour_enlarged_;
}
void ncm_qvesselitem_appearance::set_colour_giant(QRgb vessel_colour)
{
  colour_giant_ = vessel_colour;
}
QRgb ncm_qvesselitem_appearance::colour_giant() const
{
  return colour_giant_;
}
void ncm_qvesselitem_appearance::set_colour_tortuous(QRgb vessel_colour)
{
  colour_tortuous_ = vessel_colour;
}
QRgb ncm_qvesselitem_appearance::colour_tortuous() const
{
  return colour_tortuous_;
}
void ncm_qvesselitem_appearance::set_colour_ramified(QRgb vessel_colour)
{
  colour_ramified_ = vessel_colour;
}
QRgb ncm_qvesselitem_appearance::colour_ramified() const
{
  return colour_giant_;
}

void ncm_qvesselitem_appearance::set_colour_auto_distal(QRgb vessel_colour)
{
  colour_auto_distal_ = vessel_colour;
}
QRgb ncm_qvesselitem_appearance::colour_auto_distal() const
{
  return colour_auto_distal_;
}
void ncm_qvesselitem_appearance::set_colour_auto_nondistal(QRgb vessel_colour)
{
  colour_auto_nondistal_ = vessel_colour;
}
QRgb ncm_qvesselitem_appearance::colour_auto_nondistal() const
{
  return colour_auto_nondistal_;
}

//
//: Opacities
void ncm_qvesselitem_appearance::set_opacity(int opacity)
{
  opacity_ = opacity;
}
int ncm_qvesselitem_appearance::opacity() const
{
  return opacity_;
}
void ncm_qvesselitem_appearance::set_opacity_selected(int opacity)
{
  opacity_selected_ = opacity;
}
int ncm_qvesselitem_appearance::opacity_selected() const
{
  return opacity_selected_;
}

//
//: Outline only mode
void ncm_qvesselitem_appearance::set_outline_when_selected(bool outline_when_selected)
{
  outline_when_selected_ = outline_when_selected;
}
bool ncm_qvesselitem_appearance::outline_when_selected() const
{
  return outline_when_selected_;
}
