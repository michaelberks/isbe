#include "ncm_qhaemorrhageitem_appearance.h"

#include <nailfold/qtools/ncm_qglobal.h> // for colour definitions

//
//  Define class member functions

//: Constructor
ncm_qhaemorrhageitem_appearance::ncm_qhaemorrhageitem_appearance()
: show_placeholder_(false),
  show_outline_(true),
  placeholder_radius_(0.02),
  line_width_(0.004),
  colour_(NcmQt::red),
  colour_selected_(NcmQt::yellow),
  opacity_(255),
  opacity_selected_(128),
  outline_when_selected_(false)
{
}

//
//: Choose which elements of the vessel to show
void ncm_qhaemorrhageitem_appearance::set_show_placeholder(bool show_placeholder)
{
  show_placeholder_ = show_placeholder;
}
bool ncm_qhaemorrhageitem_appearance::placeholder_is_visible() const
{
  return show_placeholder_;
}
void ncm_qhaemorrhageitem_appearance::set_show_outline(bool show_outline)
{
  show_outline_ = show_outline;
}
bool ncm_qhaemorrhageitem_appearance::outline_is_visible() const
{
  return show_outline_;
}


//
//: Placeholder radius
void ncm_qhaemorrhageitem_appearance::set_placeholder_radius(double placeholder_radius)
{
  placeholder_radius_ = placeholder_radius;
}
double ncm_qhaemorrhageitem_appearance::placeholder_radius() const
{
  return placeholder_radius_;
}

//
//: Line widths
void ncm_qhaemorrhageitem_appearance::set_line_width(double line_width)
{
  line_width_ = line_width;
}
double ncm_qhaemorrhageitem_appearance::line_width() const
{
  return line_width_;
}

//
//: Colours
void ncm_qhaemorrhageitem_appearance::set_colour(QRgb vessel_colour)
{
  colour_ = vessel_colour;
}
QRgb ncm_qhaemorrhageitem_appearance::colour() const
{
  return colour_;
}
void ncm_qhaemorrhageitem_appearance::set_colour_selected(QRgb vessel_colour)
{
  colour_selected_ = vessel_colour;
}
QRgb ncm_qhaemorrhageitem_appearance::colour_selected() const
{
  return colour_selected_;
}

//
//: Opacities
void ncm_qhaemorrhageitem_appearance::set_opacity(int opacity)
{
  opacity_ = opacity;
}
int ncm_qhaemorrhageitem_appearance::opacity() const
{
  return opacity_;
}
void ncm_qhaemorrhageitem_appearance::set_opacity_selected(int opacity)
{
  opacity_selected_ = opacity;
}
int ncm_qhaemorrhageitem_appearance::opacity_selected() const
{
  return opacity_selected_;
}

//
//: Outline only mode
void ncm_qhaemorrhageitem_appearance::set_outline_when_selected(bool outline_when_selected)
{
  outline_when_selected_ = outline_when_selected;
}
bool ncm_qhaemorrhageitem_appearance::outline_when_selected() const
{
  return outline_when_selected_;
}