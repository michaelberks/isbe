//:
// \file
// \brief QGraphicsItem-derived class for creating a visualization of the two
//        points either side of the apex (during annotation).
// \author Phil Tresadern

#include "ncm_qapexitem_appearance.h"

#include <nailfold/qtools/ncm_qglobal.h>

//: Constructor
ncm_qapexitem_appearance::ncm_qapexitem_appearance()
: show_apex_(true),
  line_width_(0.004),
  width_is_relative_to_image_size_(true),
  colour_(NcmQt::green),
  colour_can_move_(NcmQt::yellow),
  colour_can_delete_(NcmQt::blue),
  opacity_(255),
  opacity_selected_(128)
{
}

void ncm_qapexitem_appearance::set_show_apex(bool show_apex)
{
  show_apex_ = show_apex;
}
bool ncm_qapexitem_appearance::apex_is_visible() const
{
  return show_apex_;
}


//
//: Line width
void ncm_qapexitem_appearance::set_line_width(double line_width)
{
  line_width_ = line_width;
}
double ncm_qapexitem_appearance::line_width() const
{
  return line_width_;
}
void ncm_qapexitem_appearance::set_width_relative_to_image_size(
    bool relative /* = true */)
{
  width_is_relative_to_image_size_ = relative;
}
bool ncm_qapexitem_appearance::width_is_relative_to_image_size()
{
  return width_is_relative_to_image_size_;
}

//
//: Colours
void ncm_qapexitem_appearance::set_colour(QRgb colour)
{
  colour_ = colour;
}
QRgb ncm_qapexitem_appearance::colour() const
{
  return colour_;
}
void ncm_qapexitem_appearance::set_colour_can_move(QRgb colour)
{
  colour_can_move_ = colour;
}
QRgb ncm_qapexitem_appearance::colour_can_move() const
{
  return colour_can_move_;
}
void ncm_qapexitem_appearance::set_colour_can_delete(QRgb colour)
{
  colour_can_delete_ = colour;
}
QRgb ncm_qapexitem_appearance::colour_can_delete() const
{
  return colour_can_delete_;
}

//
//: Opacities
void ncm_qapexitem_appearance::set_opacity(int opacity)
{
  opacity_ = opacity;
}
int ncm_qapexitem_appearance::opacity() const
{
  return opacity_;
}
void ncm_qapexitem_appearance::set_opacity_selected(int opacity)
{
  opacity_selected_ = opacity;
}
int ncm_qapexitem_appearance::opacity_selected() const
{
  return opacity_selected_;
}
