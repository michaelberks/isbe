#ifndef __ncm_qapexitem_appearance_h__
#define __ncm_qapexitem_appearance_h__

#include <QRgb>

//:
// \file
// \brief Class that holds the appearance variables for drawing an apex in a
//        scene.
//
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

class ncm_qapexitem_appearance 
{
//  INTERFACE

public:
  //  No member variables here, please
  ncm_qapexitem_appearance();

  void set_show_apex(bool show_apex);
  bool apex_is_visible() const;

  //: Set/get appearance properties of all vessels as shown onscreen
  void set_line_width(double line_width);
  double line_width() const;
  void set_width_relative_to_image_size(bool relative = true);
  bool width_is_relative_to_image_size();

  void set_colour(QRgb colour);
  QRgb colour() const;
  void set_colour_can_move(QRgb colour);
  QRgb colour_can_move() const;
  void set_colour_can_delete(QRgb colour);
  QRgb colour_can_delete() const;
  
  void set_opacity(int opacity);
  int opacity() const;
  void set_opacity_selected(int opacity);
  int opacity_selected() const;
  

//  IMPLEMENTATION

protected:

private:

  // Drawing properties for apices
  // - width
  // - colour (arterial, venous)
  // - opacity (active, inactive)

  bool show_apex_;

  double line_width_;
  bool width_is_relative_to_image_size_;

  QRgb colour_;
  QRgb colour_can_move_;
  QRgb colour_can_delete_;

  int opacity_;
  int opacity_selected_;
};

#endif // __ncm_qapexitem_appearance_h__