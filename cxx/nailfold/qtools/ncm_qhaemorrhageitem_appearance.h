#ifndef __ncm_qhaemorrhageitem_appearance_h__
#define __ncm_qhaemorrhageitem_appearance_h__

#include <QRgb>

//:
// \file
// \brief Class that holds the appearance variables for drawing a vessel in a
//        scene.
//
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

class ncm_qhaemorrhageitem_appearance 
{
//  INTERFACE

public:
  //  No member variables here, please
  ncm_qhaemorrhageitem_appearance();

  //: Choose which elements of the haemorrhage to show
  void set_show_placeholder(bool show_placeholder);
  bool placeholder_is_visible() const;
  void set_show_outline(bool show_outline);
  bool outline_is_visible() const;

  //: Set/get appearance properties of all vessels as shown onscreen
  void set_placeholder_radius(double placeholder_radius);
  double placeholder_radius() const;

  void set_line_width(double line_width);
  double line_width() const;
  
  void set_colour(QRgb vessel_colour);
  QRgb colour() const;
  void set_colour_selected(QRgb vessel_colour);
  QRgb colour_selected() const;
  
  void set_opacity(int opacity);
  int opacity() const;
  void set_opacity_selected(int opacity);
  int opacity_selected() const;
  
  void set_outline_when_selected(bool outline_when_selected);
  bool outline_when_selected() const;


//  IMPLEMENTATION

protected:

private:

  // Drawing properties should be the same for all vessels so make them
  // static
  // - radius of the placeholder
  // - width
  // - colour (arterial, venous)
  // - opacity (active, inactive)
  // - whether to draw only the outline when selected (rendering it almost
  //   invisible

  bool show_placeholder_;
  bool show_outline_;

  double placeholder_radius_;

  double line_width_;

  QRgb colour_;
  QRgb colour_selected_;

  int opacity_;
  int opacity_selected_;

  bool outline_when_selected_;
};

#endif // __ncm_qhaemorrhageitem_appearance_h__