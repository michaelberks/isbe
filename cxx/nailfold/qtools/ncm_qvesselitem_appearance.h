#ifndef __ncm_qvesselitem_appearance_h__
#define __ncm_qvesselitem_appearance_h__

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

class ncm_qvesselitem_appearance 
{
//  INTERFACE

public:
  //  No member variables here, please
  ncm_qvesselitem_appearance();

  //: Choose which elements to display
  void set_show_placeholder(bool show_placeholder);
  bool placeholder_is_visible() const;
  void set_show_path(bool show_path);
  bool path_is_visible() const;

  //: Set/get appearance properties of all vessels as shown onscreen
  void set_placeholder_radius(double placeholder_radius);
  double placeholder_radius() const;

  void set_line_width(double line_width);
  double line_width() const;
  void set_scale_enlarged(double enlarged_scale);
  double scale_enlarged() const;
  void set_scale_giant(double giant_scale);
  double scale_giant() const;
  
  void set_colour_selected(QRgb vessel_colour);
  QRgb colour_selected() const;
  void set_colour_undefined(QRgb vessel_colour);
  QRgb colour_undefined() const;
  void set_colour_normal(QRgb vessel_colour);
  QRgb colour_normal() const;

  void set_colour_enlarged(QRgb vessel_colour);
  QRgb colour_enlarged() const;
  void set_colour_giant(QRgb vessel_colour);
  QRgb colour_giant() const;
  
  void set_colour_tortuous(QRgb vessel_colour);
  QRgb colour_tortuous() const;
  void set_colour_ramified(QRgb vessel_colour);
  QRgb colour_ramified() const;

	void set_colour_auto_distal(QRgb vessel_colour);
  QRgb colour_auto_distal() const;
  void set_colour_auto_nondistal(QRgb vessel_colour);
  QRgb colour_auto_nondistal() const;

  void set_opacity(int opacity);
  int opacity() const;
  void set_opacity_selected(int opacity);
  int opacity_selected() const;
  
  void set_outline_when_selected(bool outline_when_selected);
  bool outline_when_selected() const;

  bool warn_on_overwrite_size;
  bool warn_on_overwrite_shape;


//  IMPLEMENTATION

protected:

private:

  // Drawing properties should be the same for all vessels so make them
  // static
  // - whether to display different elements of the vessel
  // - radius of the placeholder
  // - width (normal, enlarged, giant)
  // - colour (arterial, venous)
  // - opacity (active, inactive)
  // - whether to draw only the outline when selected (rendering it almost
  //   invisible

  bool show_placeholder_;
  bool show_path_;

  double placeholder_radius_;

  double line_width_;
  double enlarged_scale_;
  double giant_scale_;

  QRgb colour_selected_;
  QRgb colour_undefined_;
  QRgb colour_normal_;
  
  QRgb colour_enlarged_;
  QRgb colour_giant_;

  QRgb colour_tortuous_;
  QRgb colour_ramified_;

	QRgb colour_auto_distal_;
	QRgb colour_auto_nondistal_;

  int opacity_;
  int opacity_selected_;

  bool outline_when_selected_;
};

#endif // __ncm_qvesselitem_appearance_h__