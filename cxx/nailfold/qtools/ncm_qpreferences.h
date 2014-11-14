#ifndef __ncm_qpreferences_h__
#define __ncm_qpreferences_h__

//:
// \file
// \brief Structure to store the preferences of the system
//
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

#include "ncm_qglobal.h"

class ncm_qpreferences
{

//  INTERFACE

public:
  //  No member variables here, please
  ncm_qpreferences();

  void set_rezoom_factor(int rezoom_factor);
  int rezoom_factor() const;

  void set_vessel_minimum_interpoint_distance(double distance);
  double vessel_minimum_interpoint_distance() const;
  void set_vessel_maximum_interpoint_distance(double distance);
  double vessel_maximum_interpoint_distance() const;


//  IMPLEMENTATION

protected:

private:

  // Put preferences here

  //: Interact with Windows registry (where applicable)
  //void read_from_registry();
  //void write_to_registry();

  //: Appearance of a vessel placeholder
  QRgb placeholder_colour_unselected_;
  QRgb placeholder_colour_selected_;
  int placeholder_opacity_unselected_;
  int placeholder_opacity_selected_;

  //: Appearance of a vessel path
  QRgb vessel_colour_;
  QRgb arterial_endpoint_colour_;
  QRgb venous_endpoint_colour_;

  //: Appearance of a vessel apex
  QRgb apex_colour_selected_;
  QRgb apex_colour_unselected_;
  int apex_opacity_selected_;
  int apex_opacity_unselected_;

  //: Apex labelling parameters
  //: By how much to zoom in when centring on a vessel (-1 = don't rezoom)
  int rezoom_factor_;
};

#endif // __ncm_qpreferences_h__