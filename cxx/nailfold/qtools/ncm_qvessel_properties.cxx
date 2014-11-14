#include "ncm_qvessel_properties.h"

//
//  Define class member functions

//: Constructor
ncm_qvessel_properties::ncm_qvessel_properties()
{
}

//: Set distal/nondistal
void ncm_qvessel_properties::setDistal(bool)
{
  set_distal(true);
}
void ncm_qvessel_properties::setNondistal(bool)
{
  set_distal(false);
}

//: Set size properties
void ncm_qvessel_properties::setSizeUndefined(bool)
{
  set_size_undefined();
}
void ncm_qvessel_properties::setSizeNormal(bool)
{
  set_size_normal();
}
void ncm_qvessel_properties::setSizeEnlarged(bool)
{
  set_size_enlarged();
}
void ncm_qvessel_properties::setSizeGiant(bool)
{
  set_size_giant();
}
void ncm_qvessel_properties::setSizeIrregular(bool)
{
  set_size_irregular();
}

//: Set shape properties
void ncm_qvessel_properties::setShapeUndefined(bool)
{
  set_shape_undefined();
}
void ncm_qvessel_properties::setShapeNormal(bool)
{
  set_shape_normal();
}
void ncm_qvessel_properties::setShapeTortuous(bool is_tortuous /* = true */)
{
  // Set tortuous flag to value defined by is_tortuous, leaving all other
  // flags intact
  set_shape_tortuous(is_tortuous);
}
void ncm_qvessel_properties::setShapeRamified(bool is_ramified /* = true */)
{
  // Set ramified flag to value defined by is_ramified, leaving all other
  // flags intact
  set_shape_ramified(is_ramified);
}
void ncm_qvessel_properties::setShapeTortuousOnly(bool)
{
  // Set tortuous flag to true; reset all other flags
  set_shape_undefined();
  set_shape_tortuous();
}
void ncm_qvessel_properties::setShapeRamifiedOnly(bool)
{
  // Set ramified flag to true; reset all other flags
  set_shape_undefined();
  set_shape_ramified();
}
