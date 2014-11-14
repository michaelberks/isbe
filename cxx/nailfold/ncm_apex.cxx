#include "nailfold/ncm_apex.h"

#include <vgl/vgl_vector_2d.h>
#include <vgl/vgl_point_2d.h>

#include "ncm_vessel.h"

//: Define class member functions

//: Constructors
ncm_apex::ncm_apex(ncm_vessel& parent_vessel,
                   double x, double y)
: vessel_(parent_vessel)
{
  set_inner_point(x, y);
  set_outer_point(x, y);
}

//
//: Parent vessel
ncm_vessel& ncm_apex::parent_vessel()
{
  return vessel_;
}
ncm_vessel const& ncm_apex::parent_vessel() const
{
  return vessel_;
}

//
//: Uncertainty in the apex position (needs work)
double ncm_apex::uncertainty() const
{
  return 0.0;
}

//
//: Check that apex properties are valid (e.g. that it crosses the vessel)
bool ncm_apex::is_valid() const
{
  return true;
}

//
//: Coordinates of the apex centre
vgl_point_2d<double> ncm_apex::centre_point() const
{
  return vessel_.anchor() + centre_offset_;
}

//
//: Set coordinates of the inner point of the apex
void ncm_apex::set_inner_point(double x, double y)
{
  inner_offset_ = vgl_point_2d<double>(x, y) - vessel_.anchor();
}

//
//: Coordinates of the inner point of the apex
vgl_point_2d<double> ncm_apex::inner_point() const
{
  return vessel_.anchor() + inner_offset_;
}

//
//: Set coordinates of the inner point of the apex
void ncm_apex::set_outer_point(double x, double y)
{
  outer_offset_ = vgl_point_2d<double>(x, y) - vessel_.anchor();
}

//
//: Coordinates of the outer point of the apex
vgl_point_2d<double> ncm_apex::outer_point() const
{
  return vessel_.anchor() + outer_offset_;
}

//
//: Width of the vessel at the apex
double ncm_apex::width() const
{
  return (inner_offset_ - outer_offset_).length();
}