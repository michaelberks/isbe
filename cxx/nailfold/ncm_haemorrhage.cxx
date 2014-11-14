#include "ncm_haemorrhage.h"

#include <vgl/vgl_vector_2d.h>
#include <vgl/vgl_point_2d.h>

#include <nailfold/ncm_annotation.h>
#include <nailfold/ncm_vessel.h>

//: Define class member functions

//: Constructors
ncm_haemorrhage::ncm_haemorrhage(ncm_annotation& annotation,
                                 double x, double y,
                                 ncm_vessel* parent_vessel /* = NULL */)
: annotation_(annotation),
  anchor_(vgl_point_2d<double>(x,y)),
  vessel_(parent_vessel),
  outline_points_(0)
{
}

void ncm_haemorrhage::clear()
{
  anchor_.set(0, 0);
  vessel_ = NULL;
  outline_points_.resize(0);
}

//
//: Parent vessel
ncm_vessel const* ncm_haemorrhage::parent_vessel()
{
  return vessel_;
}

//
//: Reference point of the haemorrhage
void ncm_haemorrhage::set_anchor(double x, double y)
{
  if (x != anchor_.x() || y != anchor_.y())
  {
    anchor_.set(x, y);
    annotation_.set_modified();
  }
}
const vgl_point_2d<double>& ncm_haemorrhage::anchor() const
{
  return anchor_;
}

//
//: Add a point to the outline
void ncm_haemorrhage::add_outline_point(const vgl_point_2d<double>& point)
{
  annotation_.set_modified();

  outline_points_.push_back(point);
}
void ncm_haemorrhage::add_outline_point(double x, double y)
{
  add_outline_point(vgl_point_2d<double>(x,y));
}


//
//: Number of points in outline
unsigned ncm_haemorrhage::n_points() const
{
  return outline_points_.size();
}

//
//: Return i'th outline point
const vgl_point_2d<double>* ncm_haemorrhage::outline_point(unsigned i) const
{
  if (i < n_points())
    return &outline_points_[i];
  else
    return NULL;
}

//
//: Whether the outline is defined
bool ncm_haemorrhage::outline_is_defined() const
{
  return (n_points() > 0);
}

//
//: Delete all points from the outline
void ncm_haemorrhage::delete_outline()
{
  outline_points_.clear();
}

//
//: Attach to a specific vessel (possibly NULL)
void ncm_haemorrhage::attach_to_vessel(ncm_vessel const* vessel)
{
  vessel_ = vessel;
}