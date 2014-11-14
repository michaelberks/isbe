#include "ncm_vessel.h"

#include <vcl_algorithm.h>

#include <vgl/vgl_vector_2d.h>

#include "ncm_annotation.h"
#include "ncm_apex.h"

//  Initialize static member variables
double ncm_vessel::minimum_inter_point_distance_ = 1.0;
double ncm_vessel::maximum_inter_point_distance_ = 99.0;
unsigned ncm_vessel::normal_smoothing_ = 1;

//
//:
void ncm_vessel::set_minimum_inter_point_distance(double distance)
{
  // throw an error if distance < 0?

  if (distance > 0.0)
    minimum_inter_point_distance_ = distance;
}

//
//:
double ncm_vessel::minimum_inter_point_distance() 
{
  return minimum_inter_point_distance_;
}

//
//:
void ncm_vessel::set_maximum_inter_point_distance(double distance)
{
  // throw an error if distance < minimum_inter_point_distance_ ?

  if (distance >= minimum_inter_point_distance_)
    maximum_inter_point_distance_ = distance;
}

//
//:
double ncm_vessel::maximum_inter_point_distance() 
{
  return maximum_inter_point_distance_;
}


// Define creator methods

//
//: Create a new vessel and return the pointer (to an ncm_annotation)
ncm_vessel* ncm_vessel_handler::new_vessel(
    const ncm_annotation& parent_annotation,
    const vgl_point_2d<double>& anchor)
{
  return new ncm_vessel(parent_annotation, anchor);
}

//
//: Create a new vessel and return the pointer (to an ncm_annotation)
ncm_vessel* ncm_vessel_handler::new_vessel(
    const ncm_annotation& parent_annotation,
    double anchor_x, double anchor_y)
{
  return new ncm_vessel(parent_annotation, anchor_x, anchor_y);
}

//
//: Destroy an existing vessel
void ncm_vessel_handler::delete_vessel(ncm_vessel* vessel)
{
  delete vessel;
}

//
//  Define class member functions

//: Constructors (every vessel must have at least one point)
ncm_vessel::ncm_vessel(const ncm_annotation& parent_annotation,
                       const vgl_point_2d<double>& anchor)
: annotation_(parent_annotation),
  anchor_(anchor),
  apices_(0),
  properties_(this)
{
  //apices_.resize(0);
}

ncm_vessel::ncm_vessel(const ncm_annotation& parent_annotation,
                       double anchor_x, double anchor_y)
: annotation_(parent_annotation),
  anchor_(vgl_point_2d<double>(anchor_x, anchor_y)),
  apices_(0),
  properties_(this)
{
}

//: Destructor
ncm_vessel::~ncm_vessel()
{
  delete_all_apices();
}

//
//: Reset vessel to blank
void ncm_vessel::clear()
{
  delete_all_apices();
  anchor_.set(0, 0);
  properties().clear();
}

void ncm_vessel::set_modified()
{
  annotation_.set_modified();
}

//
//: Get the index of this vessel in the corresponding annotation
int ncm_vessel::index() const
{
  return annotation_.index_of_vessel(this);
}

//
//: Set the anchor point of the vessel
void ncm_vessel::set_anchor_to(const vgl_point_2d<double>& anchor)
{
  if (anchor_ != anchor)
  {
    anchor_ = anchor;
    set_modified();
  }
}

//
//: Return the anchor point of the vessel
const vgl_point_2d<double>& ncm_vessel::anchor() const
{
  return anchor_;
}

//
//: Set vessel properties
ncm_vessel_properties& ncm_vessel::properties() 
{
  return properties_;
}
const ncm_vessel_properties& ncm_vessel::properties() const
{
  return properties_;
}

//
//: Whether path is defined
bool ncm_vessel::path_is_defined() const
{
  return (n_points() > 0);
}

//
//: Whether the path is 'suspicious' (i.e. the vessel anchor to which it is
//  closest is one other than itself).
bool ncm_vessel::path_is_suspicious() const
{
  double minimum_distance = 1e9;
  ncm_vessel const* nearest_vessel = NULL;

  for (unsigned p = 0; p < n_points(); ++p)
  {
    // Find vessel closest to this path point and note the distance
    double distance = 1e9;
    ncm_vessel const* nearby_vessel = 
        annotation_.vessel_nearest_to(*point(p), &distance);
  
    if (distance < minimum_distance)
    {
      nearest_vessel = nearby_vessel;
      minimum_distance = distance;
    }
  }

  // Return true if there is another vessel closer to the path than this one
  return (nearest_vessel != NULL && 
          nearest_vessel != this);
}

//
//: Delete all points from the vessel path
void ncm_vessel::delete_path()
{
  if (path_is_defined())
  {
    points_.clear();
    set_modified();
  }
}

//
//: Trim off the venous end of the vessel from the point closest to (x,y)
void ncm_vessel::trim_at(double x, double y)
{
  if (path_is_defined())
  {
    int cut_index = point_index_nearest_to(x,y);
    points_.resize(cut_index+1);
    set_modified();
  }
}

//
//: Add point to arterial limb
void ncm_vessel::add_arterial_point(const vgl_point_2d<double>& point,
                                    double width /* = 0.0 */)
{
  bool point_is_valid = false;

  if (path_is_defined())
  {
    const double distance_to_endpoint = (*arterial_endpoint() - point).length();
    const bool too_close =
        (distance_to_endpoint < minimum_inter_point_distance_);
    const bool too_far = 
        (distance_to_endpoint > maximum_inter_point_distance_);

    point_is_valid = (!too_close && !too_far);
  }
  else
  {
    point_is_valid = true;
  }

  if (point_is_valid)
  {
    points_.insert(points_.begin(), 1, point);
    widths_.insert(widths_.begin(), 1, width);

    // Annotation is not updated here as set_modified also writes to the backup
    // file, which has a noticeable effect when tracing out a vessel
    // Instead, set_modified() is called by the GUI when the mouse button is 
    // released
    //set_modified();
  }
}

//
//: Add point to arterial limb
void ncm_vessel::add_arterial_point(double x, double y,
                                    double width /* = 0.0 */)
{
  add_arterial_point(vgl_point_2d<double>(x, y), width);
}

//
//: Add point to venous limb
void ncm_vessel::add_venous_point(const vgl_point_2d<double>& point, 
                                  double width /* = 0.0 */)
{
  bool point_is_valid = false;

  if (path_is_defined())
  {
    const double distance_to_endpoint = (*venous_endpoint() - point).length();
    const bool too_close =
        (distance_to_endpoint < minimum_inter_point_distance_);
    const bool too_far = 
        (distance_to_endpoint > maximum_inter_point_distance_);

    point_is_valid = (!too_close && !too_far);
  }
  else
  {
    point_is_valid = true;
  }

  if (point_is_valid)
  {
    points_.push_back(point);
    widths_.push_back(width);

    // Annotation is not updated here as set_modified also writes to the backup
    // file, which has a noticeable effect when tracing out a vessel
    // Instead, set_modified() is called by the GUI when the mouse button is 
    // released
    //set_modified();
  }
}

//
//: Add point to venous limb
void ncm_vessel::add_venous_point(double x, double y,
                                  double width /* = 0.0 */)
{
  add_venous_point(vgl_point_2d<double>(x,y),width);
}


//
//: Return number of points in the path
unsigned ncm_vessel::n_points() const
{
  return points_.size();
}

//
//: Pointer to i'th point on the path
const vgl_point_2d<double>* ncm_vessel::point(unsigned i) const
{
  if (i < n_points())
    return &points_[i];
  else
    return NULL;
}

//
//: Return arterial endpoint
const vgl_point_2d<double>* ncm_vessel::arterial_endpoint() const
{
  return point(0);
}

//
//: Return venous endpoint
const vgl_point_2d<double>* ncm_vessel::venous_endpoint() const
{
  return point(n_points()-1);
}

//
//: Return the index of the point nearest to (x,y)
int ncm_vessel::point_index_nearest_to(double x, double y) const
{
  if (!path_is_defined())
    return -1;

  const unsigned n_points = points_.size();

  vgl_point_2d<double> ref_point(x,y);

  double min_dist = (ref_point - points_[0]).length();
  unsigned min_index = 0;

  for (unsigned i = 1; i < n_points; ++i)
  {
    double dist = (ref_point - points_[i]).length();
    if (dist < min_dist)
    {
      min_dist = dist;
      min_index = i;
    }
  }

  return min_index;
}

//
//: Return the point nearest to (x,y)
const vgl_point_2d<double>* ncm_vessel::point_nearest_to(double x, double y) const
{
  if (!path_is_defined())
    return NULL;

  return &points_[point_index_nearest_to(x, y)];
}

//
//: Normal to the vessel at a given point
vgl_vector_2d<double> ncm_vessel::normal_at(double x, double y) const
{
  if (!path_is_defined())
    return vgl_vector_2d<double>(0,0); // undefined directional vector

  assert( normal_smoothing_ > 0);

  const unsigned n_points = points_.size();
  const unsigned nearest_index = point_index_nearest_to(x, y);
  
  // index of first point to use
  unsigned first_index = 0;
  if (normal_smoothing_ <= nearest_index )
    first_index = nearest_index - normal_smoothing_;

  // index of last point to use
  unsigned last_index = n_points-1;
  if (nearest_index + normal_smoothing_ < n_points)
    last_index = nearest_index + normal_smoothing_;

  // vector (approximately) parallel to vessel
  vgl_vector_2d<double> v(0,0);
  for (unsigned i = first_index; i < last_index; ++i)
    v += (points_[i+1] - points_[i]);

  vgl_vector_2d<double> normal(-v.y(), v.x());
  return normalized(normal);
}

//
//: The width of the vessel at every point
const vcl_vector<double>& ncm_vessel::widths() const
{
  return widths_;
}

//
//: How many apices
unsigned ncm_vessel::n_apices() const
{
  return apices_.size();
}

//
//: The i'th apex
//  Returns NULL if i is out of range
const ncm_apex* ncm_vessel::apex(unsigned i) const
{
  if (i < n_apices())
    return apices_[i];
  else
    return NULL;
}

ncm_apex* ncm_vessel::apex(unsigned i)
{
  if (i < n_apices())
    return apices_[i];
  else
    return NULL;
}

//
//: Add new apex to this vessel, starting at (x,y)
ncm_apex* ncm_vessel::add_apex_at(double x, double y)
{
  ncm_apex* new_apex = new ncm_apex(*this, x, y);
  apices_.push_back(new_apex);
  set_modified();
  return new_apex;
}

//
//: Delete apex
void ncm_vessel::delete_apex(unsigned i)
{
  vcl_vector<ncm_apex*>::iterator apex = apices_.begin() + i;
  if (apex < apices_.end())
  {
    delete *apex;
    apices_.erase(apex);
    set_modified();
  }
}

//
//: Delete apex
void ncm_vessel::delete_apex(ncm_apex* const apex_address)
{
  vcl_vector<ncm_apex*>::iterator it = 
    vcl_find(apices_.begin(), apices_.end(), apex_address);

  if (it != apices_.end())
  {
    delete *it;
    apices_.erase(it);
    set_modified();
  }
}

//
//: Equality operator
bool ncm_vessel::operator==(const ncm_vessel& rhs)
{
  if (&rhs == this)
    return true;

  return false;
}

//
//: Inequality (requires equality to be defined)
bool ncm_vessel::operator!= (ncm_vessel const& rhs)
{ 
  return !(*this==rhs); 
}

//
//  Private methods
//

//
//: Delete all apices
void ncm_vessel::delete_all_apices()
{
  if (n_apices() > 0)
  {
    // delete apices associated with this vessel
    for (unsigned i = 0; i < apices_.size(); ++i)
      delete apices_[i];

    apices_.resize(0);
    set_modified();
  }
}
