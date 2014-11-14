#include "nailfold/ncm_annotation.h"

#include <vcl_ctime.h>
#include <vcl_algorithm.h>
#include <vcl_iostream.h>
#include <vcl_sstream.h>
#include <vcl_iomanip.h>

#include <vul/vul_file.h>

#include <mbl/mbl_index_sort.h>

#include <nailfold/ncm_encrypt.h>
#include <nailfold/ncm_vessel.h>
#include <nailfold/ncm_vessel_properties.h>
#include <nailfold/ncm_haemorrhage.h>

//: Define static class members here
// type ncm_annotation::member_ = 0;

//: Define static const class members here
// type ncm_annotation::const_member_ = 0;

//
//  Static member variables
//

bool ncm_annotation::vessels_are_sorted_ = true;

//
//  Static methods
//

//
//: Choose whether to sort vessels by x-coordinate
void ncm_annotation::set_vessel_sorting(bool vessels_are_sorted /* = true */)
{
  vessels_are_sorted_ = vessels_are_sorted;
}
bool ncm_annotation::vessels_are_sorted()
{
  return vessels_are_sorted_;
}

//
// Public methods
//

//
//: Default constructor
ncm_annotation::ncm_annotation()
: snapshot_on_modify_(true),
  image_grade_(this),
	version_(2)
{
  clear();
}

//
//: Destructor
ncm_annotation::~ncm_annotation()
{
  set_snapshots(false);
  clear();
}

//
//: Clear all data
void ncm_annotation::clear()
{
  delete_all_vessels();
  delete_all_haemorrhages();

  image_grade().set_value(ncm_image_grade::GradeUndefined);
  observer_name_ = "-";
  datetime_ = static_cast<vcl_time_t>(0);
  time_taken_ = 0.0;

  is_modified_ = false;
}

//
//: Image filename
void ncm_annotation::set_image_filename(const vcl_string& image_filename)
{
  image_filename_ = image_filename;
}

//
//: Markup filename
void ncm_annotation::set_filename(const vcl_string& filename)
{
  filename_ = filename;
}
vcl_string ncm_annotation::filename() const
{
  return filename_;
}

//
//: Generate a tagged copy of the filename in the format:
//  YYMMDD-hhmmss_<observer_name>#<filename>
vcl_string ncm_annotation::tagged_filename() const
{
  // Create timestamp
  time_t rawtime = time(NULL);
  struct tm* ptm = gmtime(&rawtime);
  vcl_stringstream ss;
  ss << ptm->tm_year+1900 
     << vcl_setw(2) << vcl_setfill('0') << ptm->tm_mon+1 // months start at 0
     << vcl_setw(2) << vcl_setfill('0') << ptm->tm_mday
     << "-"
     << vcl_setw(2) << vcl_setfill('0') << ptm->tm_hour 
     << vcl_setw(2) << vcl_setfill('0') << ptm->tm_min 
     << vcl_setw(2) << vcl_setfill('0') << ptm->tm_sec;
  vcl_string timestamp = ss.str();

  return timestamp + "_" + 
         observer_name_ + "#" + 
         vul_file::strip_directory(filename_);
}

bool ncm_annotation::save()
{
  if (filename_.empty())
    return false;

  // else
  stop_timing();
  vcl_ofstream ofs(filename_.c_str());
  t_write(ofs);
  ofs.close();

  // No longer modified since last save
  is_modified_ = false;

  return true;
}

void ncm_annotation::set_snapshots(const bool snapshot_on_modify /* = true */)
{
  snapshot_on_modify_ = snapshot_on_modify;
}
bool ncm_annotation::snapshot() const
{
  if (image_filename_.empty())
    return false;

  // else
  vcl_ofstream ofs(snapshot_filename().c_str());
  t_write(ofs);
  ofs.close();
  return true;
}

vcl_string ncm_annotation::snapshot_filename() const
{
  return vul_file::strip_extension(image_filename_) + "_markup.txt~";
}

//
//:
ncm_image_grade& ncm_annotation::image_grade()
{
  return image_grade_;
}
ncm_image_grade const& ncm_annotation::image_grade() const
{
  return image_grade_;
}

//
//: Set the observer's name
void ncm_annotation::set_observer_name(const vcl_string& observer_name)
{
  if (observer_name_ != observer_name)
  {
    observer_name_ = observer_name;
    set_modified();
  }
}

//
//: Reset clock timer (used for storing annotation time)
void ncm_annotation::start_timing()
{
  vcl_time(&datetime_);

  // Don't reset time_taken_ as we will add the time spent annotating to the
  // stored time (if any) to reflect the total time spent marking up the image

  // Don't use this as a cue to set_modified() - instead, we'll wait until the 
  // user makes an observable change to the markup.
}

//
//: Store length of time taken to annotate image
void ncm_annotation::stop_timing()
{
  // Update total time spent so far marking up this image
  time_taken_ += vcl_difftime(vcl_time(NULL), datetime_);
  //set_modified();
}

//
//: Return date and time as a string
vcl_string ncm_annotation::date_string() const
{
  char* date_string = vcl_asctime(vcl_gmtime(&datetime_));

  if (date_string != NULL)
    return date_string;
  else
    return "-";
}

//
//: Annotation has been modified (or force reset)
void ncm_annotation::set_modified(const bool is_modified /* = true */) const
{
  // is_modified_ is mutable and therefore allowed to change
  is_modified_ = is_modified;

  if (is_modified && snapshot_on_modify_)
    snapshot();
}

//
//: Whether the annotation has been modified
bool ncm_annotation::is_modified() const
{
  return is_modified_;
}

//
//  Vessel methods
//

//
//: Create new vessel, add to vector and return pointer
ncm_vessel* ncm_annotation::create_vessel_at(
    double x, double y,
    ncm_vessel_properties* properties /* = NULL */)
{
  ncm_vessel* new_vessel = ncm_vessel_handler::new_vessel(*this, x, y);

  if (properties != NULL)
    new_vessel->properties() = *properties;

  vessels_.push_back(new_vessel);

  update_indices();

  set_modified();
  return new_vessel;
}

//
//: Delete just the selected vessel (if it is part of this annotation)
//  Returns true if vessel found and deleted successfully, false otherwise
bool ncm_annotation::delete_vessel(ncm_vessel* vessel)
{
  vcl_vector<ncm_vessel*>::iterator it;
  it = vcl_find(vessels_.begin(), vessels_.end(), vessel);

  if (it != vessels_.end())
  {
    ncm_vessel_handler::delete_vessel(*it);
    vessels_.erase(it);
    update_indices();
    set_modified();
    return true;
  }
  else
    return false;
}

//
//: Delete all vessels that we've created
void ncm_annotation::delete_all_vessels()
{
  if (vessels_.empty())
    return;

  vcl_vector<ncm_vessel*>::iterator it;
  for (it = vessels_.begin(); it != vessels_.end(); ++it)
  {
    ncm_vessel_handler::delete_vessel(*it);
  }
  vessels_.clear();
  set_modified();

  index_of_vessel_.clear();
}

//
//: Number of vessels
unsigned ncm_annotation::n_vessels() const
{
  return vessels_.size();
}
unsigned ncm_annotation::n_distal() const
{
  unsigned n_distal = 0;
  for (unsigned i = 0; i < vessels_.size(); ++i)
  {
    if (vessels_[i]->properties().is_distal())
      ++n_distal;
  }
  return n_distal;
}
//
//: Get pointer to vessel number i
ncm_vessel* ncm_annotation::vessel(unsigned i)
{
  if (i < n_vessels())
    return vessels_[ index_of_vessel_[i] ];
  else
    return NULL;
}

const ncm_vessel* ncm_annotation::vessel(unsigned i) const
{
  if (i < n_vessels())
    return vessels_[ index_of_vessel_[i] ];
  else
    return NULL;
}

int ncm_annotation::index_of_vessel(ncm_vessel const* vessel) const
{
  vcl_vector<ncm_vessel*>::const_iterator it;
  it = vcl_find(vessels_.begin(), vessels_.end(), vessel);

  if (it == vessels_.end())
    return -1; // vessel not found
  else
  {
    int unsorted_index = it - vessels_.begin(); // offset

    // find element in index_of_vessel_ for this vessel
    vcl_vector<unsigned>::const_iterator index_it;
    index_it = vcl_find(index_of_vessel_.begin(), index_of_vessel_.end(),
                        unsorted_index);

    if (index_it == index_of_vessel_.end())
      return -1; // vessel not found (should never happen).
    else
    {
      int sorted_index = index_it - index_of_vessel_.begin();
      return sorted_index;
    }
  }
}

//
//: Find the nearest *distal* vessel to the given coordinates, pos.
//  Store the distance in minimum_distance if a double pointer is given.
ncm_vessel const* ncm_annotation::vessel_nearest_to(
    const vgl_point_2d<double>& pos,
    double* minimum_distance /* = NULL */) const
{
  ncm_vessel const* nearest_vessel = NULL;
  double min_distance = 1e9;

  for (unsigned v = 0; v < n_distal(); ++v)
  {
    ncm_vessel const* test_vessel = vessel(v);

    const vgl_vector_2d<double> offset = pos - test_vessel->anchor();
    const double distance = offset.length();

    if (distance < min_distance)
    {
      nearest_vessel = test_vessel;
      min_distance = distance;
    }
  }

  // Export distance value also
  if (minimum_distance != NULL)
    *minimum_distance = min_distance;

  return nearest_vessel;
}

//
//  Haemorrhage methods
//
ncm_haemorrhage* ncm_annotation::create_haemorrhage_at(double x, double y)
{
  ncm_haemorrhage* new_haemorrhage = new ncm_haemorrhage(*this, x, y);
  haemorrhages_.push_back(new_haemorrhage);
  set_modified();
  return new_haemorrhage;
}

bool ncm_annotation::delete_haemorrhage(ncm_haemorrhage* haemorrhage)
{
  vcl_vector<ncm_haemorrhage*>::iterator it;
  it = vcl_find(haemorrhages_.begin(), haemorrhages_.end(), haemorrhage);

  if (it != haemorrhages_.end())
  {
    delete *it;
    haemorrhages_.erase(it);
    set_modified();
    return true;
  }

  return false;
}

void ncm_annotation::delete_all_haemorrhages()
{
  if (haemorrhages_.empty())
    return;

  vcl_vector<ncm_haemorrhage*>::iterator it;
  for (it = haemorrhages_.begin(); it != haemorrhages_.end(); ++it)
  {
    delete *it;
  }
  haemorrhages_.clear();
  set_modified();
}

unsigned ncm_annotation::n_haemorrhages() const
{
  return haemorrhages_.size();
}

ncm_haemorrhage* ncm_annotation::haemorrhage(unsigned i)  
{
  if (i < n_haemorrhages())
    return haemorrhages_[i];
  else
    return NULL;
}

ncm_haemorrhage const* ncm_annotation::haemorrhage(unsigned i) const
{
  if (i < n_haemorrhages())
    return haemorrhages_[i];
  else
    return NULL;
}

 
//: Return union of missing data flags to indicate missing annotations.
unsigned ncm_annotation::missing_data() const
{
  unsigned return_value = MissingData_None;

  if (image_grade_.value() == ncm_image_grade::GradeUndefined)
    return_value = return_value | MissingData_Grade;

  if (n_vessels() == 0)
    return_value = return_value | MissingData_Vessels;
  
  for (unsigned i = 0; i < n_distal(); ++i)
  {
    if (vessel(i)->properties().is_size_undefined())
      return_value = return_value | MissingData_Sizes;

    if (vessel(i)->properties().is_shape_undefined())
      return_value = return_value | MissingData_Shapes;

    if (vessel(i)->n_apices() == 0)
      return_value = return_value | MissingData_Apices;

    if (vessel(i)->n_points() == 0)
      return_value = return_value | MissingData_Paths;
  }

  if (n_haemorrhages() == 0)
    return_value = return_value | MissingData_Haemorrhages;

  return return_value;
}

//
//  Private methods
//

void ncm_annotation::update_indices()
{
  vcl_vector<unsigned> index_of_vessel;

  if (vessels_are_sorted_)
  {
    // Create vector of x-coordinates first
    vcl_vector<double> x_coords(n_vessels());
    for (unsigned i = 0; i < n_vessels(); ++i)
    {
      // must use vessels_[i] rather than vessel(i) here
      x_coords[i] = vessels_[i]->anchor().x();
    }

    mbl_index_sort(x_coords, index_of_vessel);
  }
  else
  {
    // This effectively applies, in the same framework, no sorting.

    index_of_vessel.resize(n_vessels());
    for (unsigned i = 0; i < n_vessels(); ++i)
    {
      // must use vessels_[i] rather than vessel(i) here
      index_of_vessel[i] = i;
    }
  }

  // Now partition the index vector such that the M distal vessels appear first,
  // and the N nondistal second (e.g. [d1, d2, ..., dM, nd1, nd2, ..., ndN]
  index_of_vessel_.resize(0);
  for (unsigned i = 0; i < index_of_vessel.size(); ++i)
  {
    if (vessels_[index_of_vessel[i]]->properties().is_distal())
      index_of_vessel_.push_back(index_of_vessel[i]);
  }
  for (unsigned i = 0; i < index_of_vessel.size(); ++i)
  {
    if (!(vessels_[index_of_vessel[i]]->properties().is_distal()))
      index_of_vessel_.push_back(index_of_vessel[i]);
  }
}
