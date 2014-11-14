#ifndef __ncm_annotation_h__
#define __ncm_annotation_h__

//:
// \file
// \brief Structure for holding the annotation of an nailfold capillaroscopy 
//        image
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_vector.h>
#include <vcl_string.h>
#include <vcl_ctime.h>

#include <vgl/vgl_point_2d.h>

//#include <nailfold/ncm_vessel.h>
#include <nailfold/ncm_image_grade.h>

class vsl_b_ostream;
class vsl_b_istream;

class ncm_vessel;
class ncm_vessel_properties;
class ncm_haemorrhage;

class ncm_annotation
{
//  INTERFACE

public:
  //  No member variables here, please

  //: Default constructor
  ncm_annotation(/*T param1 = def1, T param2 = def2, ...*/);

  //: Destructor
  ~ncm_annotation();

  //  Copy constructor and operator= are private (see below)
  //ncm_annotation(ncm_annotation const& rhs);
  //ncm_annotation& operator=(ncm_annotation const& rhs);

  //: Clear all annotation
  void clear();

  //: Define location of image (used for temporary autosave file)
  void set_image_filename(const vcl_string& image_filename);

  //: Define location for saving the data
  void set_filename(const vcl_string& filename);
  vcl_string filename() const;

  //: Generate a tagged copy of the filename in the format:
  //  YYMMDD-hhmmss_<observer_name>#<filename>
  vcl_string tagged_filename() const;

  //: Transmit the saved file to a remote server
  bool sendToServer() const;

  //: Save to defined filename (return false if filename not defined)
  bool save();

  //: Save to temporary backup file (return false if not done)
  bool snapshot() const;

  //: Save to temporary backup filename
  vcl_string snapshot_filename() const;

  //: Define the image class
  ncm_image_grade& image_grade();
  ncm_image_grade const& image_grade() const;

  //: Set the observer's name
  void set_observer_name(const vcl_string& observer_name);

  //: Enable/disable taking snapshots
  void set_snapshots(const bool snapshot_on_modify = true);

  //: Reset clock timer (used for storing annotation time)
  //  Also starts timing from this moment
  void start_timing();

  //: Store length of time taken to annotate image
  void stop_timing();

  //: Return date and time as a readable string
  vcl_string date_string() const;

  //: Annotation has been modified
  void set_modified(const bool is_modified = true) const;

  //: Whether the annotation has been modified
  bool is_modified() const;

  //: Choose whether to sort vessels
  static void set_vessel_sorting(bool vessels_are_sorted = true);
  static bool vessels_are_sorted();


  //  Vessel methods

  //: Create new vessel, add to vector and return pointer
  ncm_vessel* create_vessel_at(double x, double y, 
                               ncm_vessel_properties* properties = NULL);

  //: Delete just the selected vessel (if it is part of this annotation)
  //  Returns true if vessel found and deleted successfully
  bool delete_vessel(ncm_vessel* vessel);

  //: Clear just the marked vessels
  void delete_all_vessels();

  //: Number of vessels
  unsigned n_vessels() const;
  unsigned n_distal() const;
  unsigned n_nondistal() const;

  //: Get pointers to i'th vessel
  //  Returns NULL if i is out of range
  ncm_vessel* vessel(unsigned i);
  ncm_vessel const* vessel(unsigned i) const;

  int index_of_vessel(ncm_vessel const* vessel) const;

  //: Vessel closest to a given image location
  ncm_vessel const* vessel_nearest_to(const vgl_point_2d<double>& pos, 
                                      double* minimum_distance = NULL) const;

  //  Haemorrhage methods (equivalent to vessel methods)

  ncm_haemorrhage* create_haemorrhage_at(double x, double y);
  bool delete_haemorrhage(ncm_haemorrhage* haemorrhage);
  void delete_all_haemorrhages();
  unsigned n_haemorrhages() const;
  ncm_haemorrhage* haemorrhage(unsigned i);
  ncm_haemorrhage const* haemorrhage(unsigned i) const;

  //: Missing data flags
  enum missing_data_flag {
    MissingData_None = 0x0000,
    MissingData_Grade = 0x0001,
    MissingData_Vessels = 0x0002,
    MissingData_Sizes = 0x0004,
    MissingData_Shapes = 0x0008,
    MissingData_Apices = 0x0010,
    MissingData_Paths = 0x0020,
    MissingData_Haemorrhages = 0x0040,
    
    MissingData_First = MissingData_None,
    MissingData_Last = MissingData_Haemorrhages
  }; 

  //: Return union of missing data flags to indicate missing annotations.
  unsigned missing_data() const;

  ////: Other operators
  //bool operator==(ncm_annotation const& rhs); // equality
  //bool operator!=(ncm_annotation const& rhs); // inequality

  //: Functions for binary IO
  short version() const;
	short latest_manual_version() const;
	short latest_auto_version() const;
  vcl_string is_a() const;
  bool is_class(const vcl_string& class_name) const;
  //ncm_annotation* clone() const;

  //  Functions for text IO
  void t_write(vcl_ostream& tfs) const;
  void t_read(vcl_istream& tfs);
  void t_read(const vcl_string& filename);
  

  // Friends
  // These functions call protected member functions (e.g. print_summary()) and
  // therefore must be 'friend's.

  //: Binary file stream output operator for class reference
  friend void vsl_b_write(vsl_b_ostream& bfs, const ncm_annotation& b);

  //: Binary file stream input operator for class reference
  friend void vsl_b_read(vsl_b_istream& bfs, ncm_annotation& b);

  //: Stream output operator for class reference
  friend vcl_ostream& operator<<(vcl_ostream& os, const ncm_annotation& b);

  //: Stream output operator for class pointer
  friend vcl_ostream& operator<<(vcl_ostream& os, const ncm_annotation* b);



  // IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes

  //: If we provide an operator<<, there is no need for a public print_summary()
  //  so make it protected (and operator<< a friend).
  void print_summary(vcl_ostream& os) const;

  //: Similarly for vsl_b_write() and vsl_b_read() - no need for b_write() and
  //  b_read() to be public
  void b_write(vsl_b_ostream& bfs) const;
  void b_read(vsl_b_istream& bfs);

private:
  //  Members and functions visible only to objects of this class

  //: I can't think of a time when you would want to copy or assign to an
  //  instance of this class so make copy constructor and operator= private
  ncm_annotation(ncm_annotation const& rhs);
  ncm_annotation& operator=(ncm_annotation const& rhs);

  //: Update list of sorted indices of vessels
  void update_indices();

  //
  //  Variables
  //

	//Version of the annotation
	short version_;

  //: Filename of the input image
  vcl_string image_filename_;

  //: Filename of the annotation
  vcl_string filename_;

  //: Whether to save a snapshot every time set_modified() is called
  bool snapshot_on_modify_;

  //: Classification of the image by the observer
  ncm_image_grade image_grade_;

  //: Vessels
  vcl_vector<ncm_vessel*> vessels_;

  //: Haemorrhages
  vcl_vector<ncm_haemorrhage*> haemorrhages_;

  //: Observer (user)name
  vcl_string observer_name_;

  //: Date and time of annotation
  vcl_time_t datetime_;

  //: Time spent (in seconds) annotating image
  double time_taken_;

  //: Whether the annotation has been modified
  //  Since this does not form part of the annotation itself, we'll make it
  //  mutable so that it can be modified even for const annotations
  mutable bool is_modified_;

  //: Statistics over zoom level?
  //: Statistics over centre of attention?

  //: List of indices of vessels, sorted by their x-coordinate
  vcl_vector<unsigned> index_of_vessel_;

  //: Whether we're using sorting by x-coordinate
  //  Since it is not specific to a particular annotation, it is static
  static bool vessels_are_sorted_;
};

//: Add to loader for polymorphic IO
//void vsl_add_to_binary_loader(const ncm_annotation& b);

#endif // __ncm_annotation_h__