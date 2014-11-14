#ifndef __ncm_qseries_grade_h__
#define __ncm_qseries_grade_h__

//:
// \file
// \brief Structure for holding the grades of an nailfold capillaroscopy series comparison
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

class ncm_qseries_grade
{
//  INTERFACE

public:
  //  No member variables here, please

	enum ImageGrade { Normal = 0,
                    Abnormal = 1,
										Ungradeable = 2};

	enum ProgressionLevel {
										MajorTop = 2,
										MinorTop = 1,
										None = 0,
										MinorBottom = -1,
										MajorBottom = -2};

		enum ProgressReason { NoReasons = -1,
                    Reason0 = 0,
                    Reason1 = 1,
                    Reason2 = 2,
                    Reason3 = 3,
                    Reason4 = 4,
                    Reason5 = 5,
                    Reason6 = 6,
										Reason7 = 7,
										Reason8 = 8,
										Reason9 = 9,
										Reason10 = 10,
										Reason11 = 11,
										Reason12 = 12,
										NumReasons = Reason12};

  //: Default constructor
  ncm_qseries_grade(/*T param1 = def1, T param2 = def2, ...*/);

  //: Destructor
  ~ncm_qseries_grade();

  //  Copy constructor and operator= are private (see below)
  //ncm_qseries_grade(ncm_qseries_grade const& rhs);
  //ncm_qseries_grade& operator=(ncm_qseries_grade const& rhs);

  //: Clear all qseries_grade
  void clear();

  //: Define location for saving the data
  void set_filename(const vcl_string& filename);
  vcl_string filename() const;

  //: Generate a tagged copy of the filename in the format:
  //  YYMMDD-hhmmss_<observer_name>#<filename>
  vcl_string tagged_filename() const;

  //: Save to defined filename (return false if filename not defined)
  bool save();

  //: Save to temporary backup file (return false if not done)
  bool snapshot() const;

  //: Save to temporary backup filename
  vcl_string snapshot_filename() const;

  //: Set the observer's name
  void set_observer_name(const vcl_string& observer_name);

  //: Enable/disable taking snapshots
  void set_snapshots(const bool snapshot_on_modify = true);

  //: Reset clock timer (used for storing grading time)
  //  Also starts timing from this moment
  void start_timing();

  //: Store length of time taken to annotate image
  void stop_timing();

  //: Return date and time as a readable string
  vcl_string date_string() const;

  //: Grade has been modified
  void set_modified(const bool is_modified = true) const;
	void set_normal(const bool is_modified = true);
	void set_reason(const bool selected, const ProgressReason reason);
	void set_progression(const ProgressionLevel level);

  //: Whether the grade has been modified
  bool is_modified() const; 
  bool is_graded() const;
	bool is_normal() const;
	bool is_reasons();
	vcl_vector< bool > progress_reasons();
	ProgressionLevel progression_level();
	bool progress_reasons(ProgressReason reason);

  ////: Other operators
  //bool operator==(ncm_qseries_grade const& rhs); // equality
  //bool operator!=(ncm_qseries_grade const& rhs); // inequality

  //: Functions for binary IO
  short version() const;
  vcl_string is_a() const;
  bool is_class(const vcl_string& class_name) const;
  //ncm_qseries_grade* clone() const;

  //  Functions for text IO
  void t_write(vcl_ostream& tfs) const;
  void t_read(vcl_istream& tfs);
  void t_read(const vcl_string& filename);
  

  // Friends
  // These functions call protected member functions (e.g. print_summary()) and
  // therefore must be 'friend's.

  //: Binary file stream output operator for class reference
  friend void vsl_b_write(vsl_b_ostream& bfs, const ncm_qseries_grade& b);

  //: Binary file stream input operator for class reference
  friend void vsl_b_read(vsl_b_istream& bfs, ncm_qseries_grade& b);

  //: Stream output operator for class reference
  friend vcl_ostream& operator<<(vcl_ostream& os, const ncm_qseries_grade& b);

  //: Stream output operator for class pointer
  friend vcl_ostream& operator<<(vcl_ostream& os, const ncm_qseries_grade* b);



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
  ncm_qseries_grade(ncm_qseries_grade const& rhs);
  ncm_qseries_grade& operator=(ncm_qseries_grade const& rhs);

  //
  //  Variables
  //

  //: Filename of the grade file
  vcl_string filename_;

  //: Whether to save a snapshot every time set_modified() is called
  bool snapshot_on_modify_;

  //: Observer (user)name
  vcl_string observer_name_;

  //: Date and time of grading
  vcl_time_t datetime_;

  //: Time spent (in seconds) annotating image
  double time_taken_;

  //: Whether the grade has been modified
  //  Since this does not form part of the grade itself, we'll make it
  //  mutable so that it can be modified even for const objects
  mutable bool is_modified_;	
	bool is_normal_;
	bool is_graded_;
	vcl_vector< bool > progress_reasons_;
	ProgressionLevel progression_level_;

};

//: Add to loader for polymorphic IO
//void vsl_add_to_binary_loader(const ncm_annotation& b);

#endif // __ncm_qseries_grade_h__