#ifndef __ncm_image_sequence_h__
#define __ncm_image_sequence_h__

//:
// \file
// \brief A image_sequence of the NCM system
// \author Mike Berks
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_vector.h>
#include <vcl_string.h>
#include <vcl_iosfwd.h>

#include "nailfold/qtools/ncm_video_frame.h"
#include "nailfold/ncm_hand.h"


#include <QString>
#include <QStringList>
#include <QDateTime>

// forward declarations
class vsl_b_ostream;
class vsl_b_istream;
class vsl_ostream;
class vsl_istream;
class QSqlQuery;
class QSqlRecord;
class QVariant;

class ncm_image_sequence
{
//  INTERFACE

public:
  //  No member variables here, please

  //: Constructors
	ncm_image_sequence(int session_id = -1);

  //  The 'Big Three'

  //: Destructor
  ~ncm_image_sequence();
  

  //: Assignment operator
  //ncm_image_sequence& operator=(ncm_image_sequence const& rhs); // assignment

	void start();
	void finalise();

	//Set functions
	void set_hand(ncm_hand::Hand hand);
	void set_digit(ncm_hand::Digit digit);
	void set_images_root_dir(QString dirpath);
	void set_images_sub_dir(QString dirpath);
	void set_images_full_dir(QString dirpath);
	void set_images_datafile(QString filename);
	void set_sequence_name(QString name);
	void set_num_frames(int n_frames);
	void set_notes(QString notes);
	void set_acceptable_quality(bool yesno);
	void set_camera_type(QString camera);
	void set_camera_id(QString id);
	void set_frame_rate(double fps);
	void set_exposure_time(double t);
	void set_gain(double g);
	void set_apt_x_id(long id);
	void set_apt_y_id(long id);
	void set_apt_z_id(long id);
	void set_x_reversed(bool rev);
	void set_y_reversed(bool rev);
	void set_z_reversed(bool rev);
	void set_preview_mosaic(QString filename);
	void set_full_mosaic(QString filename);
	void set_preview_count_map(QString filename);
	void set_full_count_map(QString filename);

	//Get functions
	ncm_hand::Hand hand() const;
	ncm_hand::Digit digit() const;
	QString images_root_dir() const;
	QString images_sub_dir() const;
	QString images_full_dir() const;
	QString images_datafile() const;
	QString sequence_name() const;
	int num_frames() const;
	QString notes() const;
	bool acceptable_quality() const;
	QString camera_type() const;
	QString camera_id() const;
	double frame_rate() const;
	double exposure_time() const;
	double gain() const;
	int sequence_id() const; //Note this has no public set
	long apt_x_id() const;
	long apt_y_id() const;
	long apt_z_id() const;
	bool x_reversed() const;
	bool y_reversed() const;
	bool z_reversed() const;
	QString preview_mosaic() const;
	QString full_mosaic() const;
	QString preview_count_map() const;
	QString full_count_map() const;

	ncm_video_frame_header frame_header(int i) const;
	
	//:Date and time the session was started/finished
	QDateTime time_started() const;
	QDateTime time_finished() const;

	//:Add details of a video frame
	void add_frame(ncm_video_frame_header frame);

  //bool operator==(ncm_image_sequence const& rhs); // equality
  //bool operator!=(ncm_image_sequence const& rhs); // inequality

  //  Functions for binary IO
  short version() const;
  vcl_string is_a() const;
  bool is_class(const vcl_string& class_name) const;
  //ncm_image_sequence* clone() const;

	//: Functions for sql IO
	//Note the write function can't be const, as after writing to SQL we need to get the
	//auto-incremented subject ID and update the subject_id member variable
  bool sql_write(QSqlQuery& query);
	bool sql_update(QSqlQuery& query, const QString& fieldname, const QVariant &value) const;
	void sql_read(QSqlRecord& record);

	//:Try and sync the sequence name and the images directory
	bool sync_name_and_dir(QSqlQuery& query);

  //: Functions for text IO
  void t_write(vcl_ostream& os) const;
  void t_read(vcl_istream& is);
  void t_read(const vcl_string& filename);

	//Display function for use during development
	void ncm_image_sequence::display_details() const;

  // Friends
  // These functions call protected member functions (e.g. print_summary()) and
  // therefore must be 'friend's.

  //: Binary file stream output operator for class reference
  friend void vsl_b_write(vsl_b_ostream& bfs, const ncm_image_sequence& b);

  //: Binary file stream input operator for class reference
  friend void vsl_b_read(vsl_b_istream& bfs, ncm_image_sequence& b);

  //: Stream output operator for class reference
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_image_sequence& b);

  //: Stream output operator for class pointer
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_image_sequence* b);



//  IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes
  void b_write(vsl_b_ostream& os) const;
  void b_read(vsl_b_istream& is);
  void print_summary(vcl_ostream& os) const;

private:
  // Members and functions visible only to objects of this class
	void clear();

	//:Can only be set internally (by SQL read)
	void set_time_started(QDateTime time);
	void set_time_finished(QDateTime time);
	void set_sequence_id(int id);

	//Helper functions for reading and writing (these may move to a new class?)
	void read_camera_properties(vcl_istream& is);
	void read_motor_properties(vcl_istream& is);
	void read_frame(vcl_istream& is, int frame_num);

	void write_camera_properties(vcl_ostream& os) const;
	void write_motor_properties(vcl_ostream& os) const;
	void write_frame(vcl_ostream& os, ncm_video_frame_header fh) const;

	//Unique identifier for this session (created by auto-increment in SQL database)
	int sequence_id_;

	//ID of the session to which this sequence belongs
	int session_id_;
	
	//Filepaths to directory containing the session images and datafile
	//Easier to store the final subdir separately so this can be separtely updated by the GUI
	QString images_root_dir_;
	QString images_sub_dir_;

	//User given name for the sequence - this does not need to be the same as
	//images sub dir, however the function sync_name_and_dir() will try and do this
	//allowing the sequence name to be set after the sequence is saved
	QString sequence_name_;

	//Filename of the text file storing sequence information
	QString images_datafile_;

	//Hand and digit imaged
	ncm_hand::Hand hand_;
	ncm_hand::Digit digit_;

	//Time and date started and finished
	QDateTime time_started_;
	QDateTime time_finished_;

	//Number of frames in the sequence
	int num_frames_;

	//List of frame headers for frames captured in this sequence
	vcl_vector<ncm_video_frame_header> frame_headers_;

	//Free text for user to make comments on this session
	QString notes_;

	//Flag if the sequence is of acceptable quality for processing
	bool acceptable_quality_;

	//Properties of the camera for this sequence
	QString camera_type_;
	QString camera_id_;
	double frame_rate_;
	double exposure_time_; //negative value -> auto exposure
	double gain_; //negative value -> auto gain

	//Properties of motors for this sequence
	long apt_x_id_;
	long apt_y_id_;
	long apt_z_id_;
	bool x_reversed_;
	bool y_reversed_;
	bool z_reversed_;

	//Filenames of the mosaics, empty if not yet set
	QString preview_mosaic_;
	QString full_mosaic_;
	QString preview_count_map_;
	QString full_count_map_;

	//List of fields in the database
	const static QStringList sql_fields_;
	const static QStringList sql_updateable_fields_;

	//Format for date strings
	const static QString DATE_FORMAT;
};

#endif // __ncm_image_sequence_h__