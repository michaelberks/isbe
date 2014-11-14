#ifndef __ncm_image_session_h__
#define __ncm_image_session_h__

//:
// \file
// \brief A image_session of the NCM system
// \author Mike Berks
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_vector.h>
#include <vcl_string.h>
#include <vcl_iosfwd.h>
#include <QDateTime>
#include <QString>

#include "nailfold/qtools/ncm_image_sequence.h"

// forward declarations
class vsl_b_ostream;
class vsl_b_istream;
class vsl_ostream;
class vsl_istream;
class QSqlRecord;
class QSqlQuery;
class QDateTime;
class QStringList;

class ncm_image_session
{
//  INTERFACE

public:
  //  No member variables here, please

  //: Constructors
	ncm_image_session(int subject_id = -1, int user_id = -1);

  //: Destructor
  ~ncm_image_session();
  
	//Reset all the fields to empty
	void reset();

  //: Assignment operator
  //ncm_image_session& operator=(ncm_image_session const& rhs); // assignment

	int session_id() const;
	vcl_vector<ncm_image_sequence>* image_sequences();
	QString notes() const;

	//Note all other 'set' functions are private, limiting the ways the member variables can
	//be set to the constructor and SQL read
	void set_notes(QString notes);

	
	//:Date and time the session was started/finished
	QDateTime time_started() const;
	QDateTime time_finished() const;

	//:Start/end the session
	void start();
	void finalise();

	//:Add a new image sequence
	void add_image_sequence(const ncm_image_sequence sequence);

	//Generate string summarise sequences belonging to this session
	QString sequences_string() const;

  //  Other operators

  //bool operator==(ncm_image_session const& rhs); // equality
  //bool operator!=(ncm_image_session const& rhs); // inequality

  //  Functions for binary IO
  short version() const;
  vcl_string is_a() const;
  bool is_class(const vcl_string& class_name) const;
  //ncm_image_session* clone() const;

  //: Functions for text IO
  void t_write(vcl_ostream& os) const;
  void t_read(vcl_istream& is);

	//: Functions for sql IO
	//Note the write function can't be const, as after writing to SQL we need to get the
	//auto-incremented subject ID and update the subject_id member variable
  bool sql_write(QSqlQuery& query);
	bool sql_update(QSqlQuery& query, const QString& fieldname, const QVariant &value) const;
	void sql_read(QSqlRecord& record);
	void sql_read(QSqlRecord& record, QSqlQuery& query);
	bool sql_delete(QSqlQuery& query);
	bool sql_add_to_process_list(QSqlQuery& query, bool process_bad);
	

	//Display info for debugging
	void display_details() const;
  

  // Friends
  // These functions call protected member functions (e.g. print_summary()) and
  // therefore must be 'friend's.

  //: Binary file stream output operator for class reference
  friend void vsl_b_write(vsl_b_ostream& bfs, const ncm_image_session& b);

  //: Binary file stream input operator for class reference
  friend void vsl_b_read(vsl_b_istream& bfs, ncm_image_session& b);

  //: Stream output operator for class reference
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_image_session& b);

  //: Stream output operator for class pointer
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_image_session* b);



//  IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes
  void b_write(vsl_b_ostream& os) const;
  void b_read(vsl_b_istream& is);
  void print_summary(vcl_ostream& os) const;

private:
  // Members and functions visible only to objects of this class

	//These should only be used by SQL read
	void set_session_id(int id);
	void set_user_id(int id);
	void set_subject_id(int id);
	void set_time_started(QDateTime time);
	void set_time_finished(QDateTime time);
	
	//Session ID: unique identifier in SQL database
	int session_id_;

	//Unique ID of user who created the session
	int user_id_;

	//Unique ID of subject being imaged
	int subject_id_;

	//Start/finish time of session
	QDateTime time_started_;
	QDateTime time_finished_;

	//Sequences captured in this session
	vcl_vector<ncm_image_sequence> sequences_;

	//Any notes made by the user about the session
	QString notes_;

	//List of fields in the database
	const static QStringList sql_fields_;
	const static QStringList sql_updateable_fields_;
	
};

#endif // __ncm_image_session_h__