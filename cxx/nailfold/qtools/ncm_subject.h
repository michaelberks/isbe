#ifndef __ncm_subject_h__
#define __ncm_subject_h__

//:
// \file
// \brief A user of the NCM system
// \author Mike Berks
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_vector.h>
#include <vcl_string.h>
#include <vcl_iosfwd.h>
#include "nailfold/ncm_hand.h"
#include "ncm_image_session.h"
#include <QString>

// forward declarations
class vsl_b_ostream;
class vsl_b_istream;
class vsl_ostream;
class vsl_istream;
class QSqlRecord;
class QSqlQuery;
class QDateTime;

class ncm_subject
{
//  INTERFACE

public:
  //  No member variables here, please

  //: Constructors
  ncm_subject();

  //  The 'Big Three'

  //: Destructor
  ~ncm_subject();
  

  //: Assignment operator
  //ncm_subject& operator=(ncm_subject const& rhs); // assignment

	//Reset all the fields to empty
	void reset();

  //: Get basic attirbutes of the subject
  QString name() const;
	int ncm_id() const;
	int study_id() const;
	QString name_in_study() const;
	ncm_hand::Hand dominant_hand() const;
	QString notes() const;
	QDateTime first_imaged() const;
	QDateTime last_imaged() const;

	//Return string compunding information on previous sessions into displayable form
	QString previous_sessions_string() const;

	//: Return pointer to live session object
	ncm_image_session* live_session();

	//: Return pointer to session selected by idx
	//ncm_image_session* select_session(int idx);

	//: Return a copy of the object
	ncm_image_session select_session(int idx);

	//: Return the session index of session with given id
	int find_session(int id);

	//: Set basic attirbutes of the subject
	void set_name(QString name);
	void set_ncm_id(int id);
	void set_study_id(int id);
	void set_name_in_study(QString name);
	void set_dominant_hand(ncm_hand::Hand hand);
	void set_notes(QString notes);

	void set_live_session(ncm_image_session session);
	void add_live_to_previous();

  //  Other operators

  //bool operator==(ncm_subject const& rhs); // equality
  //bool operator!=(ncm_subject const& rhs); // inequality

  //  Functions for binary IO
  short version() const;
  vcl_string is_a() const;
  bool is_class(const vcl_string& class_name) const;
  //ncm_subject* clone() const;

  //: Functions for text IO
  void t_write(vcl_ostream& os) const;
  void t_read(vcl_istream& is);

	//: Functions for sql IO
	//Note the write function can't be const, as after writing to SQL we need to get the
	//auto-incremented subject ID and update the subject_id member variable
  bool sql_write(QSqlQuery& query); 
  void sql_read(QSqlRecord& record);
	void sql_read(QSqlRecord& record, QSqlQuery& query);
  bool sql_update(QSqlQuery& query, const QString& fieldname, const QVariant &value) const;

  // Friends
  // These functions call protected member functions (e.g. print_summary()) and
  // therefore must be 'friend's.

  //: Binary file stream output operator for class reference
  friend void vsl_b_write(vsl_b_ostream& bfs, const ncm_subject& b);

  //: Binary file stream input operator for class reference
  friend void vsl_b_read(vsl_b_istream& bfs, ncm_subject& b);

  //: Stream output operator for class reference
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_subject& b);

  //: Stream output operator for class pointer
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_subject* b);



//  IMPLEMENTATION

protected:
  // Members and functions visible to objects of this class and derived classes
  void b_write(vsl_b_ostream& os) const;
  void b_read(vsl_b_istream& is);
  void print_summary(vcl_ostream& os) const;

private:
  // Members and functions visible only to objects of this class
 
  //Name:
	QString name_;
	int ncm_id_;
	int study_id_;
	QString name_in_study_;
	ncm_hand::Hand dominant_hand_;

	vcl_vector<ncm_image_session> image_sessions_;
	ncm_image_session live_session_;
	QString notes_;

	//List of fields in the database
	const static QStringList sql_fields_;
	const static QStringList sql_updateable_fields_;
};

#endif // __ncm_user_h__