#ifndef __ncm_user_h__
#define __ncm_user_h__

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

#include <QString>
#include <QStringList>

// forward declarations
class vsl_b_ostream;
class vsl_b_istream;
class vsl_ostream;
class vsl_istream;
class QSqlRecord;
class QSqlQuery;

class ncm_user
{
//  INTERFACE

public:
  //  No member variables here, please

  //: Constructors
  ncm_user(int id = 1);

  //  The 'Big Three'

  //: Destructor
  ~ncm_user();
  

  //: Assignment operator
  //ncm_user& operator=(ncm_user const& rhs); // assignment

  //: Name of the user
  QString name() const;
	QString host() const;
	int user_id() const;

	void set_name(QString name);
	void set_host(QString host);

  //  Other operators

  //bool operator==(ncm_user const& rhs); // equality
  //bool operator!=(ncm_user const& rhs); // inequality

	/*
  //  Functions for binary IO
  short version() const;
  vcl_string is_a() const;
  bool is_class(const vcl_string& class_name) const;
  //ncm_user* clone() const;

	//: Functions for sql IO
  //bool sql_write(QSqlQuery& query); - users created offline? 
  void sql_read(QSqlRecord& record);

  //: Functions for text IO
  void t_write(vcl_ostream& os) const;
  void t_read(vcl_istream& is);
  

  // Friends
  // These functions call protected member functions (e.g. print_summary()) and
  // therefore must be 'friend's.

  //: Binary file stream output operator for class reference
  friend void vsl_b_write(vsl_b_ostream& bfs, const ncm_user& b);

  //: Binary file stream input operator for class reference
  friend void vsl_b_read(vsl_b_istream& bfs, ncm_user& b);

  //: Stream output operator for class reference
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_user& b);

  //: Stream output operator for class pointer
  friend vcl_ostream& operator<<(vcl_ostream& os,const ncm_user* b);
	*/


//  IMPLEMENTATION

protected:
	/*
  // Members and functions visible to objects of this class and derived classes
  void b_write(vsl_b_ostream& os) const;
  void b_read(vsl_b_istream& is);
  void print_summary(vcl_ostream& os) const;
	*/
private:
  // Members and functions visible only to objects of this class
 
  //Name:
	QString name_;

	//Host
	QString host_;

	//ID:
	int user_id_;

	/*
	//List of fields in the database
	const static QStringList sql_fields_;
	const static QStringList sql_updateable_fields_;
	*/
};

#endif // __ncm_user_h__