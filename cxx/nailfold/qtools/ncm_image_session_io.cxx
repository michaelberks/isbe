//:
// \file
// \brief Functions for Input/Output of ncm_image_session:
//          version()
//          is_a()
//          is_class()
//          vsl_b_write()
//          vsl_b_read()
//          t_write()
//          t_read()
//          operator<<
//          b_write()
//          b_read()
//          print_summary()
//          add_to_binary_loader()
// \author Phil Tresadern

#include "ncm_image_session.h"

#include <vsl/vsl_binary_io.h>
#include <vsl/vsl_vector_io.h>
#include <vsl/vsl_indent.h>

#include <vul/vul_string.h>

#include <mbl/mbl_parse_block.h>
#include <mbl/mbl_read_props.h>
#include <mbl/mbl_parse_tuple.h>
#include <mbl/mbl_parse_sequence.h>

#include <QStringList>
#include <QVariant>
#include <QtSql/QSqlQuery>
#include <QtSql/QSqlError>
#include <QtSql/QSqlRecord>

//Initialize the list of fieldmames in the SQL database
const QStringList ncm_image_session::sql_fields_ = QStringList() << 
	"subject_id" <<
	"user_id" <<
	"time_started" <<
	"time_finished" <<
	"session_comments";

const QStringList ncm_image_session::sql_updateable_fields_ = QStringList() << 
	"time_finished" <<
	"session_comments";

//
// Public functions

//
//: Version number of this class's binary IO
short ncm_image_session::version() const
{ 
  return 1; 
}

//
//: The class name
vcl_string ncm_image_session::is_a() const
{
  return "ncm_image_session"; 
}

//
//: Class name equivalence
bool ncm_image_session::is_class(const vcl_string& class_str) const
{ 
  return (is_a() == class_str); 
}

//
//: Binary write
void vsl_b_write(vsl_b_ostream& bfs, const ncm_image_session& b)
{
  b.b_write(bfs);
}

//
//: Binary read
void vsl_b_read(vsl_b_istream& bfs, ncm_image_session& b)
{
  b.b_read(bfs);
}

//
//: Text write
void ncm_image_session::t_write(vcl_ostream& tfs) const
{
  tfs << vsl_indent() << is_a() << ": {\n";
  vsl_indent_inc(tfs);

  tfs << vsl_indent() << "version: " << version() << '\n';
	tfs << vcl_flush;
}

//
//: Text read
void ncm_image_session::t_read(vcl_istream& tfs)
{
  // Read only the bit between the braces - the calling function should already
  // have read the identifying token (otherwise, how would it know to call 
  // this function?)

  mbl_read_props_type props = mbl_read_props(tfs);

  short version = 
      vul_string_atoi(props.get_required_property("version"));

  switch (version)
  {
    case 1:
      {
        break;
      }

    default:
      // probably a bit extreme to abort() but will leave it for now
      vcl_cerr << is_a() << "::t_read() ";
      vcl_cerr << "Unexpected version number " << version << vcl_endl;
      abort();
  }
}

//
//: SQL write: sessions get written when they are first created.
// This means we only have values for the subject, study and user ids and the time started
bool ncm_image_session::sql_write(QSqlQuery& query)
{
	query.prepare(
		"INSERT INTO session (subject_id, user_id, time_started) "
    "VALUES (:subject_id, :user_id, :time_started)");
  query.bindValue(":subject_id", subject_id_);
	query.bindValue(":user_id", user_id_);
	query.bindValue(":time_started", time_started_);
  bool success = query.exec();

	//If successful write, get back the auto-incremented ID (this is why func can't be const)
	if (success)
		session_id_ = query.lastInsertId().toInt();
	else
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;

	return success;
}

//
//: SQL update
bool ncm_image_session::sql_update(QSqlQuery& query, const QString& fieldname, const QVariant &value) const
{
	//Can only process if sequence has already been written
	if (session_id_ < 0)
		return false; //Or be stronger and assert?

	//and if the fieldname is one of the updateable fields
	else if (sql_updateable_fields_.indexOf(fieldname) < 0)
		return false; //Or be stronger and assert?

	QString query_str = "UPDATE session SET " + fieldname + " = :" + fieldname + " WHERE session_id = " 
		+ QString::number(session_id_);
	query.prepare(query_str);
	query.bindValue(":" + fieldname, value);
  bool success = query.exec();

	//If not success pritn out the error
	if (!success)
	{
		vcl_cout << "Query: " << query_str.toStdString() << vcl_endl;
		vcl_cout << "Value: " << value.toString().toStdString() << vcl_endl;
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;
	}
	return success;
}

//
//: SQL read
void ncm_image_session::sql_read(QSqlRecord& record)
{
	//session id
	set_session_id(record.value("session_id").toInt());
	set_user_id(record.value("subject_id").toInt());
	set_subject_id(record.value("user_id").toInt());
	set_time_started(record.value("time_started").toDateTime());
	set_time_finished(record.value("time_finished").toDateTime());
	set_notes(record.value("session_comments").toString());
}

void ncm_image_session::sql_read(QSqlRecord& record, QSqlQuery& query)
{ 
	sql_read(record); //Load in subject fields
	query.prepare("SELECT * FROM image_sequence WHERE image_session_id = " + QString::number(session_id_));
	bool success = query.exec();

	if (success)
	{
		//Make sessions list is clear to start with so we don't duplicate
		sequences_.clear();
		while (query.next())
		{
			ncm_image_sequence sequence;
			sequence.sql_read(query.record());
			sequences_.push_back(sequence);
		}
	}
	else
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;
}

bool ncm_image_session::sql_delete(QSqlQuery& query)
{
	query.prepare(
		"DELETE FROM session WHERE session_id = " + QString::number(session_id_));
	bool success = query.exec();

	if (!success)
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;
	
	return success;

}

bool ncm_image_session::sql_add_to_process_list(QSqlQuery& query, bool process_bad)
{
	query.prepare(
		"INSERT INTO sequences_to_process (sequence_id) "
    "VALUES (:sequence_id)");
	
	//Loop through each sequence
	bool success = false;
	for(size_t i = 0, size = sequences_.size(); i != size; ++i)
	{
		//If we're processing bad sequences or if the sequence is good, write it to the process list
		if (process_bad || sequences_[i].acceptable_quality())
		{
			query.bindValue(":sequence_id", sequences_[i].sequence_id());
			success = query.exec();

			if (!success) //If any SQL write fails we may as well bug out now
				break;
		}

	}

	if (!success)
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;
	
	return success;

}
	

//
//: Display details of the object to standard screen output to help development
void ncm_image_session::display_details() const
{
	vcl_cout << "NCM image session object" << vcl_endl;
	if (!time_started_.isNull())
		vcl_cout << "Started: " << time_started_.toString("hh:mm:ss dd/MM/yyyy").toStdString() << vcl_endl;
	if (!time_finished_.isNull())
		vcl_cout << "Finished: " << time_finished_.toString("hh:mm:ss dd/MM/yyyy").toStdString() << vcl_endl;
	vcl_cout << "Number of sequences in session: " << sequences_.size() << vcl_endl;
}

//
//: Text output by reference
vcl_ostream& operator<<(vcl_ostream& os, const ncm_image_session& b)
{
  os << b.is_a() << ": ";
  vsl_indent_inc(os);
  b.print_summary(os);
  vsl_indent_dec(os);
  return os;
}

//
//: Text output by pointer
vcl_ostream& operator<<(vcl_ostream& os, const ncm_image_session* b)
{
  if (b)  
    return os << *b;
  else
    return os << "No " << b->is_a() << " defined.";
}

//
// Protected functions

//
//: Binary write
void ncm_image_session::b_write(vsl_b_ostream& bfs) const
{
  vsl_b_write(bfs,version());

  //// version 2
  //vsl_b_write(bfs,new_variable_);

  //// version 1
  //vsl_b_write(bfs,points_);
  //vsl_b_write(bfs,widths_);
  //vsl_b_write(bfs,apex_indices_);
  //vsl_b_write(bfs,vessel_traits_);
}

//
//: Binary read
void ncm_image_session::b_read(vsl_b_istream& bfs)
{
  short version;
  vsl_b_read(bfs,version);
  //switch (version)
  //{
  //  //case (2): 
  //  //  vsl_b_read(bfs,new_variable_);
  //  //  // no break

  //  case (1):
  //    vsl_b_read(bfs,points_);
  //    vsl_b_read(bfs,widths_);
  //    vsl_b_read(bfs,apex_indices_);
  //    vsl_b_read(bfs,vessel_traits_);
  //    break;

  //  default:
  //    vcl_cerr << "ncm_image_session::b_read() ";
  //    vcl_cerr << "Unexpected version number " << version << vcl_endl;
  //    abort();
  //}
}

//
//: Print summary to output stream os
void ncm_image_session::print_summary(vcl_ostream& os) const
{
  //os << "Points: " << points_.size() << '\n';
  //os << "Apices: " << apex_indices_.size() << '\n';
  //os << "Traits: ";
  //if (is_enlarged()) 
  //  os << "Enlarged ";
  //if (is_giant()) 
  //  os << "Giant ";
  //if (is_bushy()) 
  //  os << "Bushy ";
  //os << '\n';

  //os << vcl_flush;
}

