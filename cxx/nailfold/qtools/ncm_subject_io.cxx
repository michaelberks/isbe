//:
// \file
// \brief Functions for Input/Output of ncm_subject:
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

#include "ncm_subject.h"

#include <vsl/vsl_binary_io.h>
#include <vsl/vsl_vector_io.h>
#include <vsl/vsl_indent.h>

#include <vul/vul_string.h>

#include <mbl/mbl_parse_block.h>
#include <mbl/mbl_read_props.h>
#include <mbl/mbl_parse_tuple.h>
#include <mbl/mbl_parse_sequence.h>

#include <QVariant>
#include <QtSql/QSqlQuery>
#include <QtSql/QSqlRecord>
#include <QtSql/QSqlError>

//Initialize the list of fieldmames in the SQL database
const QStringList ncm_subject::sql_fields_ = QStringList() << 
	"subject_study_id" <<
	"subject_study_name" <<
	"dominant_hand" <<
	"notes";

const QStringList ncm_subject::sql_updateable_fields_ = QStringList() << 
	"dominant_hand" <<
	"notes";

//
// Public functions

//
//: Version number of this class's binary IO
short ncm_subject::version() const
{ 
  return 1; 
}

//
//: The class name
vcl_string ncm_subject::is_a() const
{
  return "ncm_subject"; 
}

//
//: Class name equivalence
bool ncm_subject::is_class(const vcl_string& class_str) const
{ 
  return (is_a() == class_str); 
}

//
//: Binary write
void vsl_b_write(vsl_b_ostream& bfs, const ncm_subject& b)
{
  b.b_write(bfs);
}

//
//: Binary read
void vsl_b_read(vsl_b_istream& bfs, ncm_subject& b)
{
  b.b_read(bfs);
}

//
//: Text write
void ncm_subject::t_write(vcl_ostream& tfs) const
{
  tfs << vsl_indent() << is_a() << ": {\n";
  vsl_indent_inc(tfs);

  tfs << vsl_indent() << "version: " << version() << '\n';
	tfs << vcl_flush;
}

//
//: Text read
void ncm_subject::t_read(vcl_istream& tfs)
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
//: SQL write
bool ncm_subject::sql_write(QSqlQuery& query)
{
	query.prepare("INSERT INTO subject (subject_study_id, subject_study_name, dominant_hand, notes) "
                   "VALUES (:subject_study_id, :subject_study_name, :dominant_hand, :notes)");
  query.bindValue(":subject_study_id", study_id_);
  query.bindValue(":subject_study_name", name_in_study_);
	query.bindValue(":dominant_hand", ncm_hand::toLetter(dominant_hand_).c_str());
	query.bindValue(":notes", notes_);
  bool success = query.exec();

	//If successful write, get back the auto-incremented ID (this is why func can't be const)
	if (success)
		ncm_id_ = query.lastInsertId().toInt();
	else
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;

	return success;
}

//
//: SQL read from record, won't try and load in previous sessions
void ncm_subject::sql_read(QSqlRecord& record)
{
	//NCM ID
	set_ncm_id( record.value("ncm_subject_id").toInt() );

	//Study ID
	set_study_id( record.value("subject_study_id").toInt() );

	//Study name
	set_name_in_study( record.value("subject_study_name").toString() );

	//Dominant hand
	set_dominant_hand( ncm_hand::fromString(record.value("dominant_hand").toString().toStdString()) );

	//Notes
	set_notes( record.value("notes").toString() );

}
//
//: SQL read from record, and load associated sessions from database
void ncm_subject::sql_read(QSqlRecord& record, QSqlQuery &query)
{
	sql_read(record); //Load in subject fields
	query.prepare("SELECT * FROM session WHERE subject_id = " + QString::number(ncm_id_));
	bool success = query.exec();

	if (success)
	{
		//Make sessions list is clear to start with so we don't duplicate
		image_sessions_.clear();
		while (query.next()) 
		{
			ncm_image_session session;
			session.sql_read(query.record());
			image_sessions_.push_back(session);
		}
	}
	else
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;

}

//
//: SQL update
bool ncm_subject::sql_update(QSqlQuery& query, const QString& fieldname, const QVariant &value) const
{
	//Can only process if ncm has already been written
	if (ncm_id_ < 0)
		return false; //Or be stronger and assert?

	//and if the fieldname is one of the updateable fields
	else if (sql_updateable_fields_.indexOf(fieldname) < 0)
		return false; //Or be stronger and assert?

	QString query_str = "UPDATE session SET " + fieldname + " = :" + fieldname + " WHERE ncm_subject_id = " 
		+ QString::number(ncm_id_);
	query.prepare(query_str);
	query.bindValue(":" + fieldname, value);
  bool success = query.exec();

	//If successful write, get back the auto-incremented ID (this is why func can't be const)
	if (!success)
	{
		vcl_cout << "Query: " << query_str.toStdString() << vcl_endl;
		vcl_cout << "Value: " << value.toString().toStdString() << vcl_endl;
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;
	}
	return success;
}


//
//: Text output by reference
vcl_ostream& operator<<(vcl_ostream& os, const ncm_subject& b)
{
  os << b.is_a() << ": ";
  vsl_indent_inc(os);
  b.print_summary(os);
  vsl_indent_dec(os);
  return os;
}

//
//: Text output by pointer
vcl_ostream& operator<<(vcl_ostream& os, const ncm_subject* b)
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
void ncm_subject::b_write(vsl_b_ostream& bfs) const
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
void ncm_subject::b_read(vsl_b_istream& bfs)
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
  //    vcl_cerr << "ncm_subject::b_read() ";
  //    vcl_cerr << "Unexpected version number " << version << vcl_endl;
  //    abort();
  //}
}

//
//: Print summary to output stream os
void ncm_subject::print_summary(vcl_ostream& os) const
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

