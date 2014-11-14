//:
// \file
// \brief Functions for Input/Output of ncm_image_sequence:
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

#include "ncm_image_sequence.h"

#include <vsl/vsl_binary_io.h>
#include <vsl/vsl_vector_io.h>
#include <vsl/vsl_indent.h>

#include <vul/vul_string.h>

#include <mbl/mbl_parse_block.h>
#include <mbl/mbl_read_props.h>
#include <mbl/mbl_parse_tuple.h>
#include <mbl/mbl_parse_sequence.h>
#include <mbl/mbl_parse_keyword_list.h>

#include <QtSql/QSqlQuery>
#include <QtSql/QSqlRecord>
#include <QtSql/QSqlError>
#include <QVariant>
#include <QDir>

//Initialize the list of fieldmames in the SQL database
const QStringList ncm_image_sequence::sql_fields_ = QStringList() << 
	"image_session_id" <<
	"hand" <<
	"digit" <<
	"num_frames" <<
	"time_started" <<
	"time_finished" <<
	"images_dir" <<
	"images_datafile" <<
	"sequence_comments" <<
	"sequence_name" <<
	"acceptable_quality" <<
	"preview_mosaic" <<
	"full_mosaic" <<
	"preview_count_map" <<
	"full_count_map";

const QStringList ncm_image_sequence::sql_updateable_fields_ = QStringList() << 
	"images_dir" <<
	"images_datafile" <<
	"sequence_comments" <<
	"sequence_name" <<
	"acceptable_quality" <<
	"preview_mosaic" <<
	"full_mosaic" <<
	"preview_count_map" <<
	"full_count_map";

//Initialize the format for pritning the date/time
const QString ncm_image_sequence::DATE_FORMAT = QString("ddd dd.MM.yyyy hh:mm:ss");


// Public functions

//
//: Version number of this class's binary IO
short ncm_image_sequence::version() const
{ 
  return 1; 
}

//
//: The class name
vcl_string ncm_image_sequence::is_a() const
{
  return "ncm_image_sequence"; 
}

//
//: Class name equivalence
bool ncm_image_sequence::is_class(const vcl_string& class_str) const
{ 
  return (is_a() == class_str); 
}

//
//: Binary write
void vsl_b_write(vsl_b_ostream& bfs, const ncm_image_sequence& b)
{
  b.b_write(bfs);
}

//
//: Binary read
void vsl_b_read(vsl_b_istream& bfs, ncm_image_sequence& b)
{
  b.b_read(bfs);
}

//:Try and sync the sequence name and the images directory
// This allows the sequence name to be selected after the images are saved
// Returns true if the operation is successful (the only reason it should fail
// is if something else has created a directory of that name before we get
// to do the sync)
bool ncm_image_sequence::sync_name_and_dir(QSqlQuery& query)
{
	QDir dir;
	bool success = dir.rename(
		images_root_dir_ + "/" + images_sub_dir_,
		images_root_dir_ + "/" + sequence_name_);

	if (success)
	{
		images_sub_dir_ = sequence_name_;
		sql_update(query, "images_dir", images_root_dir_ + "/" + images_sub_dir_);
		sql_update(query, "sequence_name", sequence_name_);
	}
	return success;
}

//
//: SQL write: sequences get written when they are finished.
// This means we have values for all fields apart from the ID which is returned from the database (and why this can't be const)
bool ncm_image_sequence::sql_write(QSqlQuery& query)
{
	query.prepare(
		"INSERT INTO image_sequence (image_session_id, hand, digit, num_frames," 
			" time_started, time_finished, images_dir, images_datafile, sequence_comments,"
			" sequence_name, acceptable_quality, preview_mosaic, full_mosaic, preview_count_map, full_count_map) "
    "VALUES(:image_session_id, :hand, :digit, :num_frames," 
			" :time_started, :time_finished, :images_dir, :images_datafile, :sequence_comments,"
			" :sequence_name, :acceptable_quality, :preview_mosaic, :full_mosaic, :preview_count_map, :full_count_map)");

	/*query.prepare(
		"INSERT INTO image_sequence (image_session_id, hand, digit, num_frames, time_started, time_finished) "
    "VALUES (:image_session_id, :hand, :digit, :num_frames, :time_started, :time_finished)");*/

  query.bindValue(":image_session_id", session_id_);
	query.bindValue(":hand", ncm_hand::toLetter(hand_).c_str());
	query.bindValue(":digit", digit_); 
	query.bindValue(":num_frames", num_frames_); 
	query.bindValue(":time_started", time_started_);
	query.bindValue(":time_finished", time_finished_);
	query.bindValue(":images_dir", images_root_dir_ + "/" + images_sub_dir_);
	query.bindValue(":images_datafile", images_datafile_);
	query.bindValue(":sequence_comments", notes_);
	query.bindValue(":sequence_name", sequence_name_);
	query.bindValue(":acceptable_quality", acceptable_quality_);
	query.bindValue(":preview_mosaic", preview_mosaic_);
	query.bindValue(":full_mosaic", full_mosaic_);
	query.bindValue(":preview_mosaic", preview_count_map_);
	query.bindValue(":full_mosaic", full_count_map_);
  bool success = query.exec();

	//If successful write, get back the auto-incremented ID (this is why func can't be const)
	if (success)
	{
		sequence_id_ = query.lastInsertId().toInt();
	}
	else
	{
		vcl_cout << query.lastQuery().toStdString() << vcl_endl;
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;
	}

	return success;
}

//
//: SQL update - can be called at anytime. However, the values for the IDs and time started can't be modified
// so only time finished and session comments can be updated
bool ncm_image_sequence::sql_update(QSqlQuery& query, const QString& fieldname, const QVariant &value) const
{
	//Can only process if sequence has already been written
	if (sequence_id_ < 0)
		return false; //Or be stronger and assert?

	//and if the fieldname is one of the updateable fields
	else if (sql_updateable_fields_.indexOf(fieldname) < 0)
		return false; //Or be stronger and assert?

	QString query_str = "UPDATE image_sequence SET " + fieldname + " = :" + fieldname + " WHERE image_sequence_id = " 
		+ QString::number(sequence_id_);
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
//: SQL read
void ncm_image_sequence::sql_read(QSqlRecord& record)
{
	sequence_id_ = record.value("image_sequence_id").toInt();
	session_id_ = record.value("image_session_id").toInt();
	hand_ = ncm_hand::fromString(record.value("hand").toString().toStdString());
	digit_ = static_cast<ncm_hand::Digit>(record.value("digit").toInt()); 
	num_frames_ = record.value("num_frames").toInt(); 
	time_started_ = record.value("time_started").toDateTime();
	time_finished_ = record.value("time_finished").toDateTime();
	set_images_full_dir( record.value("images_dir").toString() );
	images_datafile_ = record.value("images_datafile").toString();
	notes_ = record.value("sequence_comments").toString();
	sequence_name_ = record.value("sequence_name").toString();
	acceptable_quality_ = record.value("acceptable_quality").toBool();
	preview_mosaic_ = record.value("preview_mosaic").toString();
	full_mosaic_ = record.value("full_mosaic").toString();
}

//
//: Text write
void ncm_image_sequence::t_write(vcl_ostream& tfs) const
{
  tfs << vsl_indent() << is_a() << ": {" << '\n';

  vsl_indent_inc(tfs);
  tfs << vsl_indent() << "version: " << version() << '\n';

  // version 1
	tfs << vsl_indent() << "sequence_id: " << sequence_id_ << '\n';
	tfs << vsl_indent() << "sequence_name: " << sequence_name_.toStdString() << '\n';
	tfs << vsl_indent() << "time_started: " << time_started().toString(DATE_FORMAT).toStdString() << '\n';
	tfs << vsl_indent() << "time_finished: " << time_finished().toString(DATE_FORMAT).toStdString() << '\n';
	tfs << vsl_indent() << "session_id: " << session_id_ << '\n';
	tfs << vsl_indent() << "hand: " << ncm_hand::toString(hand_) << '\n';
	tfs << vsl_indent() << "digit: " << digit_ << '\n';
	tfs << vsl_indent() << "num_frames: " << num_frames_ << '\n';
	tfs << vsl_indent() << "image_dir: " << images_full_dir().toStdString() << '\n'; //This *should* always be the dir where the file is located but we won't assume this
	tfs << vsl_indent() << "acceptable_quality: " << acceptable_quality_ << '\n';
  tfs << vsl_indent() << "ncm_camera_properties: {" << '\n';
	vsl_indent_inc(tfs);
  write_camera_properties(tfs);
	vsl_indent_dec(tfs);
	tfs << vsl_indent() << "} // ncm_camera_properties" << '\n';

  tfs << vsl_indent() << "frames: {" << '\n';
  vsl_indent_inc(tfs);

  for (vcl_vector<ncm_video_frame_header>::const_iterator it = frame_headers_.begin();
       it != frame_headers_.end(); ++it)
	{
			 write_frame(tfs, (*it));
	}

  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // frames" << '\n';

	//tfs << vsl_indent() << "notes: " << notes_.toStdString() << '\n'; //MB: What if notes are multiline??

  vsl_indent_dec(tfs);
  tfs << vsl_indent() << "} // " << is_a() << '\n';
  tfs << vcl_flush;
	
}

//
//: Text read
void ncm_image_sequence::t_read(vcl_istream& tfs)
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
				sequence_id_ = vul_string_atoi(props.get_required_property("sequence_id"));
				sequence_name_ =  props.get_required_property("sequence_name").c_str();
				vcl_string time_started = props.get_required_property("time_started");
				time_started_ = QDateTime::fromString(time_started.c_str(), DATE_FORMAT);

				vcl_string time_finished = props.get_required_property("time_finished");
				time_finished_ = QDateTime::fromString(time_finished.c_str(), DATE_FORMAT);

				session_id_ = vul_string_atoi(props.get_required_property("session_id"));
				hand_ = ncm_hand::fromString(props.get_required_property("hand"));
				digit_ = ncm_hand::Digit(vul_string_atoi(props.get_required_property("digit")));
				num_frames_ = vul_string_atoi(props.get_required_property("num_frames"));
				set_images_full_dir(props.get_required_property("image_dir").c_str()); //This *should* always be the dir where the file is located but we won't assume this
				acceptable_quality_ = vul_string_to_bool(props.get_required_property("acceptable_quality"));

				//Get camera properties
				vcl_string camera_data = props.get_required_property("ncm_camera_properties");
				vcl_stringstream camera_stream(camera_data);        
        read_camera_properties(camera_stream);

				//Get motor properties
				//vcl_string motor_data = props.get_required_property("ncm_motor_properties");
				//vcl_stringstream motor_stream(camera_data);        
        //read_motor_properties(motor_stream);

        // Read frames
        vcl_string frames_data = props.get_required_property("frames");
        vcl_stringstream frames_stream(frames_data);
        vcl_vector<vcl_string> frame_strings;
        const bool discard_comments = true;
        mbl_parse_keyword_list(frames_stream, "ncm_frame:", frame_strings,
                               discard_comments);

				int frame_num = 0;
				for (vcl_vector<vcl_string>::iterator it = frame_strings.begin();
             it != frame_strings.end(); ++it)
        {
					read_frame(vcl_stringstream(*it), ++frame_num);
        }

				//QString notes_ =  props.get_required_property("notes").c_str(); //MB: What if notes are multiline??

				#if _DEBUG
					display_details();
				#endif

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
//: Text read (from file, with identifying token)
void ncm_image_sequence::t_read(const vcl_string& filename)
{
	clear();
	//Set the image datafile from filename
	images_datafile_ = QString(filename.c_str()).section("/",-1);
  vcl_ifstream tfs(filename.c_str());

  mbl_read_props_type props = mbl_read_props(tfs);
  vcl_string token = is_a();
  vcl_string sequence_string = props.get_required_property(token);
  vcl_stringstream sequence_stream(sequence_string);
  t_read(sequence_stream);

  tfs.close();
}

void ncm_image_sequence::read_camera_properties(vcl_istream& tfs)
{
  // Read only the bit between the braces - the calling function should already
  // have read the identifying token (otherwise, how would it know to call 
  // this function?)

  //clear();

  mbl_read_props_type props = mbl_read_props(tfs);

	camera_type_ = QString(props.get_required_property("camera_type").c_str());
	camera_id_ = QString(props.get_required_property("camera_id").c_str());
	frame_rate_ = vul_string_atof(props.get_required_property("frame_rate"));
	exposure_time_ = vul_string_atof(props.get_required_property("exposure_time"));
	gain_ = vul_string_atof(props.get_required_property("gain"));
}

void ncm_image_sequence::read_motor_properties(vcl_istream& tfs)
{
  // Read only the bit between the braces - the calling function should already
  // have read the identifying token (otherwise, how would it know to call 
  // this function?)

  //clear();

  mbl_read_props_type props = mbl_read_props(tfs);

	apt_x_id_ = vul_string_atoi(props.get_required_property("apt_x_id"));
	apt_y_id_ = vul_string_atoi(props.get_required_property("apt_y_id"));
	apt_z_id_ = vul_string_atoi(props.get_required_property("apt_z_id"));
	x_reversed_ = vul_string_to_bool(props.get_required_property("x_reversed"));
	y_reversed_ = vul_string_to_bool(props.get_required_property("y_reversed"));
	z_reversed_ = vul_string_to_bool(props.get_required_property("z_reversed"));
}

void ncm_image_sequence::read_frame(vcl_istream& tfs, int frame_num)
{
  // Read only the bit between the braces - the calling function should already
  // have read the identifying token (otherwise, how would it know to call 
  // this function?)

  //clear();

  mbl_read_props_type props = mbl_read_props(tfs);

  vcl_string frame_name = props.get_required_property("frame_name");
	double frame_number = vul_string_atoi(props.get_required_property("frame_number"));
	vcl_string time_str = props.get_required_property("time");
	QTime time = QTime::fromString(time_str.c_str(), "hhmmsszzz");
	double motor_x = vul_string_atof(props.get_required_property("motor_x"));
	double motor_y = vul_string_atof(props.get_required_property("motor_y"));
	double motor_z = vul_string_atof(props.get_required_property("motor_z"));
	double sharpness = vul_string_atof(props.get_required_property("sharpness"));

	ncm_video_frame_header fh(time);
  fh.frame_name_ = frame_name,
	fh.frame_num_ = frame_num;
	fh.motor_position_[0] = motor_x;
	fh.motor_position_[1] = motor_y;
	fh.motor_position_[2] = motor_z;
	fh.sharpness_ = sharpness;
	add_frame(fh);

}

void ncm_image_sequence::write_camera_properties(vcl_ostream& os) const
{
	os << vsl_indent() << "camera_type: " << camera_type_.toStdString() << '\n';
	os << vsl_indent() << "camera_id: " << camera_id_.toStdString() << '\n';
	os << vsl_indent() << "frame_rate: " << frame_rate_ << '\n';
	os << vsl_indent() << "exposure_time: " << exposure_time_ << '\n';
	os << vsl_indent() << "gain: " << gain_ << '\n';
}

void ncm_image_sequence::write_motor_properties(vcl_ostream& os) const
{
	os << vsl_indent() << "apt_x_id: " << apt_x_id_ << '\n';
	os << vsl_indent() << "apt_y_id: " << apt_y_id_ << '\n';
	os << vsl_indent() << "apt_z_id: " << apt_z_id_ << '\n';
	os << vsl_indent() << "x_reversed: " << x_reversed_ << '\n';
	os << vsl_indent() << "y_reversed: " << y_reversed_ << '\n';
	os << vsl_indent() << "z_reversed: " << z_reversed_ << '\n';
}

void ncm_image_sequence::write_frame(vcl_ostream& os, ncm_video_frame_header fh) const
{
	os << vsl_indent() << "ncm_frame: {" << '\n';
	vsl_indent_inc(os);

	os << vsl_indent() << "frame_name: " << fh.frame_name_ << '\n';
	os << vsl_indent() << "frame_number: " << fh.frame_num_ << '\n';
	os << vsl_indent() << "time: " << fh.frame_time_.toString("hhmmsszzz").toStdString() << '\n';
	os << vsl_indent() << "motor_x: " << fh.motor_position_[0] << '\n';
	os << vsl_indent() << "motor_y: " << fh.motor_position_[1] << '\n';
	os << vsl_indent() << "motor_z: " << fh.motor_position_[2] << '\n';
	os << vsl_indent() << "sharpness: " << fh.sharpness_ << '\n';

	vsl_indent_dec(os);
	os << vsl_indent() << "} // ncm_frame" << '\n';
}

//
//: Display details of the object to standard screen output to help development
void ncm_image_sequence::display_details() const
{ 
	vcl_cout << "NCM image sequence object" << vcl_endl;
  vcl_cout << "Version: " << version() << vcl_endl;
	vcl_cout << "Sequence ID: " << sequence_id_ << vcl_endl;
	vcl_cout << "Sequence name: " << sequence_name_.toStdString() << vcl_endl;
	vcl_cout << "Session ID: " << session_id_ << vcl_endl;
	vcl_cout << "Hand: " << ncm_hand::toLetter(hand_) << vcl_endl;
	vcl_cout << "Digit: " << digit_ << vcl_endl;
	vcl_cout << "Num frames: " << num_frames_ << vcl_endl;
	vcl_cout << "Acceptable quality: " << acceptable_quality_ << vcl_endl;
	vcl_cout << "Camera type: " << camera_type_.toStdString() << vcl_endl;
	vcl_cout << "Camera ID: " << camera_id_.toStdString() << vcl_endl;
	vcl_cout << "Frame rate: " << frame_rate_ << vcl_endl;
	vcl_cout << "Exposure time: " << exposure_time_ << vcl_endl;
	vcl_cout << "Gain: " << gain_ << vcl_endl;

	if (!time_started_.isNull())
		vcl_cout << "Started: " << time_started_.toString("hh:mm:ss dd/MM/yyyy").toStdString() << vcl_endl;
	if (!time_finished_.isNull())
		vcl_cout << "Finished: " << time_finished_.toString("hh:mm:ss dd/MM/yyyy").toStdString() << vcl_endl;

	vcl_cout << "Image directory: " << images_full_dir().toStdString() << vcl_endl;
	vcl_cout << "Images datafile: " << images_datafile_.toStdString() << vcl_endl;
	vcl_cout << "Comments: " << notes_.toStdString() << vcl_endl;

}

//
//: Text output by reference
vcl_ostream& operator<<(vcl_ostream& os, const ncm_image_sequence& b)
{
  os << b.is_a() << ": ";
  vsl_indent_inc(os);
  b.print_summary(os);
  vsl_indent_dec(os);
  return os;
}

//
//: Text output by pointer
vcl_ostream& operator<<(vcl_ostream& os, const ncm_image_sequence* b)
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
void ncm_image_sequence::b_write(vsl_b_ostream& bfs) const
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
void ncm_image_sequence::b_read(vsl_b_istream& bfs)
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
  //    vcl_cerr << "ncm_image_sequence::b_read() ";
  //    vcl_cerr << "Unexpected version number " << version << vcl_endl;
  //    abort();
  //}
}

//
//: Print summary to output stream os
void ncm_image_sequence::print_summary(vcl_ostream& os) const
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

