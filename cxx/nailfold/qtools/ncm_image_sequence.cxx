#include "ncm_image_sequence.h"

//: Define class member functions

//: Constructors
ncm_image_sequence::ncm_image_sequence(int session_id)//const ncm_image_session &parent_session
	: session_id_(session_id),
	  hand_(ncm_hand::NotKnown),
	  digit_(ncm_hand::None),
		images_root_dir_(""),
		images_sub_dir_(""),
		sequence_name_(""),
		num_frames_(0),
		frame_headers_(0),
		notes_(""),
		sequence_id_(-1),
		acceptable_quality_(true),
		camera_type_(""),
		camera_id_(""),
		frame_rate_(0),
		exposure_time_(-1),
		gain_(-1),
		apt_x_id_(-1),
		apt_y_id_(-1),
		apt_z_id_(-1),
		x_reversed_(true),
		y_reversed_(false),
		z_reversed_(false),
		preview_mosaic_(""),
		full_mosaic_("")
{
}

ncm_image_sequence::~ncm_image_sequence()
{
}

void ncm_image_sequence::start()
{
	if (time_started_.isNull()) //Can only be called once
		time_started_ = QDateTime::currentDateTime();
}

void ncm_image_sequence::finalise()
{
	if (time_finished_.isNull()) //Can only be called once
		time_finished_ = QDateTime::currentDateTime();
}

void ncm_image_sequence::clear()
{
	session_id_ = -1;
	hand_ = ncm_hand::NotKnown;
	digit_ = ncm_hand::None;
	images_datafile_ = "undefined";
	num_frames_ = 0;
	frame_headers_.clear();
	notes_ = "";
	sequence_name_ = "";
	sequence_id_ = -1;
	acceptable_quality_ = true;

	images_root_dir_ = "";
	images_sub_dir_ = "";
	sequence_name_ = "";

	time_started_ = QDateTime();
	time_finished_ = QDateTime();

}

//Set functions
void ncm_image_sequence::set_hand(ncm_hand::Hand hand)
{
	hand_ = hand;
}

void ncm_image_sequence::set_digit(ncm_hand::Digit digit)
{
	digit_ = digit;
}

void ncm_image_sequence::set_images_root_dir(QString dirpath)
{
	images_root_dir_ = dirpath;
}

void ncm_image_sequence::set_images_sub_dir(QString dirpath)
{
	images_sub_dir_ = dirpath;
}

void ncm_image_sequence::set_images_full_dir(QString dirpath)
{
	images_root_dir_ = dirpath.section("/",0,-2);
	images_sub_dir_ = dirpath.section("/", -1);
}

void ncm_image_sequence::set_images_datafile(QString filename)
{
	images_datafile_ = filename;
}

void ncm_image_sequence::set_sequence_name(QString name)
{
	sequence_name_ = name;
}

void ncm_image_sequence::set_num_frames(int n_frames)
{
	num_frames_ = n_frames;
}
void ncm_image_sequence::set_notes(QString notes)
	{
	notes_ = notes;
}

void ncm_image_sequence::set_acceptable_quality(bool yesno)
{
	acceptable_quality_ = yesno;
}

void ncm_image_sequence::set_camera_type(QString camera)
{
	camera_type_ = camera;
}
void ncm_image_sequence::set_camera_id(QString id)
{
	camera_id_ = id;
}
void ncm_image_sequence::set_frame_rate(double fps)
{
	frame_rate_ = fps;
}
void ncm_image_sequence::set_exposure_time(double t)
{
	exposure_time_ = t;
}
void ncm_image_sequence::set_gain(double g)
{
	gain_ = g;
}

void ncm_image_sequence::set_apt_x_id(long id)
{
	apt_x_id_ = id;
}
void ncm_image_sequence::set_apt_y_id(long id)
{
	apt_y_id_ = id;
}
void ncm_image_sequence::set_apt_z_id(long id)
{
	apt_z_id_ = id;
}
void ncm_image_sequence::set_x_reversed(bool rev)
{
	x_reversed_ = rev;
}
void ncm_image_sequence::set_y_reversed(bool rev)
{
	y_reversed_ = rev;
}
void ncm_image_sequence::set_z_reversed(bool rev)
{
	z_reversed_ = rev;
}

void ncm_image_sequence::set_preview_mosaic(QString filename)
{
	preview_mosaic_ = filename;
}

void ncm_image_sequence::set_full_mosaic(QString filename)
{
	full_mosaic_ = filename;
} 

void ncm_image_sequence::set_preview_count_map(QString filename)
{
	preview_count_map_ = filename;
} 

void ncm_image_sequence::set_full_count_map(QString filename)
{
	full_count_map_ = filename;
} 

//Get functions
ncm_hand::Hand ncm_image_sequence::hand() const
{
	return hand_;
}

ncm_hand::Digit ncm_image_sequence::digit() const
{
	return digit_;
}

QString ncm_image_sequence::notes() const
	{
	return notes_;
}

int ncm_image_sequence::num_frames() const
{
	return num_frames_;
}

QString ncm_image_sequence::images_root_dir() const
{
	return images_root_dir_;
}

QString ncm_image_sequence::images_sub_dir() const
{
	return images_sub_dir_;
}

QString ncm_image_sequence::images_full_dir() const
{
	return images_root_dir_ + "/" + images_sub_dir_;
}

QString ncm_image_sequence::images_datafile() const
	{
	return images_datafile_;
}

QString ncm_image_sequence::sequence_name() const
{
	return sequence_name_;
}

bool ncm_image_sequence::acceptable_quality() const
{
	return acceptable_quality_;
}

QString ncm_image_sequence::camera_type() const
{
	return camera_type_;
}

QString ncm_image_sequence::camera_id() const
{
	return camera_id_;
}

double ncm_image_sequence::frame_rate() const
{
	return frame_rate_;
}

double ncm_image_sequence::exposure_time() const
{
	return exposure_time_;
}

double ncm_image_sequence::gain() const
{
	return gain_;
}

int ncm_image_sequence::sequence_id() const
{
	return sequence_id_;
}

long ncm_image_sequence::apt_x_id() const
{
	return apt_x_id_;
}
long ncm_image_sequence::apt_y_id() const
{
	return apt_y_id_;
}
long ncm_image_sequence::apt_z_id() const
{
	return apt_z_id_;
}
bool ncm_image_sequence::x_reversed() const
{
	return x_reversed_;
}
bool ncm_image_sequence::y_reversed() const
{
	return y_reversed_;
}
bool ncm_image_sequence::z_reversed() const
{
	return z_reversed_;
}

QString ncm_image_sequence::preview_mosaic() const
{
	return preview_mosaic_;
}
QString ncm_image_sequence::full_mosaic() const
{
	return full_mosaic_;
}

QString ncm_image_sequence::preview_count_map() const
{
	return preview_count_map_;
}
QString ncm_image_sequence::full_count_map() const
{
	return full_count_map_;
}

ncm_video_frame_header ncm_image_sequence::frame_header(int i) const
{
	if (i < num_frames_)
    return frame_headers_[ i ];
  else
	{
		//Return an empty header
		ncm_video_frame_header fh;
    return fh;
	}
}

//
//: Return time session was started/finished
QDateTime ncm_image_sequence::time_started() const
{
	return time_started_;
}

QDateTime ncm_image_sequence::time_finished() const
{
	return time_finished_;
}

void ncm_image_sequence::add_frame(ncm_video_frame_header frame)
{
	frame_headers_.push_back(frame);
	num_frames_ = frame_headers_.size();
}

//: Private functions
void ncm_image_sequence::set_time_started(QDateTime time)
{
	if (time_started_.isNull())
		time_started_ = time;
}

void ncm_image_sequence::set_time_finished(QDateTime time)
{
	if (time_finished_.isNull())
		time_finished_ = time;
}

void ncm_image_sequence::set_sequence_id(int id)
{
	sequence_id_ = id;
}