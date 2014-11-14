#include "ncm_image_session.h"
#include "nailfold/ncm_hand.h"

//: Define class member functions

ncm_image_session::ncm_image_session(int subject_id, int user_id)
	: subject_id_(subject_id),
		user_id_(user_id),
		sequences_(0),
		notes_("")
{
}

ncm_image_session::~ncm_image_session()
{
}

//Reset all the fields to zero
void ncm_image_session::reset()
{
	subject_id_ = -1;
	user_id_ = -1;
	sequences_.clear();
	notes_ = "";
	time_started_  = QDateTime();
	time_finished_ = QDateTime();
}

//Return the session ID
int ncm_image_session::session_id() const
{
	return session_id_;
}

//Return the image sequences captured for this session
vcl_vector<ncm_image_sequence>* ncm_image_session::image_sequences()
{
	return &sequences_;
}

void ncm_image_session::set_notes(QString notes)
	{
	notes_ = notes;
}

QString ncm_image_session::notes() const
	{
	return notes_;
}

//
//: Return time session was started/finished
QDateTime ncm_image_session::time_started() const
{
	return time_started_;
}

QDateTime ncm_image_session::time_finished() const
{
	return time_finished_;
}

void ncm_image_session::start()
{
	//This should only be called once, so if time_started isn't NULL return
	if (!time_started_.isNull())
		return; //Could be stronger and assert this to catch miscontrol in the GUI??

	time_started_ = QDateTime::currentDateTime();
}

void ncm_image_session::finalise()
{
	//This should only be called once, so if time_finished isn't NULL return
	if (!time_finished_.isNull())
		return; //Could be stronger and assert this to catch miscontrol in the GUI??

	time_finished_ = QDateTime::currentDateTime();
}

void ncm_image_session::add_image_sequence(const ncm_image_sequence sequence)
{
	sequences_.push_back(sequence);
}

QString ncm_image_session::sequences_string() const
{
	QString sequence_string = "<p>";
	for ( vcl_vector<ncm_image_sequence>::const_iterator itr = sequences_.begin(),
		end = sequences_.end(); itr != end; ++itr )
	{
		QString hand = QString(ncm_hand::toLetter((*itr).hand()).c_str());
		QString digit = QString::number((*itr).digit());
		QString nframes = QString::number((*itr).num_frames());
		QString time = (*itr).time_finished().toString("hh:mm:ss");
		if ((*itr).acceptable_quality())
			sequence_string +=  ("<b>" + hand + digit + ": " + nframes + " frames, " + time + " (Ok) </b> <br>");
			//sequence_string +=  ("<P>" + hand + digit + ": " + nframes + " frames, " + time + " (Ok) </P>");
		else
			sequence_string +=  (hand + digit + ": " + nframes + " frames, " + time + " (Bad) <br>");
			//sequence_string +=  ("<P><FONT COLOR='#ff0000'>" + hand + digit + ": " + nframes + " frames, " + time + " (Bad) </P>");
	}	
	sequence_string += "</p>";
	return sequence_string;
}

//Private functions
void ncm_image_session::set_session_id(int id)
{
	session_id_ = id;
}
void ncm_image_session::set_user_id(int id)
{
	user_id_ = id;
}
void ncm_image_session::set_subject_id(int id)
{
	subject_id_ = id;
}
void ncm_image_session::set_time_started(QDateTime time)
{
	time_started_ = time;
}
void ncm_image_session::set_time_finished(QDateTime time)
{
	time_finished_ = time;
}