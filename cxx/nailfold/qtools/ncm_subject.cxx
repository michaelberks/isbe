#include "ncm_subject.h"
#include <vcl_iostream.h>
//: Define class member functions

//: Constructors
ncm_subject::ncm_subject()
: name_(""),
	ncm_id_(-1),
	study_id_(-1),
	name_in_study_(""),
	dominant_hand_(ncm_hand::NotKnown),
	image_sessions_(0),
	notes_("")
{
}

ncm_subject::~ncm_subject()
{
	//We may want to do some tidying up of the database?
}

//Reset all the fields back to default state
void ncm_subject::reset()
{
	set_name("");
	set_ncm_id(-1);
	set_study_id(-1);
	set_name_in_study("");
	set_dominant_hand(ncm_hand::NotKnown);
	set_notes("");
	live_session_.reset();
	image_sessions_.clear();
}

//
//: Return subject name
QString ncm_subject::name() const
{
  return name_;
}

int ncm_subject::ncm_id() const{
  return ncm_id_;
}



int ncm_subject::study_id() const{
  return study_id_;
}

QString ncm_subject::name_in_study() const
{
  return name_in_study_;
}

ncm_hand::Hand ncm_subject::dominant_hand() const
{
  return dominant_hand_;
}
QString ncm_subject::notes() const
{
  return notes_;
}

QDateTime ncm_subject::first_imaged() const
{
	if (image_sessions_.empty())
	{
		//Return a NULL date time
		QDateTime dt;
		return dt;
	}
	else
		return image_sessions_.front().time_started();
}

QDateTime ncm_subject::last_imaged() const
{
  if (image_sessions_.empty())
	{
		//Return a NULL date time
		QDateTime dt;
		return dt;
	}
	else
		return image_sessions_.back().time_started();
}

//Return string compunding information on previous sessions into displayable form
QString ncm_subject::previous_sessions_string() const
{
	QString previous_sessions;
	for ( vcl_vector<ncm_image_session>::const_iterator itr = image_sessions_.begin(),
		end = image_sessions_.end(); itr != end; ++itr )
	{
		previous_sessions += ((*itr).time_started().toString("ddd MMMM d yyyy, hh:mm") + "\n");
	}	
	return previous_sessions;
}

//: Return pointer to live session object
ncm_image_session* ncm_subject::live_session()
{
	return &live_session_;
}

//: Return pointer to session selected by idx
/*ncm_image_session* ncm_subject::select_session(int idx)
{
	int i = 0;
	for ( vcl_vector<ncm_image_session>::iterator itr = image_sessions_.begin(),
		end = image_sessions_.end(); itr != end; ++itr )
	{
		if (i == idx)
			return &(*itr);
		else
			++i;
	}
	return NULL;
}*/

ncm_image_session ncm_subject::select_session(int idx)
{
	assert(idx < image_sessions_.size());
	return image_sessions_[idx];
}

//: Return the session index of session with given id. If session id not found returns -1
int ncm_subject::find_session(int id)
{
	int session_idx = 0;
	for ( vcl_vector<ncm_image_session>::const_iterator itr = image_sessions_.begin(),
		end = image_sessions_.end(); itr != end; ++itr )

	{
		if (itr->session_id() == id)
			return session_idx;

		++session_idx;
	}

	//If we got this far the session wasn't found
	return -1;
}

//Set subject name
void ncm_subject::set_name(QString name)
{
	name_ = name;
}

//Set subject study name
void ncm_subject::set_name_in_study(QString name)
{
	name_in_study_ = name;
}

void ncm_subject::set_ncm_id(int id) 
{
  ncm_id_ = id;
}

void ncm_subject::set_study_id(int id) 
{
  study_id_ =id;
}

//Set subject's dominant hand
void ncm_subject::set_dominant_hand(ncm_hand::Hand hand)
{
	dominant_hand_ = hand;
}

//Set subject notes
void ncm_subject::set_notes(QString notes)
{
	notes_ = notes;
}


void ncm_subject::set_live_session(ncm_image_session session)
{
	live_session_ = session;
}

void ncm_subject::add_live_to_previous()
{
	image_sessions_.push_back(live_session_);
}


