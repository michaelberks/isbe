#include "ncm_qdata_saver.h"

#include <vcl_iostream.h>
#include <vcl_iomanip.h>
#include <vcl_sstream.h>
#include <vcl_iosfwd.h>

#include <vil/vil_save.h>

#include <nailfold/qtools/ncm_video_frame.h>
#include <nailfold/qtools/ncm_video_frame_queue.h>

// Initialize all static variables.
// These are shared among all data streams.
vcl_string ncm_qdata_saver::root_dir_ = "";
vcl_string ncm_qdata_saver::patient_name_ = "";
vcl_string ncm_qdata_saver::sequence_name_ = "";
vcl_string ncm_qdata_saver::hand_ = "Left";
vcl_string ncm_qdata_saver::digit_ = "4";
vcl_string ncm_qdata_saver::mag_ = "x300";
	

ncm_qdata_saver::ncm_qdata_saver()
:	sub_dir_(""),
  frame_prefix_("frame_"),
  frame_suffix_(""),
  queue_(NULL),
  is_saving_(false),
  properties_as_subdir_(false)
{
}
 
void ncm_qdata_saver::setRootDir( vcl_string root_dir )
{
	root_dir_ = root_dir;
}

void ncm_qdata_saver::setPatientName( vcl_string patient_name )
{
	patient_name_ = patient_name;
}

void ncm_qdata_saver::setSequenceName( vcl_string sequence_name )
{
	sequence_name_ = sequence_name;
}

void ncm_qdata_saver::setHand( vcl_string hand )
{
	hand_ = hand;
}

void ncm_qdata_saver::setDigit( vcl_string digit )
{
	digit_ = digit;
}

void ncm_qdata_saver::setMag( vcl_string mag )
{
	mag_ = mag;
}

void ncm_qdata_saver::setSubDir( vcl_string sub_dir )
{
	sub_dir_ = sub_dir;
}

void ncm_qdata_saver::setFramePrefix( vcl_string frame_prefix )
{
	frame_prefix_ = frame_prefix;
}

void ncm_qdata_saver::setFrameSuffix( vcl_string frame_suffix )
{
	frame_suffix_ = frame_suffix;
}

void ncm_qdata_saver::set_properties_as_subdir( 
  bool properties_as_subdir /* = true */)
{
	properties_as_subdir_ = properties_as_subdir;
}


//: Define the queue from which to pop frames.
void ncm_qdata_saver::set_queue(ncm_video_frame_queue* queue)
{
  queue_ = queue;
}

 
//: Create the output folder if need be.
vcl_string ncm_qdata_saver::makeSaveDir()
{
	// Check if base_name is a directory and if not create it
	QString dir_name = QString::fromStdString( output_dir() );

	QDir dir;
	if ( !dir.exists( dir_name ) )
		dir.mkpath( dir_name );

	return output_dir();
}
 
//: Return true if there are more frames to save.
bool ncm_qdata_saver::is_saving() const
{
  return is_saving_;
}
 
//
// Public slots
//
 
bool ncm_qdata_saver::save_frame(int frame_num)
{
  if (queue_ == NULL)
    return false;

  is_saving_ = true;

	vcl_stringstream ss;

	QSharedPointer<ncm_video_frame> vid_frame = queue_->dequeue();

  // Create a timestamp string.
	QTime frame_time = vid_frame->frame_time();
	ss << vcl_setw(2) << vcl_setfill('0') << frame_time.hour() << "_"
		 << vcl_setw(2) << vcl_setfill('0') << frame_time.minute() << "_"
		 << vcl_setw(2) << vcl_setfill('0') << frame_time.second() << "_" 
		 << vcl_setw(3) << vcl_setfill('0') << frame_time.msec() << "_"
		 << vcl_setw(4) << vcl_setfill('0') << frame_num;
  vcl_string timestamp = ss.str();

  // Define full filename for the image.
	ss << output_dir()
     << hand_ << "_"
     << "D" << digit_ << "_"
     << "M" << mag_;
  
  if (properties_as_subdir_)
    ss << "/";
  else
    ss << "_";

  ss << frame_prefix_
     << frame_suffix_ 
     << timestamp << ".png";

  vcl_string full_filename = ss.str();

	// Save image stored in save queue.
	vil_save( *vid_frame->frame(), full_filename.c_str() );

  emit frame_saved(frame_num);

  if (vid_frame->is_last_frame())
  {
    is_saving_ = false;
    emit saving_finished();
  }

  return true;
}
 
//
// Private methods
//
 
//: Return the string containing the output folder.
vcl_string ncm_qdata_saver::output_dir()
{
  vcl_stringstream ss;

  ss << root_dir_ << "/";

  if (!patient_name_.empty())
    ss << patient_name_ << "/";

  if (!sequence_name_.empty())
    ss << sequence_name_ << "/";

  // Now in filename of every frame
  //ss << hand_ << "."
  //   << "Digit" << digit_ << "."
  //   << mag_ << "/";

  if (!sub_dir_.empty())
    ss << sub_dir_ << "/";

  return ss.str();
}
