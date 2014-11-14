#include "ncm_qcapture_saver.h"

#include <vcl_iostream.h>
#include <vcl_iomanip.h>
#include <vcl_sstream.h>
#include <vcl_iosfwd.h>
#include <vil/vil_save.h>

#include <nailfold/qtools/ncm_video_frame.h>
#include "ncm_qcapture_data_manager.h"

ncm_qcapture_saver::ncm_qcapture_saver()
:	root_dir_(""),
	hand_("Left"),
	digit_("4"),
	mag_("x300"),
	study_name_(""),
	subject_name_(""),
	session_name_(""),
	sequence_name_(""),
  frame_prefix_("frame"),
  is_saving_(false)
{
	update_frame_name();
}
 
void ncm_qcapture_saver::setHand( vcl_string hand )
{
	hand_ = hand;
	update_frame_name();
}

void ncm_qcapture_saver::setDigit( vcl_string digit )
{
	digit_ = digit;
	update_frame_name();
}

void ncm_qcapture_saver::setMag( vcl_string mag )
{
	mag_ = mag;
	update_frame_name();
}

void ncm_qcapture_saver::setRootDir( vcl_string root_dir )
{
	root_dir_ = root_dir;
	update_frame_name();
}

void ncm_qcapture_saver::setStudyName( vcl_string name )
{
	study_name_ = name;
	update_frame_name();
}

void ncm_qcapture_saver::setSubjectName( vcl_string name )
{
	subject_name_ = name;
	update_frame_name();
}

void ncm_qcapture_saver::setSessionName( vcl_string name )
{
	session_name_ = name;
	update_frame_name();
}

void ncm_qcapture_saver::setSequenceName( vcl_string name )
{
	sequence_name_ = name;
	update_frame_name();
}

vcl_string ncm_qcapture_saver::makeSaveDir()
{
	//Check if frame_dir_ is a directory and if not create it
	QString dir_name = QString::fromStdString( frame_dir_ );
	QDir dir;
	if ( !dir.exists( dir_name ) )
		dir.mkpath( dir_name );

	return frame_dir_;
}
 
//: Return true if there are more frames to save.
bool ncm_qcapture_saver::is_saving() const
{
  return is_saving_;
}
 
//
// Public slots
//
 
void ncm_qcapture_saver::save_frame(int frame_num)
{
  is_saving_ = true;

	// Pop frame from main queue
	QSharedPointer< ncm_video_frame > vid_frame = 
      ncm_qcapture_data_manager::Instance()->remove_from_save_queue();

	/*QTime frame_time = vid_frame->frame_time();

	// Make name for this frame
	vcl_stringstream ss;
	ss    << frame_dir_ << frame_name_ << "_"
				<< vcl_setw( 2 ) << vcl_setfill( '0' ) << frame_time.hour() << "_"
				<< vcl_setw( 2 ) << vcl_setfill( '0' ) << frame_time.minute() << "_"
				<< vcl_setw( 2 ) << vcl_setfill( '0' ) << frame_time.second() << "_" 
				<< vcl_setw( 3 ) << vcl_setfill( '0' ) << frame_time.msec() << "_"
				<< vcl_setw( 4 ) << vcl_setfill( '0' ) << frame_num << ".png";*/

	// Save image stored in mem_buffer
	vcl_string frame_name = frame_dir_ + "/" + vid_frame->frame_name();
	vil_save( *vid_frame->frame(), frame_name.c_str() );
  emit frame_saved(frame_num);

  if (vid_frame->is_last_frame())
  {
    is_saving_ = false;
    emit saving_finished();
  }
}
 
//
// Private methods
//
void ncm_qcapture_saver::update_frame_name()
{
	frame_dir_ = root_dir_;
	frame_dir_.append( "/" );
	frame_dir_.append( study_name_ );
	frame_dir_.append( "/" );
	frame_dir_.append( subject_name_ );
	frame_dir_.append( "/" );
	frame_dir_.append( session_name_ );
	frame_dir_.append( "/" );
	frame_dir_.append( sequence_name_ );
	
}