#include "ncm_qcapture2_saver.h"


ncm_qcapture2_saver::ncm_qcapture2_saver()
:	root_dir_(""),
	sequence_name_(""),
  hand_("Left"),
	digit_("4"),
	mag_("x300"),
	sub_dir_(""),
	patient_name_("PatientX"),
	frame_prefix_("frame"),
	camera_suffix_(""),
  is_saving_(false)
{
	update_frame_name();
}

ncm_qcapture2_saver::~ncm_qcapture2_saver()
{

}

void ncm_qcapture2_saver::setHand( vcl_string hand )
{
	vcl_cout << "Set hand activated" << vcl_endl;
	hand_ = hand;
	//update_frame_name();
}
void ncm_qcapture2_saver::setDigit( vcl_string digit )
{
	digit_ = digit;
	//update_frame_name();
}
void ncm_qcapture2_saver::setMag( vcl_string mag )
{
	mag_ = mag;
	//update_frame_name();
}
void ncm_qcapture2_saver::setSequenceName( vcl_string sequence_name )
{
	sequence_name_ = sequence_name;
	//update_frame_name();
}
void ncm_qcapture2_saver::setRootDir( vcl_string root_dir )
{
	root_dir_ = root_dir;
	//update_frame_name();
}

void ncm_qcapture2_saver::setSubDir( vcl_string sub_dir )
{
	sub_dir_ = sub_dir;
	//update_frame_name();
}

void ncm_qcapture2_saver::setPatientName( vcl_string patient_name )
{
	patient_name_ = patient_name;
	//update_frame_name();
}

void ncm_qcapture2_saver::setCameraSuffix( vcl_string camera_suffix )
{
	camera_suffix_ = camera_suffix;
	//update_frame_name();
}

vcl_string ncm_qcapture2_saver::makeSaveDir()
{
	//Check if frame_dir_ is a directory and if not create it
	QString dir_name = QString::fromStdString( frame_dir_ );
	QDir dir;
	if ( !dir.exists( dir_name ) )
		dir.mkpath( dir_name );

	return frame_dir_;
}
 
//: Return true if there are more frames to save.
bool ncm_qcapture2_saver::is_saving() const
{
  return is_saving_;
}

//
// Public slots
//

void ncm_qcapture2_saver::save_frame(int frame_num, int cam_num)
{
  is_saving_ = true;

	// Pop frame from main queue
	QSharedPointer< ncm_video_frame > vid_frame = 
      ncm_qcapture2_data_manager::Instance()->remove_from_save_queue( cam_num );

	QTime frame_time = vid_frame->frame_time();

	// Make name for this frame
	vcl_stringstream ss;
	ss    << frame_dir_ << frame_name_ << "_"
				<< vcl_setw( 2 ) << vcl_setfill( '0' ) << frame_time.hour() << "_"
				<< vcl_setw( 2 ) << vcl_setfill( '0' ) << frame_time.minute() << "_"
				<< vcl_setw( 2 ) << vcl_setfill( '0' ) << frame_time.second() << "_" 
				<< vcl_setw( 3 ) << vcl_setfill( '0' ) << frame_time.msec() << "_"
				<< vcl_setw( 4 ) << vcl_setfill( '0' ) << frame_num << ".png";

	// Save image stored in mem_buffer
	vil_save( *vid_frame->frame(), ss.str().c_str() );
  emit frame_saved(frame_num, cam_num);

  if (vid_frame->is_last_frame())
  {
    is_saving_ = false;
    emit saving_finished( cam_num );
  }
}

//
// Private methods
//
void ncm_qcapture2_saver::update_frame_name()
{
	frame_dir_ = root_dir_;
	frame_dir_.append( "/" );
	frame_dir_.append( patient_name_ );
	frame_dir_.append( "/" );
	frame_dir_.append( sequence_name_ );
	frame_dir_.append( "/" );

	frame_name_ = frame_prefix_;
	frame_name_.append( "_" );
	frame_name_.append( hand_ );
	frame_name_.append( "_D" );
	frame_name_.append( digit_ );
	frame_name_.append( "_M" );
	frame_name_.append( mag_ );
	if ( !camera_suffix_.empty() ) {
		frame_name_.append( "_C" );
		frame_name_.append( camera_suffix_ );
	}
	
	/*
	frame_dir_.append( "/" );
  if (sub_dir_.length() > 0)
  {
	  frame_dir_.append( sub_dir_ );
	  frame_dir_.append( "/" );
  }*/

	vcl_cout << "Frame base name set to: " << frame_dir_ << frame_name_ << vcl_endl;
}