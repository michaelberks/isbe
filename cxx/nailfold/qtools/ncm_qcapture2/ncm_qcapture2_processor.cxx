#include "ncm_qcapture2_processor.h"

//  Class that pops an image from the processing queue, performs any
//  necessary processing on the frame, and places it on the save queue for
//  the saver to deal with.
 
//
// Public methods
//

ncm_qcapture2_processor::ncm_qcapture2_processor()
: num_frames_to_save_(2,0),
	n_processed_(2,0),
  sharpness_(-1.0),
  lower_sharpness_threshold_(-1e9),
  upper_sharpness_threshold_(1e9),
	diff_frame_ready_(2, false),
	diff_frame_fliph_(2, false),
	diff_frame_flipv_(2, false)
{
}

//: Set number of frames to save.
void ncm_qcapture2_processor::set_num_frames_to_save(int num_frames_to_save, int cam_num )
{
	num_frames_to_save_[cam_num] = num_frames_to_save;

  // Ensure that is_processing() returns false.
  n_processed_[cam_num] = num_frames_to_save;
}
 
//: Return the number of frames to be saved.
int ncm_qcapture2_processor::num_frames_to_save( int cam_num ) const
{
  return num_frames_to_save_[ cam_num ];
}
 
//: Begin transferring images from the main queue to the save queue.
void ncm_qcapture2_processor::start_processing( int cam_num )
{
	n_processed_[cam_num] = 0;
}

//: Return true if there are still frames to process.
bool ncm_qcapture2_processor::is_processing( int cam_num ) const
{
  return (n_processed_[cam_num] < num_frames_to_save_[cam_num]);
}
 
//: Return the sharpness of the current image.
double ncm_qcapture2_processor::image_sharpness() const
{
  return sharpness_;
}
 
//: Set the lower threshold on sharpness.
void ncm_qcapture2_processor::set_lower_sharpness_threshold(double lower)
{
  lower_sharpness_threshold_ = lower;
}
 
//: Set the lower threshold on sharpness.
void ncm_qcapture2_processor::set_upper_sharpness_threshold(double upper)
{
  upper_sharpness_threshold_ = upper;
}
 
//
// Public slots
//
 
//: Pop a frame from the main queue, process it, and push it onto the save 
//  queue.
void ncm_qcapture2_processor::process_frame(int cam_num , bool diff_frame, bool fliph, bool flipv)
{
	// Pop frame from main queue
	QSharedPointer< ncm_video_frame > vid_frame = 
      ncm_qcapture2_data_manager::Instance()->remove_from_main_queue( cam_num );

  // Check if we're waiting to process more frames.
	if (is_processing( cam_num ))
	{
    // Tag the frame with its sharpness
    // May slow down the frame rate so might revisit this later.
    double sharpness = sharpness_evaluator_.sharpness(*vid_frame->frame());
    vid_frame->set_sharpness(sharpness);

    // Push to the save queue.
		ncm_qcapture2_data_manager::Instance()->add_to_save_queue(cam_num, vid_frame);

    ++n_processed_[cam_num];

		if (cam_num==0)
			emit frame_to_save1(n_processed_[0], 0);
		else
			emit frame_to_save2(n_processed_[1], 1);

    if (n_processed_[cam_num] == num_frames_to_save_[cam_num])
    {
      vid_frame->set_last_frame(true);
      //ncm_qcapture2_data_manager::Instance()->set_qdiffframe( vid_frame );
      stop_processing();
    }
	}

  // Only draw the frame if it was last in the queue.
  const bool that_was_last_frame = 
      ncm_qcapture2_data_manager::Instance()->main_queue( cam_num )->is_empty();

  if (that_was_last_frame)
  {
	  ncm_qcapture2_data_manager::Instance()->set_qimage(cam_num, vid_frame);

	  emit frame_to_draw( cam_num );

    //update_sharpness(cam_num);
  }

	if ( diff_frame ) {
		if (!cam_num) {
			diff_frame1_ = vid_frame;
			diff_frame_ready_[0] = true;
			diff_frame_fliph_[0] = fliph;
			diff_frame_flipv_[0] = flipv;
		}
		else {
			diff_frame2_ = vid_frame;
			diff_frame_ready_[1] = true;
			diff_frame_fliph_[1] = fliph;
			diff_frame_flipv_[1] = flipv;
		}
		if (diff_frame_ready_[0] && diff_frame_ready_[1]) {
			ncm_qcapture2_data_manager::Instance()->set_qdiffframe( diff_frame1_, diff_frame2_, diff_frame_fliph_, diff_frame_flipv_ );
			emit diff_to_draw();
		}
			
	}
			
}

//
// Private methods
//
 
//: Emit a signal to indicate that all frames have been processed.
void ncm_qcapture2_processor::stop_processing()
{
	if ( (n_processed_[0] >= num_frames_to_save_[0]) && (n_processed_[1] >= num_frames_to_save_[1]) ) {
		n_processed_[0] = num_frames_to_save_[0];
		n_processed_[1] = num_frames_to_save_[1];
		emit processing_finished();
	}
}
 
//: Update the sharpness score and check for threshold crossings.
void ncm_qcapture2_processor::update_sharpness( int cam_num )
{
	QImage* qimage = ncm_qcapture2_data_manager::Instance()->get_qimage( cam_num );

  if (!qimage->isNull())
  {
    vil_image_view<vxl_byte> vxl_image;
    qcore_convert_image(vxl_image, *qimage);

    sharpness_ = sharpness_evaluator_.sharpness(vxl_image);

    if (sharpness_ > upper_sharpness_threshold_)
      emit upper_sharpness_exceeded();
    else if (sharpness_ < lower_sharpness_threshold_)
      emit lower_sharpness_exceeded();
  }
  else
  {
    sharpness_ = -1.0;
  }
}
/*
ncm_qcapture2_processor::ncm_qcapture2_processor()
	:		saving1_(false),
			saving2_(false),
			curr_frame1_(0),
			curr_frame2_(0)
{
}

ncm_qcapture2_processor::~ncm_qcapture2_processor()
{

}

void ncm_qcapture2_processor::set_num_frames_to_save1(int num_frames_to_save)
{
	num_frames_to_save1_ = num_frames_to_save;
}
void ncm_qcapture2_processor::set_num_frames_to_save2(int num_frames_to_save)
{
	num_frames_to_save2_ = num_frames_to_save;
}

//: Return the number of frames to be saved.
int ncm_qcapture2_processor::num_frames_to_save1() const
{
  return num_frames_to_save1_;
}
int ncm_qcapture2_processor::num_frames_to_save2() const
{
  return num_frames_to_save2_;
}

//: Return true if there are still frames to process.
bool ncm_qcapture2_processor::is_processing() const
{
  return (n_processed1_ < num_frames_to_save1_);
}

void ncm_qcapture2_processor::start_processing( bool use1, bool use2 )
{
	n_processed1_ = 0;
	saving1_ = use1;
	saving2_ = use2;
}

void ncm_qcapture2_processor::process_frame1()
{
	//Pop frame from main queue
	QSharedPointer< ncm_video_frame > _vid_frame = ncm_qcapture2_data_manager::Instance()->remove_from_main_queue1();

	//Check if we're currently saving frames
	if ( saving1_ )
	{
		if (curr_frame1_ < num_frames_to_save1_)
		{
			//Increment frame count
			curr_frame1_++;

			//Add to the save queue and emit a save ready signal
			ncm_qcapture2_data_manager::Instance()->add_to_save_queue1( _vid_frame );
			//emit frame_to_save(curr_frame1_, 1);
			emit frame_to_save1(curr_frame1_);
		}
		else
		{
			//Stop saving and reset the current frame count
			saving1_ = false;
			curr_frame1_ = 0;
			if ( !saving2_ )
				emit save_finished();
		}
	}

	//Set QImage from this frame and emit frame to draw signal
	ncm_qcapture2_data_manager::Instance()->set_qimage1( _vid_frame );
	emit frame_to_draw1();

}

void ncm_qcapture2_processor::process_frame2()
{
	//Pop frame from main queue
	QSharedPointer< ncm_video_frame > _vid_frame = ncm_qcapture2_data_manager::Instance()->remove_from_main_queue2();

	//Check if we're currently saving frames
	if ( saving2_ )
	{
		if (curr_frame2_ < num_frames_to_save2_)
		{
			//Increment frame count
			curr_frame2_++;

			//Add to the save queue and emit a save ready signal
			ncm_qcapture2_data_manager::Instance()->add_to_save_queue2( _vid_frame );
			//emit frame_to_save(curr_frame2_, 2);
			emit frame_to_save2(curr_frame2_);
		}
		else
		{
			//Stop saving and reset the current frame count
			saving2_ = false;
			curr_frame2_ = 0;
			if ( !saving1_ )
				emit save_finished();
		}
	}

	//Set QImage from this frame and emit frame to draw signal
	ncm_qcapture2_data_manager::Instance()->set_qimage2( _vid_frame );
	emit frame_to_draw2();
	//vcl_cout << "Camera 2 frame ready to draw!!" << vcl_endl;

}*/