#include "ncm_qcapture_processor.h"
#include "ncm_qcapture_data_manager.h"

#include <QGraphicsPixmapItem>
#include <QGraphicsScene>
#include <QVector>

#include <vcl_iostream.h>
#include <vcl_iomanip.h>
#include <vcl_string.h>
#include <vcl_sstream.h>
#include <vcl_iosfwd.h>
#include <vil/vil_save.h>

#include <qcore/qcore_convert_image.h>
 
 
//  Class that pops an image from the processing queue, performs any
//  necessary processing on the frame, and places it on the save queue for
//  the saver to deal with.
 
//
// Public methods
//
 
ncm_qcapture_processor::ncm_qcapture_processor(ncm_qcapture_preferences & preferences)
: preferences_(preferences),
	num_frames_to_save_(0),
	n_processed_(0),
  sharpness_array_size_(preferences_.n_sharpness_),
  sharpness_array_(preferences_.n_sharpness_, 0.0),
  sharpness_array_index_(0),
  lower_sharpness_threshold_(-1e9),
  upper_sharpness_threshold_(1e9),
	is_paused_(false),
	is_stopped_(false)
{
}
 
//: Set number of frames to save.
void ncm_qcapture_processor::set_num_frames_to_save(int num_frames_to_save)
{
	if (num_frames_to_save < 0)
		num_frames_to_save_ = preferences_.max_frames_to_save_;
	else
		num_frames_to_save_ = num_frames_to_save;

  // Ensure that is_processing() returns false.
  n_processed_ = num_frames_to_save;
}
 
//: Return the number of frames to be saved.
int ncm_qcapture_processor::num_frames_to_save() const
{
  return num_frames_to_save_;
}
 
//: Begin transferring images from the main queue to the save queue.
void ncm_qcapture_processor::start_processing()
{
	n_processed_ = 0;
}

//: Stop saving - we need this now we save until we're told to stop
//rather than knowing in advance how many frames to process
void ncm_qcapture_processor::stop_saving()
{
	is_stopped_ = true;
}

//: Return true if there are still frames to process.
bool ncm_qcapture_processor::is_processing() const
{
  return (n_processed_ < num_frames_to_save_);
}

//: Return true if saving is paused.
bool ncm_qcapture_processor::is_paused() const
{
  return is_paused_;
}

//: Return true if saving is paused.
bool ncm_qcapture_processor::is_stopped() const
{
  return is_stopped_;
}
 
//: Return the average sharpness over the ring buffer.
double ncm_qcapture_processor::image_sharpness() const
{
  double sharpness_sum = 0.0;

  for (unsigned i = 0; i < sharpness_array_size_; ++i)
    sharpness_sum += sharpness_array_[i];
  
  return sharpness_sum / sharpness_array_size_;
}
 
//: Set paused status
void ncm_qcapture_processor::set_paused(bool paused)
{
  is_paused_ = paused;
}

//: Set the lower threshold on sharpness.
void ncm_qcapture_processor::set_lower_sharpness_threshold(double lower)
{
  lower_sharpness_threshold_ = lower;
}
 
//: Return the lower threshold on sharpness.
double ncm_qcapture_processor::lower_sharpness_threshold() const
{
  return lower_sharpness_threshold_;
}
 
//: Set the upper threshold on sharpness.
void ncm_qcapture_processor::set_upper_sharpness_threshold(double upper)
{
  upper_sharpness_threshold_ = upper;
}
 
//: Return the upper threshold on sharpness.
double ncm_qcapture_processor::upper_sharpness_threshold() const
{
  return upper_sharpness_threshold_;
}
 
//
// Public slots
//
 
//: Pop a frame from the main queue, process it, and push it onto the save 
//  queue.
void ncm_qcapture_processor::process_frame()
{
	// Pop frame from main queue
	QSharedPointer< ncm_video_frame > vid_frame = 
      ncm_qcapture_data_manager::Instance()->remove_from_main_queue();

  // Check if we're waiting to save more frames.
	if (is_processing())
	{
		bool align_frame;

		if (!is_paused()) //When we're paused we don't save frames, but still need to check for being stopped
		{
			++n_processed_;

			if (preferences_.compute_sharpness_)
			{
				// Tag the frame with its sharpness
				// May slow down the frame rate so might revisit this later.
				double sharpness = sharpness_evaluator_.sharpness(*vid_frame->frame());
				vid_frame->set_sharpness(sharpness);		

				// Add sharpness to the ring buffer.
				store_sharpness(sharpness);
			}

			//Tag the frame with a number and name
			vid_frame->set_frame_num(n_processed_);
			vcl_stringstream ss;
			ss    << preferences_.frame_prefix_
						<< vcl_setw( 5 ) << vcl_setfill( '0' ) << (n_processed_) << ".png";
			vid_frame->set_frame_name(ss.str());

			// Push to the save queue.
			ncm_qcapture_data_manager::Instance()->add_to_save_queue(vid_frame);
			

			//Check if we need to align this frame
			align_frame = (preferences_.align_nth_frame_>0) && 
				(n_processed_ % preferences_.align_nth_frame_ == 0);
			if (align_frame)
			{
				ncm_qcapture_data_manager::Instance()->add_to_align_queue(vid_frame);
				emit frame_to_align(n_processed_);
			}
			emit frame_to_save(vid_frame->get_frame_header());
		}

    if (is_stopped_ || n_processed_ == num_frames_to_save_)
    {
			//Tag this as the last frame
      vid_frame->set_last_frame(true);
      ncm_qcapture_data_manager::Instance()->set_qlastframe( vid_frame );

			if ((preferences_.align_nth_frame_>0) && !align_frame)
			{
				//If we're aligning, always align the last frame
				ncm_qcapture_data_manager::Instance()->add_to_align_queue(vid_frame);
				emit frame_to_align(n_processed_);
			}

			//We're done
      stop_processing();
    }
	}

  // Only draw the frame if it was last in the queue.
  const bool that_was_last_frame = 
      ncm_qcapture_data_manager::Instance()->main_queue()->is_empty();

  if (that_was_last_frame)
  {
	  ncm_qcapture_data_manager::Instance()->set_qimage(vid_frame);

	  emit frame_to_draw();

		if (preferences_.compute_sharpness_)
		{
			double sharpness = sharpness_evaluator_.sharpness(*vid_frame->frame());
			store_sharpness(sharpness);
			update_sharpness();
		}
  }
}

//
// Private methods
//
 
//: Emit a signal to indicate that all frames have been processed.
void ncm_qcapture_processor::stop_processing()
{
	int total_processed = n_processed_;
  num_frames_to_save_ = n_processed_;
	is_paused_ = false;
	is_stopped_= false;
  emit processing_finished(total_processed);
}
 
void ncm_qcapture_processor::resize_sharpness_array()
{
	if (preferences_.n_sharpness_ == sharpness_array_size_)
		return;

	else
		resize_sharpness_array(preferences_.n_sharpness_);
}

//: Resize the sharpness array.
void ncm_qcapture_processor::resize_sharpness_array(int new_size)
{
  // Resize only if 0 < new_size <= 100.
  if ( (0 < new_size) && 
           (new_size <= 100) )
  {
    // Record the current average.
    const double current_average = image_sharpness();

    // Resize the array and fill with copies of the average.
    sharpness_array_size_ = new_size;
    sharpness_array_.resize(new_size, current_average);
    sharpness_array_index_ = 0;
  }
}

//: Add a sharpness value to the ring buffer.
void ncm_qcapture_processor::store_sharpness(double sharpness)
{
  sharpness_array_[sharpness_array_index_] = sharpness;
  
  ++sharpness_array_index_;

  if (sharpness_array_index_ == sharpness_array_size_)
    sharpness_array_index_ = 0;
}

//: Update the sharpness score and check for threshold crossings.
void ncm_qcapture_processor::update_sharpness()
{
	//QImage* qimage = ncm_qcapture_data_manager::Instance()->get_qimage();

  //if (!qimage->isNull())
  //{
    //vil_image_view<vxl_byte> vxl_image;
    //qcore_convert_image(vxl_image, *qimage);

    const double sharpness = image_sharpness();

    if (sharpness > upper_sharpness_threshold_)
      emit upper_sharpness_exceeded();
    else if (sharpness < lower_sharpness_threshold_)
      emit lower_sharpness_exceeded();
  //}
  //else
  //{
  //  sharpness_ = -1.0;
  //}
}