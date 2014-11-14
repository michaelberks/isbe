#include "ncm_qdata_processor.h"

#include <QGraphicsPixmapItem>
#include <QGraphicsScene>

#include <vcl_iostream.h>
#include <vcl_iomanip.h>
#include <vcl_string.h>
#include <vcl_sstream.h>
#include <vcl_iosfwd.h>
#include <vil/vil_save.h>

#include <qcore/qcore_convert_image.h>

#include <nailfold/qtools/ncm_video_frame.h>
#include <nailfold/qtools/convert_ncm_video_frame_to_qimage.h>
 
 
//  Class that pops an image from the processing queue, performs any
//  necessary processing on the frame, and places it on the save queue for
//  the saver to deal with.
 
//
// Public methods
//
 
ncm_qdata_processor::ncm_qdata_processor()
: input_queue_(NULL),
  output_queue_(NULL),
  num_frames_to_save_(0),
	n_processed_(0),
  sharpness_(-1.0),
  lower_sharpness_threshold_(-1e9),
  upper_sharpness_threshold_(1e9)
{
}
 
//: Set number of frames to save.
void ncm_qdata_processor::set_num_frames_to_save(int num_frames_to_save)
{
	num_frames_to_save_ = num_frames_to_save;

  // Ensure that is_processing() returns false.
  n_processed_ = num_frames_to_save;
}
 
//: Return the number of frames to be saved.
int ncm_qdata_processor::num_frames_to_save() const
{
  return num_frames_to_save_;
}
 
//: Begin transferring images from the main queue to the save queue.
void ncm_qdata_processor::start_processing()
{
	n_processed_ = 0;
}
 
//: Return true if there are still frames to process.
bool ncm_qdata_processor::is_processing() const
{
  return (n_processed_ < num_frames_to_save_);
}
 
//: Return the sharpness of the current image.
double ncm_qdata_processor::image_sharpness() const
{
  return sharpness_;
}
 
//: Set the lower threshold on sharpness.
void ncm_qdata_processor::set_lower_sharpness_threshold(double lower)
{
  lower_sharpness_threshold_ = lower;
}
 
//: Set the lower threshold on sharpness.
void ncm_qdata_processor::set_upper_sharpness_threshold(double upper)
{
  upper_sharpness_threshold_ = upper;
}
 
//: Define the input queue.
void ncm_qdata_processor::set_input_queue(ncm_video_frame_queue* input_queue)
{
  input_queue_ = input_queue;
}
 
//: Define the output queue.
void ncm_qdata_processor::set_output_queue(ncm_video_frame_queue* output_queue)
{
  output_queue_ = output_queue;
}
 
//: Return the frame to display.
QImage* ncm_qdata_processor::qimage()
{
	return &qimage_;
}
 
//: Return the last frame in the sequence.
QImage* ncm_qdata_processor::qlastframe()
{
	return &qlastframe_;
}
 
//
// Public slots
//
 
//: Pop a frame from the main queue, process it, and push it onto the save 
//  queue.
void ncm_qdata_processor::process_frame()
{
  if (input_queue_->is_empty())
    return;

	// Pop frame from main queue
  QSharedPointer<ncm_video_frame> vid_frame = input_queue_->dequeue();

  // Check if we're waiting to process more frames.
	if (is_processing())
	{
    // Tag the frame with its sharpness
    // May slow down the frame rate so might revisit this later.
    double sharpness = sharpness_evaluator_.sharpness(*vid_frame->frame());
    vid_frame->set_sharpness(sharpness);

    // Push to the save queue.
    output_queue_->enqueue(vid_frame);

    ++n_processed_;

    emit frame_to_save(n_processed_);

    if (n_processed_ == num_frames_to_save_)
    {
      vid_frame->set_last_frame(true);
      set_qlastframe(vid_frame);
      stop_processing();
    }
	}

  // Only draw the frame if it was last in the queue.
  const bool that_was_last_frame = input_queue_->is_empty();

  if (that_was_last_frame)
  {
	  set_qimage(vid_frame);
	  emit frame_to_draw();
    update_sharpness_of(*vid_frame);
  }
}

//
// Private methods
//
 
//: Emit a signal to indicate that all frames have been processed.
void ncm_qdata_processor::stop_processing()
{
  n_processed_ = num_frames_to_save_;

  emit processing_finished();
}
 
//: Update the sharpness score and check for threshold crossings.
void ncm_qdata_processor::update_sharpness_of(
  const ncm_video_frame& video_frame)
{
  if (video_frame.frame() != NULL)
  {    
    sharpness_ = sharpness_evaluator_.sharpness(*video_frame.frame());

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
 
//: Make a copy of the frame to display.
void ncm_qdata_processor::set_qimage(
  QSharedPointer< ncm_video_frame > vid_frame)
{
	qimage_mutex_.lock();
	convert_ncm_video_frame_to_qimage(qimage_, vid_frame);
	qimage_mutex_.unlock();
}
 
//: Make a copy of the last frame in the sequence.
void ncm_qdata_processor::set_qlastframe( 
  QSharedPointer< ncm_video_frame > vid_frame)
{
	qlastframe_mutex_.lock();
	convert_ncm_video_frame_to_qimage(qlastframe_, vid_frame);
	qlastframe_mutex_.unlock();
}
