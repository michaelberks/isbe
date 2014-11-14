#include "ncm_qcapture_aligner.h"
 
#include <nailfold/qtools/ncm_video_frame.h>
#include <nailfold/qtools/ncm_qcapture/ncm_qcapture_data_manager.h>
 
//
// Public methods
//
 
/*ncm_qcapture_aligner::ncm_qcapture_aligner()
	: MAX_LENGTH(1e4)
{
}*/

ncm_qcapture_aligner::ncm_qcapture_aligner(int max_length /*= 1e4*/)
	: pixels_per_mm_(1160),
		MAX_LENGTH(max_length),
		use_motor_positions_(true),
		motor_x0_(0),
		motor_y0_(0)
{
}

//: Set the pixel resolution of the motors (i.e. how many pixels does 1mm of motor movement
// equate too - this may be slightly different to 'real' world mm because of inaccuracies in the motors
// but can be calibrated in the capture software)
void ncm_qcapture_aligner::setPixelsPerMm(double pmm)
{
	pixels_per_mm_ = pmm;
}

//: Set the number of frames that should be aligned.
void ncm_qcapture_aligner::setNumToAlign(
  unsigned n_to_align)
{
	if (n_to_align < MAX_LENGTH)
		n_to_align_ = n_to_align;
	else
		n_to_align_ = MAX_LENGTH;

	resetDisplacements();
}
 
//: Return to initial state where destination is first image.
void ncm_qcapture_aligner::reset()
{
	n_aligned_ = 0;
	alignNextFrame();
}
 
//
// Public slots
//
 
//: Align the current frame with the previous one.
void ncm_qcapture_aligner::alignNextFrame()
{
	if (n_aligned_ < n_to_align_)
  {
		/*// Pop frame from the align queue and add it to the mosaic queue
		QSharedPointer< ncm_video_frame > vid_frame = 
      ncm_qcapture_data_manager::Instance()->remove_from_align_queue();
		ncm_qcapture_data_manager::Instance()->add_to_mosaic_queue( vid_frame );*/

		// Get frame from the align queue and add it to the mosaic queue, frame
		// will later be removed from the align queue after it is aligned
		QSharedPointer< ncm_video_frame > vid_frame = 
			ncm_qcapture_data_manager::Instance()->align_queue()->head();
		ncm_qcapture_data_manager::Instance()->add_to_mosaic_queue( vid_frame );

		//Take copy of image data in frame
		vil_image_view<vxl_byte> img = *(vid_frame->frame());
		
		//If first frame set up the aligner
		if (n_aligned_ == 0)
		{
			motor_x0_ = -vid_frame->get_frame_header().motor_position_[0];
			motor_y0_ = vid_frame->get_frame_header().motor_position_[1];
			ncm_qframe_aligner::reset(img);
		}
		else //Otherwise align this frame
		{
			bool is_last_frame = vid_frame->is_last_frame();
			is_last_frame = is_last_frame || (n_aligned_+1 == n_to_align_);

			if (use_motor_positions_)
			{
				//Get motor positions for this frame
				double mi = -vid_frame->get_frame_header().motor_position_[0];
				double mj = vid_frame->get_frame_header().motor_position_[1];

				//Compute difference to previous motor positions and convert to pixels
				//to compute displacements
				int di = int((mi-motor_x0_)*pixels_per_mm_);
				int dj = int((mj-motor_y0_)*pixels_per_mm_);

				setDisplacement(n_aligned_, QPoint(di, dj));

				//Save motor positions for next iteration
				motor_x0_ = mi;
				motor_y0_ = mj;

				//vcl_cout << "Displacement set to x = " << di << ", y = " << dj << vcl_endl;  
			}
			ncm_qframe_aligner::alignNextFrame(img, is_last_frame);
		}
	}
}
 
//
// Private methods
//
