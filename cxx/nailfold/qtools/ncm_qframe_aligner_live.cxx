#include "ncm_qframe_aligner_live.h"
 
#include <vil/vil_load.h>
#include <nailfold/qtools/ncm_video_frame.h>
#include <nailfold/qtools/ncm_qcapture/ncm_qcapture_data_manager.h>
 
//
// Public methods
//
 
/*ncm_qframe_aligner_live::ncm_qframe_aligner_live()
	: MAX_LENGTH(1e4)
{
}*/

ncm_qframe_aligner_live::ncm_qframe_aligner_live(int max_length /*= 1e4*/)
	: MAX_LENGTH(max_length)
{
}
 
//: Set the number of frames that should be aligned.
void ncm_qframe_aligner_live::setNumToAlign(
  unsigned n_to_align)
{
  n_to_align_ = MAX_LENGTH;
}
 
//: Return to initial state where destination is first image.
void ncm_qframe_aligner_live::reset()
{
	// Pop frame from the align queue
	QSharedPointer< ncm_video_frame > vid_frame = 
    ncm_qcapture_data_manager::Instance()->remove_from_align_queue();

	vil_image_view<vxl_byte> img = vid_frame->frame();
	ncm_qframe_aligner::reset(img);

}
 
//
// Public slots
//
 
//: Align the current frame with the previous one.
void ncm_qframe_aligner_live::alignNextFrame()
{
  if (n_aligned_ < n_to_align_)
  {
		// Pop frame from the align queue
		QSharedPointer< ncm_video_frame > vid_frame = 
      ncm_qcapture_data_manager::Instance()->remove_from_align_queue();

    vil_image_view<vxl_byte> img = vid_frame->frame();
		ncm_qframe_aligner::alignNextFrame(img, vid_frame->is_last_frame());
	}
}
 
//
// Private methods
//
