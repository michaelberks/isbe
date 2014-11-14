#include "ncm_dmk_grab_frame_filter.h"

using namespace DShowLib;

ncm_dmk_grab_frame_filter::ncm_dmk_grab_frame_filter()
{
	//Register mem_buffer class so it can be queued by QT
	qRegisterMetaType< QSharedPointer< vil_image_view<vxl_byte> > >("QSharedPointer< vil_image_view<vxl_byte> >");
}

/*
 *	This method returns filter information.
 *	This is not really needed here, because the filter will never be loaded from a .ftf file.
 *
 *	However, the presence of getStaticFilterInfo is required because
 *	we derived from FrameFilterImpl.
 */
FilterInfo ncm_dmk_grab_frame_filter::getStaticFilterInfo()
{
	// Return a filter name and declare the filter as eFC_INTERNAL.
	FilterInfo fi = { L"Binarization", L"", eFC_INTERNAL };
	return fi;
}

/*
 *	This method fills the FrameTypeInfoArray arr with the frame types this filter
 *	accepts as input.
 *
 *	For the binarization filter, only the gray color formats eY800 and eRGB8 are accepted.
 */
void ncm_dmk_grab_frame_filter::getSupportedInputTypes( DShowLib::FrameTypeInfoArray& arr ) const
{
	// This filter works for 8-bit-gray images only
	arr.push_back( eRGB8 );
	arr.push_back( eY800 );
}

/*
 *	This method returns the output frame type for a given input frame type.
 *
 *	The binarization filter does not change size or color format,
 *	so the only output frame type is the input frame type.
 */
bool ncm_dmk_grab_frame_filter::getTransformOutputTypes( const DShowLib::FrameTypeInfo& in_type,
												   DShowLib::FrameTypeInfoArray& out_types ) const
{
	// We don't change the image type, output = input
	out_types.push_back( in_type );

	return true;
}

/*
 *	This method is called to copy image data from the src frame to the dest 
 *  frame.
 *
 *	In our case, we make our own copy of the data and leave the dest frame 
 *  blank, interrupting the default DMK stream of data to the ring buffer sink
 */
bool ncm_dmk_grab_frame_filter::transform( const DShowLib::IFrame& src, 
                                           DShowLib::IFrame& dest )
{
	//Get current time
	QTime curr_time = QTime::currentTime();
	
	//Copy src frame into new video frame object
	BYTE* dmk_im = src.getPtr();
	SIZE dims = src.getFrameType().dim;

	ncm_video_frame* frame = 
      new ncm_video_frame(unsigned (dims.cx), unsigned (dims.cy), curr_time);

	vcl_memcpy(frame->frame()->top_left_ptr(), dmk_im, 
             src.getFrameType().buffersize );

	//Create a shared pointer to the new video frame object
	QSharedPointer< ncm_video_frame > qs_frame = 
      QSharedPointer< ncm_video_frame >( frame );

	//Add the shared pointer to the queue
	queue_->enqueue( qs_frame );

	//Emit frame ready signal
	emit frame_ready();

	return true;
}

void ncm_dmk_grab_frame_filter::set_queue( ncm_video_frame_queue* _queue )
{
	queue_ = _queue;
}

