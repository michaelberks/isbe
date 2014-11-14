#include "ncm_qcapture_data_manager.h"

ncm_qcapture_data_manager* ncm_qcapture_data_manager::_instance = 0;

ncm_qcapture_data_manager* ncm_qcapture_data_manager::Instance()
{
	if (_instance == 0 )
		_instance = new ncm_qcapture_data_manager;

  return _instance;
}

//: Default constructor.
ncm_qcapture_data_manager::ncm_qcapture_data_manager()
{
}

//: Return pointer to the main queue (frame buffer).
ncm_video_frame_queue* ncm_qcapture_data_manager::main_queue()
{
	return &main_queue_;
}

//: Add an image to the frame buffer.
void ncm_qcapture_data_manager::add_to_main_queue( QSharedPointer< ncm_video_frame > _qs)
{
	main_queue_.enqueue( _qs );
}

//: Pop an image from the frame buffer.
QSharedPointer< ncm_video_frame > ncm_qcapture_data_manager::remove_from_main_queue( )
{
	return main_queue_.dequeue();
}

//: Empty the main queue buffer
void ncm_qcapture_data_manager::empty_main_queue( )
{
	main_queue_.empty();
}

//: Return pointer to the save queue (frame buffer).
ncm_video_frame_queue* ncm_qcapture_data_manager::save_queue()
{
	return &save_queue_;
}

//: Add an image to the 'save' queue for writing to disk.
void ncm_qcapture_data_manager::add_to_save_queue( QSharedPointer< ncm_video_frame > _qs)
{
	save_queue_.enqueue( _qs );
}

//: Pop an image from the save queue.
QSharedPointer< ncm_video_frame > ncm_qcapture_data_manager::remove_from_save_queue( )
{
	return save_queue_.dequeue();
}

//: Return pointer to the align queue (frame buffer).
ncm_video_frame_queue* ncm_qcapture_data_manager::align_queue()
{
	return &align_queue_;
}

//: Add an image to the 'align' queue for writing to disk.
void ncm_qcapture_data_manager::add_to_align_queue( QSharedPointer< ncm_video_frame > _qs)
{
	align_queue_.enqueue( _qs );
}

//: Pop an image from the align queue.
QSharedPointer< ncm_video_frame > ncm_qcapture_data_manager::remove_from_align_queue( )
{
	return align_queue_.dequeue();
}

//: Return pointer to the mosaic queue (frame buffer).
ncm_video_frame_queue* ncm_qcapture_data_manager::mosaic_queue()
{
	return &mosaic_queue_;
}

//: Add an image to the 'mosaic' queue for writing to disk.
void ncm_qcapture_data_manager::add_to_mosaic_queue( QSharedPointer< ncm_video_frame > _qs)
{
	mosaic_queue_.enqueue( _qs );
}

//: Pop an image from the mosaic queue.
QSharedPointer< ncm_video_frame > ncm_qcapture_data_manager::remove_from_mosaic_queue( )
{
	return mosaic_queue_.dequeue();
}

//: Get the 
QImage* ncm_qcapture_data_manager::get_qimage()
{
	qimage_mutex_.lock();
	QImage *_qimage = &qimage_;
	qimage_mutex_.unlock();
	return _qimage;
}

QImage* ncm_qcapture_data_manager::get_qlastframe()
{
	qlastframe_mutex_.lock();
	QImage *_qimage = &qlastframe_;
	qlastframe_mutex_.unlock();
	return _qimage;
}

void ncm_qcapture_data_manager::set_qimage( QSharedPointer< ncm_video_frame > _vid_frame )
{
	qimage_mutex_.lock();
	convert_qimage(qimage_, _vid_frame); 
	qimage_mutex_.unlock();
}

void ncm_qcapture_data_manager::set_qlastframe( QSharedPointer< ncm_video_frame > _vid_frame )
{
	qlastframe_mutex_.lock();
	convert_qimage(qlastframe_, _vid_frame); 
	qlastframe_mutex_.unlock();
}

//
// Private methods
//

void ncm_qcapture_data_manager::convert_qimage(QImage &qimage, QSharedPointer< ncm_video_frame > _vid_frame )
{
	int src_cols = _vid_frame->frame()->ni();
	int src_rows = _vid_frame->frame()->nj();

	if ( qimage.isNull() || 
      (qimage.width() != src_cols) || 
      (qimage.height() != src_rows) ) 
	{
		
		// Initialise qimage
		qimage = QImage( src_cols, src_rows, QImage::Format_Indexed8 );
		
		// Set colour values
		const unsigned n_colours = 256;
		qimage.setNumColors( n_colours );

		for (unsigned i = 0; i < n_colours; ++i)
		{
			qimage.setColor(i,qRgb(i,i,i));
		}

		vcl_cout << "Initialised qimage: " 
             << qimage.width() << " x " << qimage.height() << vcl_endl;
	}

	// Copy pixel data from frame to qimage
	vcl_memcpy( qimage.bits(), 
              _vid_frame->frame()->top_left_ptr(), 
              src_cols*src_rows );
}