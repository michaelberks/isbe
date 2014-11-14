#include "ncm_qdata_manager.h"

ncm_qdata_manager* ncm_qdata_manager::_instance = 0;

ncm_qdata_manager* ncm_qdata_manager::Instance()
{
	if (_instance == 0 )
		_instance = new ncm_qdata_manager;

  return _instance;
}

//: Default constructor.
ncm_qdata_manager::ncm_qdata_manager()
{
}

//: Return pointer to the main queue (frame buffer).
ncm_video_frame_queue* ncm_qdata_manager::main_queue()
{
	return &main_queue_;
}

//: Add an image to the frame buffer.
void ncm_qdata_manager::add_to_main_queue( QSharedPointer< ncm_video_frame > _qs)
{
	main_queue_.enqueue( _qs );
}

//: Pop an image from the frame buffer.
QSharedPointer< ncm_video_frame > ncm_qdata_manager::remove_from_main_queue( )
{
	return main_queue_.dequeue();
}

//: Return pointer to the main queue (frame buffer).
ncm_video_frame_queue* ncm_qdata_manager::save_queue()
{
	return &save_queue_;
}

//: Add an image to the 'save' queue for writing to disk.
void ncm_qdata_manager::add_to_save_queue( QSharedPointer< ncm_video_frame > _qs)
{
	save_queue_.enqueue( _qs );
}

//: Pop an image from the save queue.
QSharedPointer< ncm_video_frame > ncm_qdata_manager::remove_from_save_queue( )
{
	return save_queue_.dequeue();
}

//: Get the 
QImage* ncm_qdata_manager::get_qimage()
{
	qimage_mutex_.lock();
	QImage *_qimage = &qimage_;
	qimage_mutex_.unlock();
	return _qimage;
}

QImage* ncm_qdata_manager::get_qlastframe()
{
	qlastframe_mutex_.lock();
	QImage *_qimage = &qlastframe_;
	qlastframe_mutex_.unlock();
	return _qimage;
}

void ncm_qdata_manager::set_qimage( QSharedPointer< ncm_video_frame > _vid_frame )
{
	qimage_mutex_.lock();
	//convert_qimage(qimage_, _vid_frame); 
	qimage_mutex_.unlock();
}

void ncm_qdata_manager::set_qlastframe( QSharedPointer< ncm_video_frame > _vid_frame )
{
	qlastframe_mutex_.lock();
	//convert_qimage(qlastframe_, _vid_frame); 
	qlastframe_mutex_.unlock();
}
