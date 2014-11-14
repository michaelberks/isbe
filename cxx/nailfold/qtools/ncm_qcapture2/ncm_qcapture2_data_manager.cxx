#include "ncm_qcapture2_data_manager.h"

ncm_qcapture2_data_manager* ncm_qcapture2_data_manager::_instance = 0;

ncm_qcapture2_data_manager* ncm_qcapture2_data_manager::Instance()
{
	if (_instance == 0 )
	{
		_instance = new ncm_qcapture2_data_manager;
	}
	return _instance;
}

ncm_qcapture2_data_manager::ncm_qcapture2_data_manager()
{
}
//-------------------------------------------------------------
ncm_video_frame_queue* ncm_qcapture2_data_manager::main_queue(int cam_num)
{
	if ( cam_num )
		return &main_queue2_;
	else
		return &main_queue1_;
}
/*
ncm_video_frame_queue* ncm_qcapture2_data_manager::main_queue1()
{
	return &main_queue1_;
}
ncm_video_frame_queue* ncm_qcapture2_data_manager::main_queue2()
{
	return &main_queue2_;
}*/

//-------------------------------------------------------------
ncm_video_frame_queue* ncm_qcapture2_data_manager::save_queue(int cam_num)
{
	if ( cam_num )
		return &save_queue2_;
	else
		return &save_queue1_;
}
/*
ncm_video_frame_queue* ncm_qcapture2_data_manager::save_queue1()
{
	return &save_queue1_;
}
ncm_video_frame_queue* ncm_qcapture2_data_manager::save_queue2()
{
	return &save_queue2_;
}
*/
//----------------------------------------------------------------
void ncm_qcapture2_data_manager::add_to_main_queue( int cam_num, QSharedPointer< ncm_video_frame > _qs)
{
	if ( cam_num )
		main_queue2_.enqueue( _qs );
	else
		main_queue1_.enqueue( _qs );
}
/*
void ncm_qcapture2_data_manager::add_to_main_queue1( QSharedPointer< ncm_video_frame > _qs)
{
	main_queue1_.enqueue( _qs );
}
void ncm_qcapture2_data_manager::add_to_main_queue2( QSharedPointer< ncm_video_frame > _qs)
{
	main_queue2_.enqueue( _qs );
}
*/
//-----------------------------------------------------------------
void ncm_qcapture2_data_manager::add_to_save_queue( int cam_num, QSharedPointer< ncm_video_frame > _qs)
{
	if ( cam_num )
		save_queue2_.enqueue( _qs );
	else
		save_queue1_.enqueue( _qs );
}
/*
void ncm_qcapture2_data_manager::add_to_save_queue1( QSharedPointer< ncm_video_frame > _qs)
{
	save_queue1_.enqueue( _qs );
}
void ncm_qcapture2_data_manager::add_to_save_queue2( QSharedPointer< ncm_video_frame > _qs)
{
	save_queue2_.enqueue( _qs );
}
*/
//-------------------------------------------------------------------
QSharedPointer< ncm_video_frame > ncm_qcapture2_data_manager::remove_from_main_queue( int cam_num )
{
	if ( cam_num )
		return main_queue2_.dequeue( );
	else
		return main_queue1_.dequeue( );
}
/*
QSharedPointer< ncm_video_frame > ncm_qcapture2_data_manager::remove_from_main_queue1( )
{
	return main_queue1_.dequeue();
}
QSharedPointer< ncm_video_frame > ncm_qcapture2_data_manager::remove_from_main_queue2( )
{
	return main_queue2_.dequeue();
}
*/
//--------------------------------------------------------------------
QSharedPointer< ncm_video_frame > ncm_qcapture2_data_manager::remove_from_save_queue( int cam_num )
{
	if ( cam_num )
		return save_queue2_.dequeue( );
	else
		return save_queue1_.dequeue( );
}
/*
QSharedPointer< ncm_video_frame > ncm_qcapture2_data_manager::remove_from_save_queue1( )
{
	return save_queue1_.dequeue();
}
QSharedPointer< ncm_video_frame > ncm_qcapture2_data_manager::remove_from_save_queue2( )
{
	return save_queue2_.dequeue();
}
*/
//-------------------------------------------------------------------------
QImage* ncm_qcapture2_data_manager::get_qimage( int cam_num )
{
	if ( cam_num )//==1
	{
		qimage_mutex2_.lock();
		QImage *_qimage = &qimage2_;
		qimage_mutex2_.unlock();
		return _qimage;
	}
	else {
		qimage_mutex1_.lock();
		QImage *_qimage = &qimage1_;
		qimage_mutex1_.unlock();
		return _qimage;
	}	
}
/*
QImage* ncm_qcapture2_data_manager::get_qimage1()
{
	qimage_mutex1_.lock();
	QImage *_qimage = &qimage1_;
	qimage_mutex1_.unlock();
	return _qimage;
}

QImage* ncm_qcapture2_data_manager::get_qimage2()
{
	qimage_mutex2_.lock();
	QImage *_qimage = &qimage2_;
	qimage_mutex2_.unlock();
	return _qimage;
}
*/
//------------------------------------------------------------
void ncm_qcapture2_data_manager::set_qimage(int cam_num, QSharedPointer< ncm_video_frame > _vid_frame )
{
	if ( cam_num ) //==1
		set_qimage( qimage2_, qimage_mutex2_, _vid_frame );
	else
		set_qimage( qimage1_, qimage_mutex1_, _vid_frame );
}
/*
void ncm_qcapture2_data_manager::set_qimage1(QSharedPointer< ncm_video_frame > _vid_frame )
{
	set_qimage( qimage1_, qimage_mutex1_, _vid_frame );
}

void ncm_qcapture2_data_manager::set_qimage2(QSharedPointer< ncm_video_frame > _vid_frame )
{
	set_qimage( qimage2_, qimage_mutex2_, _vid_frame );
}
*/
//---------------------------------------------------------------
QImage* ncm_qcapture2_data_manager::get_qdiffframe()
{
	qdiffframe_mutex_.lock();
	QImage *_qimage = &qdiffframe_;
	qdiffframe_mutex_.unlock();
	return _qimage;
}

void ncm_qcapture2_data_manager::set_qdiffframe( QSharedPointer< ncm_video_frame > _vid_frame1,  QSharedPointer< ncm_video_frame > _vid_frame2, QVector<bool> fliph, QVector<bool> flipv)
{
	qdiffframe_mutex_.lock();
	diff_qimage(qdiffframe_, _vid_frame1, _vid_frame2, fliph, flipv); 
	qdiffframe_mutex_.unlock();
}

//
// Private methods:
//
void ncm_qcapture2_data_manager::set_qimage( QImage &qimage, QMutex &qimage_mutex, QSharedPointer< ncm_video_frame > _vid_frame )
{
	int src_cols = _vid_frame->frame()->ni();
	int src_rows = _vid_frame->frame()->nj();

	//Check if we need to initialise the Qimage, otherwise we can just copy across
	//values assuming width, height, format etc are unchanged
	qimage_mutex.lock();
	if ( qimage.isNull() || (qimage.width() != src_cols) || (qimage.height() != src_rows) ) 
	{
		
		//Initialise qimage
		qimage = QImage( src_cols, src_rows, QImage::Format_Indexed8 );
		
		// set colour values
		unsigned int n_colours= 256;
		qimage.setNumColors( n_colours );
		for (unsigned int i=0;i<n_colours;++i)
		{
			qimage.setColor(i,qRgb(i,i,i));
		}
		vcl_cout << "Initialised qimage: " << qimage.width() << " x " << qimage.height() << vcl_endl;
	}
	
	//Copy pixel data from frame to qimage
	vcl_memcpy( qimage.bits(), _vid_frame->frame()->top_left_ptr(), src_cols*src_rows );
	qimage_mutex.unlock();

}

void ncm_qcapture2_data_manager::convert_qimage(QImage &qimage, QSharedPointer< ncm_video_frame > _vid_frame )
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

void ncm_qcapture2_data_manager::diff_qimage(QImage &qimage, QSharedPointer< ncm_video_frame > _vid_frame1,  QSharedPointer< ncm_video_frame > _vid_frame2, QVector<bool> fliph, QVector<bool> flipv)
{
	//convert_qimage( qimage, _vid_frame2);
	int src_cols = _vid_frame1->frame()->ni();
	int src_rows = _vid_frame1->frame()->nj();

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
	
	// Loop through each pixel, subtracting vid_frame2 from vid_frame1 and storing in the qimage
	vxl_byte *_q_im = qimage.bits();
	vxl_byte *_vf1 = _vid_frame1->frame()->top_left_ptr();
	vxl_byte *_vf2 = _vid_frame2->frame()->top_left_ptr();

	int i_q = 0;
	int i_1 = 0;
	int i_2 = 0;
	
	for (int i_r = 0; i_r < src_rows; ++i_r) {

		if (flipv[0]) //vid_frame 1 v flipped, set i_1 from bottom row
			i_1 = (src_rows - i_r) * src_cols;
		else
			i_1 = i_q;

		if (flipv[1]) //vid_frame 2 flipped, set i_2 from bottom row
			i_2 = (src_rows - i_r) * src_cols;
		else
			i_2 = i_q;

		if (fliph[0]) //vid_frame 1 flipped, set i_1 to end of row
			i_1+=(src_cols - 1);	

		if (fliph[1]) //vid_frame 2 flipped, set i_2 to end of row
			i_2+=(src_cols - 1);

		for (int i_c = 0; i_c < src_cols; ++i_c) {
			//Take difference of frames and map to 0-255
			_q_im[i_q] = ( 128 - (_vf2[i_2]>>1) ) + (_vf1[i_1]>>1);
			
			//Update iterators
			i_q++;

			if (fliph[0])
				i_1--; //Move backward along the row
			else
				i_1++; //Move forwards

			if (fliph[1])
				i_2--; //Move backward along the row
			else
				i_2++; //Move forwards

		}
	}

	/*
	for (int i = 0; i < src_cols*src_rows; ++i) {
		_q_im[i] = ( 128 - (_vf2[i]>>1) ) + (_vf1[i]>>1);
	}
	*/
}