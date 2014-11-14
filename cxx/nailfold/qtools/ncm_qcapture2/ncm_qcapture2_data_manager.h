#ifndef NCM_QCAPTURE2_DATA_MANAGER_H
#define NCM_QCAPTURE2_DATA_MANAGER_H

#include <QtCore>
#include <QImage>

#include <vcl_cstring.h>
#include <vcl_sstream.h>
#include <vcl_iostream.h>
#include <vil/vil_image_view.h>

#include <nailfold/qtools/ncm_video_frame.h>
#include <nailfold/qtools/ncm_video_frame_queue.h>

class ncm_qcapture2_data_manager
{

	public:
		//Controls access to unique instance for all clients
		static ncm_qcapture2_data_manager* Instance();

		//Rest of interface

		//:Return pointers to main/save queue
		ncm_video_frame_queue* main_queue( int cam_num );
		ncm_video_frame_queue* save_queue( int cam_num );

		/*ncm_video_frame_queue* main_queue1();
		ncm_video_frame_queue* main_queue2();
		ncm_video_frame_queue* save_queue1();
		ncm_video_frame_queue* save_queue2();*/

		//:Add a video frame to main queue
		void add_to_main_queue( int cam_num, QSharedPointer< ncm_video_frame > );
		/*void add_to_main_queue1( QSharedPointer< ncm_video_frame > );
		void add_to_main_queue2( QSharedPointer< ncm_video_frame > );*/
		
		//: Pop a video frame from the main queue
		QSharedPointer< ncm_video_frame > remove_from_main_queue( int cam_num );
		/*QSharedPointer< ncm_video_frame > remove_from_main_queue1( );
		QSharedPointer< ncm_video_frame > remove_from_main_queue2( );*/

		//: Add a video frame to the save queue
		void add_to_save_queue( int cam_num, QSharedPointer< ncm_video_frame > );
		/*void add_to_save_queue1( QSharedPointer< ncm_video_frame > );
		void add_to_save_queue2( QSharedPointer< ncm_video_frame > );*/

		//: Pop a video frame from the main queue
		QSharedPointer< ncm_video_frame > remove_from_save_queue( int cam_num );
		/*QSharedPointer< ncm_video_frame > remove_from_save_queue1( );
		QSharedPointer< ncm_video_frame > remove_from_save_queue2( );*/

		//: Get the current qimage
		QImage* get_qimage( int cam_num );
		/*QImage* get_qimage1();
		QImage* get_qimage2();*/

		//: Set the current qimage
		void set_qimage( int cam_num, QSharedPointer< ncm_video_frame > _vid_frame );
		/*void set_qimage1( QSharedPointer< ncm_video_frame > _vid_frame );
		void set_qimage2( QSharedPointer< ncm_video_frame > _vid_frame );*/

		QImage* get_qdiffframe();
		void set_qdiffframe( QSharedPointer< ncm_video_frame > _vid_frame1, QSharedPointer< ncm_video_frame > _vid_frame2, QVector<bool> fliph, QVector<bool> flipv);

	protected:
		ncm_qcapture2_data_manager();

private:
		static ncm_qcapture2_data_manager* _instance;

		void set_qimage( QImage &qimage, QMutex &qimage_mutex, QSharedPointer< ncm_video_frame > _vid_frame );

		void convert_qimage(QImage &qimage, QSharedPointer< ncm_video_frame > _vid_frame );
		void diff_qimage(QImage &qimage, QSharedPointer< ncm_video_frame > _vid_frame1, QSharedPointer< ncm_video_frame > _vid_frame2, QVector<bool> fliph, QVector<bool> flipv);

		//Queue of shared points to frames that need displaying
		//QVector<ncm_video_frame_queue> main_queues_; 
		ncm_video_frame_queue main_queue1_;
		ncm_video_frame_queue main_queue2_;

		//Queue of shared points to frames that need saving
		//QVector<ncm_video_frame_queue> save_queues_; 
		
		ncm_video_frame_queue save_queue1_;
		ncm_video_frame_queue save_queue2_;
		
		//QVector <QImage> qimages_;
		QImage qimage1_;
		QImage qimage2_;
		QImage qdiffframe_;

		//QVector <QMutex> qimage_mutexes_;
		QMutex qimage_mutex1_;
		QMutex qimage_mutex2_;
		QMutex qdiffframe_mutex_;

};

#endif //NCM_QCAPTURE2_SAVER