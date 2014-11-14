#ifndef NCM_QCAPTURE_DATA_MANAGER_H
#define NCM_QCAPTURE_DATA_MANAGER_H

#include <QtCore>
#include <QImage>

#include <vcl_cstring.h>
#include <vcl_sstream.h>
#include <vcl_iostream.h>
#include <vil/vil_image_view.h>

#include <nailfold/qtools/ncm_video_frame.h>
#include <nailfold/qtools/ncm_video_frame_queue.h>

class ncm_qcapture_data_manager
{

// INTERFACE

public:

  //Controls access to unique instance for all clients
	static ncm_qcapture_data_manager* Instance();

  // void add_stream(ncm_qcapture_data_stream stream);
  // ncm_qcapture_data_stream stream(int i);


	//: Return pointer to the main queue.
	ncm_video_frame_queue* main_queue();

	void add_to_main_queue( QSharedPointer< ncm_video_frame > );
	QSharedPointer< ncm_video_frame > remove_from_main_queue();
	void empty_main_queue();

  //: Return pointer to save queue.
	ncm_video_frame_queue* save_queue();

	void add_to_save_queue( QSharedPointer< ncm_video_frame > );
	QSharedPointer< ncm_video_frame > remove_from_save_queue();

	//: Return pointer to align queue.
	ncm_video_frame_queue* align_queue();

	void add_to_align_queue( QSharedPointer< ncm_video_frame > );
	QSharedPointer< ncm_video_frame > remove_from_align_queue();

	//: Return pointer to mosaic queue.
	ncm_video_frame_queue* mosaic_queue();

	void add_to_mosaic_queue( QSharedPointer< ncm_video_frame > );
	QSharedPointer< ncm_video_frame > remove_from_mosaic_queue();

	QImage* get_qimage();
	void set_qimage( QSharedPointer< ncm_video_frame > _vid_frame );

	QImage* get_qlastframe();
	void set_qlastframe( QSharedPointer< ncm_video_frame > _vid_frame );



//  IMPLEMENTATION

protected:

  ncm_qcapture_data_manager();


private:

  // vcl_vector<ncm_qcapture_data_stream> streams_;

	void convert_qimage( QImage &qimage, 
                       QSharedPointer< ncm_video_frame > _vid_frame );

  static ncm_qcapture_data_manager* _instance;
	
	//: Queue of shared points to frames that need displaying
  ncm_video_frame_queue main_queue_; 

	//: Queue of shared points to frames that need saving
  ncm_video_frame_queue save_queue_;

	//: Queue of shared points to frames that need aligning in mosaic
  ncm_video_frame_queue align_queue_;

	//: Queue of shared points to frames that need aligning in mosaic
  ncm_video_frame_queue mosaic_queue_;
	
  QImage qimage_;
	QImage qlastframe_;

	QMutex qimage_mutex_;
	QMutex qlastframe_mutex_;
};

#endif //NCM_QCAPTURE_SAVER