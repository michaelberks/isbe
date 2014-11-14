#ifndef ncm_qdata_manager_H
#define ncm_qdata_manager_H

#include <QtCore>
#include <QImage>

#include <vcl_cstring.h>
#include <vcl_sstream.h>
#include <vcl_iostream.h>
#include <vil/vil_image_view.h>

#include <nailfold/qtools/ncm_video_frame.h>
#include <nailfold/qtools/ncm_video_frame_queue.h>

class ncm_qdata_manager
{

// INTERFACE

public:

  //Controls access to unique instance for all clients
	static ncm_qdata_manager* Instance();

  // void add_stream(ncm_qcapture_data_stream stream);
  // ncm_qcapture_data_stream stream(int i);


	//: Return pointer to the main queue.
	ncm_video_frame_queue* main_queue();

	void add_to_main_queue( QSharedPointer< ncm_video_frame > );
	QSharedPointer< ncm_video_frame > remove_from_main_queue();

  //: Return pointer to save queue.
	ncm_video_frame_queue* save_queue();

	void add_to_save_queue( QSharedPointer< ncm_video_frame > );
	QSharedPointer< ncm_video_frame > remove_from_save_queue();

	QImage* get_qimage();
	void set_qimage( QSharedPointer< ncm_video_frame > _vid_frame );

	QImage* get_qlastframe();
	void set_qlastframe( QSharedPointer< ncm_video_frame > _vid_frame );



//  IMPLEMENTATION

protected:

  ncm_qdata_manager();


private:

  // vcl_vector<ncm_qcapture_data_stream> streams_;

	void convert_qimage( QImage &qimage, 
                       QSharedPointer< ncm_video_frame > _vid_frame );

  static ncm_qdata_manager* _instance;
	
	//: Queue of shared points to frames that need displaying
  ncm_video_frame_queue main_queue_; 

	//: Queue of shared points to frames that need saving
  ncm_video_frame_queue save_queue_;
	
  QImage qimage_;
	QImage qlastframe_;

	QMutex qimage_mutex_;
	QMutex qlastframe_mutex_;
};

#endif //NCM_QCAPTURE_SAVER