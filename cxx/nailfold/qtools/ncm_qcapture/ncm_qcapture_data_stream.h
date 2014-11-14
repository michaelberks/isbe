#ifndef __ncm_qcapture_data_stream_h__
#define __ncm_qcapture_data_stream_h__

class ncm_qcapture_data_stream
{
public: // Methods

  

public slots:

signals:

protected:

private: // Methods

  //: Convert image from ncm_video_frame to QImage.
	void convert_qimage( QImage &qimage, 
                       QSharedPointer< ncm_video_frame > _vid_frame );

private slots:

private: // Variables

	//: Camera object
	ncm_dmk_camera camera_;

  //: Saver object
  ncm_qcapture_saver saver_;

  //: Processor object
  ncm_qcapture_processor processor_;

	//: Queue of shared points to frames that need displaying
  ncm_video_frame_queue main_queue_; 

	//: Queue of shared points to frames that need saving
  ncm_video_frame_queue save_queue_;
	
  QImage qimage_;
	QImage qlastframe_;

	QMutex qimage_mutex_;
	QMutex qlastframe_mutex_;
};

#endif // __ncm_qcapture_data_stream_h__
