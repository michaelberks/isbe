#ifndef __ncm_qdata_stream_h__
#define __ncm_qdata_stream_h__

#include <nailfold/qtools/ncm_cameras/ncm_dmk_camera.h>

#include <nailfold/qtools/ncm_video_frame_queue.h>
#include <nailfold/qtools/ncm_qdata_processor.h>
#include <nailfold/qtools/ncm_qdata_saver.h>

class ncm_qdata_stream
{
public: // Methods

  //: Constructor.
  ncm_qdata_stream();

  //: Return pointer to camera.
	ncm_dmk_camera& camera();
	const ncm_dmk_camera& camera() const;

	//: Return pointer to the main queue.
	ncm_video_frame_queue& main_queue();
	const ncm_video_frame_queue& main_queue() const;

  //: Return pointer to processor.
	ncm_qdata_processor& processor();
	const ncm_qdata_processor& processor() const;

  //: Return pointer to save queue.
	ncm_video_frame_queue& save_queue();
	const ncm_video_frame_queue& save_queue() const;

  //: Return pointer to saver.
	ncm_qdata_saver& saver();
	const ncm_qdata_saver& saver() const;

public slots:

signals:

protected:

private: // Methods

private slots:

private: // Variables

	//: Camera object
	ncm_dmk_camera camera_;

	//: Queue of shared points to frames that need displaying
  ncm_video_frame_queue main_queue_; 

  //: Processor object
  ncm_qdata_processor processor_;

	//: Queue of shared points to frames that need saving
  ncm_video_frame_queue save_queue_;
	
  //: Saver object
  ncm_qdata_saver saver_;
};

#endif // __ncm_qdata_stream_h__
