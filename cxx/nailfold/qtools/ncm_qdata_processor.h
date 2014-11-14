#ifndef __ncm_qdata_processor_h__
#define __ncm_qdata_processor_h__

#include <QObject>
#include <QImage>

#include <nailfold/ncm_sharpness_evaluator.h>

#include <nailfold/qtools/ncm_video_frame_queue.h>


class ncm_qdata_processor : public QObject
{
  Q_OBJECT

//  INTERFACE

public:
	ncm_qdata_processor();
	
	//: Set number of frames to save.
  void set_num_frames_to_save(int num_frames_to_save);

  //: Return the number of frames to be saved.
  int num_frames_to_save() const;

  //: Begin transferring images from the main queue to the save queue.
  void start_processing();
  
  //: Return true if there are still frames to process.
  bool is_processing() const;

  //: Return the sharpness of the current image.
  double image_sharpness() const;

  //: Set the lower threshold on sharpness.
  void set_lower_sharpness_threshold(double lower);

  //: Set the lower threshold on sharpness.
  void set_upper_sharpness_threshold(double lower);

  //: Define the input queue.
  void set_input_queue(ncm_video_frame_queue* input_queue);

  //: Define the output queue.
  void set_output_queue(ncm_video_frame_queue* output_queue);

  QImage* qimage();
	
	QImage* qlastframe();

signals:

	void frame_to_draw();
	void frame_to_save( int frame_num );
	void processing_finished();
  void upper_sharpness_exceeded();
  void lower_sharpness_exceeded();

public slots:

  //: Pop a frame from the main queue, process it, and push it onto the save 
  //  queue.
	void process_frame();


//  IMPLEMENTATION

private: // Methods

  //: Emit a signal to indicate that all frames have been processed.
  void stop_processing();

  //: Update the sharpness score and check for threshold crossings.
  void update_sharpness_of(const ncm_video_frame& video_frame);

  void set_qimage(QSharedPointer<ncm_video_frame> vid_frame );

  void set_qlastframe(QSharedPointer<ncm_video_frame> vid_frame );
	
private: // Variables

  ncm_video_frame_queue* input_queue_;

  ncm_video_frame_queue* output_queue_;

	//: Number of frames to save given current save sequence params
	int num_frames_to_save_;

  //: Number of frames left to process.
	int n_processed_;

  //: Class that computes sharpness of an image.
  ncm_sharpness_evaluator sharpness_evaluator_;

  //: Sharpness of current image.
  double sharpness_;

  //: Lower threshold for sharpness.
  double lower_sharpness_threshold_;

  //: Upper threshold for sharpness.
  double upper_sharpness_threshold_;

  //: Current frame
  QImage qimage_;

  //: Last frame in the sequence
	QImage qlastframe_;

  //: Mutex objects to prevent multiple accesses to data members.
	QMutex qimage_mutex_;
	QMutex qlastframe_mutex_;
};

#endif // __ncm_qdata_processor_h__