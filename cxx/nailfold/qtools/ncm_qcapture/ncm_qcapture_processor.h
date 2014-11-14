#ifndef NCM_QCAPTURE_PROCESSOR
#define NCM_QCAPTURE_PROCESSOR

#include <vcl_vector.h>

#include <QObject>

#include <nailfold/ncm_sharpness_evaluator.h>
#include <nailfold/qtools/ncm_video_frame.h>
#include "ncm_qcapture_preferences.h"

class ncm_qcapture_processor : public QObject
{
  Q_OBJECT

//  INTERFACE

public:
	ncm_qcapture_processor(ncm_qcapture_preferences & preferences);
	
	//: Set number of frames to save.
  void set_num_frames_to_save(int num_frames_to_save);

  //: Return the number of frames to be saved.
  int num_frames_to_save() const;

  //: Begin transferring images from the main queue to the save queue.
  void start_processing();

	//: Begin transferring images from the main queue to the save queue.
  void stop_saving();
  
  //: Return true if there are still frames to process.
  bool is_processing() const;

	//: Return true saving has been paused
  bool is_paused() const;

	//: Return true saving has been stopped
  bool is_stopped() const;

  //: Resize the sharpness array.
	void resize_sharpness_array();
  void resize_sharpness_array(int new_size);

  //: Return the sharpness of the current image.
  double image_sharpness() const;

	//:Set the paused flag
	void set_paused(bool paused);

  //: Set the lower threshold on sharpness.
  void set_lower_sharpness_threshold(double lower);

  //: Return the lower threshold on sharpness.
  double lower_sharpness_threshold() const;

  //: Set the lower threshold on sharpness.
  void set_upper_sharpness_threshold(double lower);

  //: Return the upper threshold on sharpness.
  double upper_sharpness_threshold() const;

signals:

	void frame_to_draw();
	void frame_to_save( ncm_video_frame_header header);
	void frame_to_align(int frame_num);
	void processing_finished(int n_processed);
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

  //: Add a sharpness value to the ring buffer.
  void store_sharpness(double sharpness);

  //: Update the sharpness score and check for threshold crossings.
  void update_sharpness();

private: // Variables

	//Reference to the main GUI's preferences
	ncm_qcapture_preferences & preferences_;

	//: Number of frames to save given current save sequence params
	int num_frames_to_save_;

  //: Number of frames left to process.
	int n_processed_;

  //: Class that computes sharpness of an image.
  ncm_sharpness_evaluator sharpness_evaluator_;

  //: Size of the sharpness array.
  //  (Must be declared before sharpness_array_ for initialization.)
  unsigned sharpness_array_size_; 

  //: Array of sharpness values that can be averaged to smooth the output.
  vcl_vector<double> sharpness_array_;

  //: Index into the sharpness buffer.
  unsigned sharpness_array_index_;

  //: Lower threshold for sharpness.
  double lower_sharpness_threshold_;

  //: Upper threshold for sharpness.
  double upper_sharpness_threshold_;

	//: Flag if saving paused
	bool is_paused_;

	//: Flag if saving stopped
	bool is_stopped_;
	
};

#endif //NCM_QCAPTURE_PROCESSOR