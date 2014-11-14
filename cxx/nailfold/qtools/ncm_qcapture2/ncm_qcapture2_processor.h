#ifndef NCM_QCAPTURE2_PROCESSOR
#define NCM_QCAPTURE2_PROCESSOR

#include <QObject>
#include <QGraphicsPixmapItem>
#include <QGraphicsScene>
#include <QVector>

#include <vcl_iostream.h>
#include <vcl_iomanip.h>
#include <vcl_string.h>
#include <vcl_sstream.h>
#include <vcl_iosfwd.h>
#include <vil/vil_save.h>

#include <qcore/qcore_convert_image.h>

#include <nailfold/ncm_sharpness_evaluator.h>
#include <nailfold/qtools/ncm_video_frame.h>
#include "ncm_qcapture2_data_manager.h"

class ncm_qcapture2_processor : public QObject
{
  Q_OBJECT

//  INTERFACE

public:
	ncm_qcapture2_processor();
	
	//: Set number of frames to save.
  void set_num_frames_to_save(int num_frames_to_save, int cam_num );

  //: Return the number of frames to be saved.
  int num_frames_to_save( int cam_num ) const;

  //: Begin transferring images from the main queue to the save queue.
  void start_processing( int cam_num );
  
  //: Return true if there are still frames to process.
  bool is_processing(int cam_num) const;

  //: Return the sharpness of the current image.
  double image_sharpness() const;

  //: Set the lower threshold on sharpness.
  void set_lower_sharpness_threshold(double lower);

  //: Set the lower threshold on sharpness.
  void set_upper_sharpness_threshold(double lower);

signals:

	void diff_to_draw();
	void frame_to_draw( int cam_num );
	void frame_to_save1( int frame_num, int cam_num );
	void frame_to_save2( int frame_num, int cam_num );
	void processing_finished();
  void upper_sharpness_exceeded();
  void lower_sharpness_exceeded();

public slots:

  //: Pop a frame from the main queue, process it, and push it onto the save 
  //  queue.
	void process_frame( int cam_num, bool diff_frame, bool fliph, bool flipv );


//  IMPLEMENTATION

private: // Methods

  //: Emit a signal to indicate that all frames have been processed.
  void stop_processing();

  //: Update the sharpness score and check for threshold crossings.
  void update_sharpness(int cam_cum);

private: // Variables

	//: Number of frames to save given current save sequence params
	QVector<int> num_frames_to_save_;

  //: Number of frames left to process.
	QVector<int> n_processed_;

  //: Class that computes sharpness of an image.
  ncm_sharpness_evaluator sharpness_evaluator_;

  //: Sharpness of current image.
  double sharpness_;

  //: Lower threshold for sharpness.
  double lower_sharpness_threshold_;

  //: Upper threshold for sharpness.
  double upper_sharpness_threshold_;

	//Shared pointers to frames to compute difference image
	QSharedPointer< ncm_video_frame > diff_frame1_;
	QSharedPointer< ncm_video_frame > diff_frame2_;

	//: Bool to check diff frames are ready
	QVector< bool > diff_frame_ready_;

	//:Bool to check if either frame needs to be flipped in the diff image
	QVector< bool > diff_frame_fliph_;
	QVector< bool > diff_frame_flipv_;

};

/*
#include <QtCore>
#include <QGraphicsPixmapItem>
#include <QGraphicsScene>
#include <QVector>

#include <vcl_iostream.h>
#include <vcl_iomanip.h>
#include <vcl_string.h>
#include <vcl_sstream.h>
#include <vcl_iosfwd.h>
#include <vil/vil_save.h>

#include <nailfold/qtools/ncm_video_frame.h>
#include "ncm_qcapture2_data_manager.h"

class ncm_qcapture2_processor : public QObject
{
    Q_OBJECT

public:
	ncm_qcapture2_processor();
	~ncm_qcapture2_processor();

	//: Set the number of frames to be saved.
	void set_num_frames_to_save1(int num_frames_to_save);
	void set_num_frames_to_save2(int num_frames_to_save);

	//: Return the number of frames to be saved.
  int num_frames_to_save1() const;
  int num_frames_to_save2() const;

	//: Begin transferring images from the main queue to the save queue.
	void start_processing( bool use1, bool use2 );

	//: Return true if there are still frames to process.
	bool is_processing() const;

public slots:
	void process_frame1();
	void process_frame2();

	signals:
	void frame_to_draw1( );
	void frame_to_draw2( );
	void frame_to_save( int frame_num, int camera_num );
	void frame_to_save1( int frame_num );
	void frame_to_save2( int frame_num );
	void save_finished();

private:
	void process_frame( );

	//Flag to say if we're saving frames
	bool saving1_;
	bool saving2_;

	//Number of frames to save given current save sequence params
	int num_frames_to_save1_;
	int num_frames_to_save2_;

	//: Number of frames left to process.
	int n_processed1_;

	//Current frame number in save sequence
	int curr_frame1_;
	int curr_frame2_;
};
*/
#endif //NCM_QCAPTURE2_PROCESSOR