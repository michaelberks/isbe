#ifndef __ncm_qdata_saver_h__
#define __ncm_qdata_saver_h__

#include <QtCore>

#include <vcl_string.h>

// Forward declarations.
class ncm_video_frame_queue;

class ncm_qdata_saver : public QObject
{
  Q_OBJECT

public: // Methods

	ncm_qdata_saver();
	
	void setRootDir( vcl_string root_dir );
  void setPatientName( vcl_string patient_name );
	void setSequenceName( vcl_string sequence_name );
	void setHand( vcl_string hand );
	void setDigit( vcl_string digit );
	void setMag( vcl_string mag );
	void setSubDir( vcl_string sub_dir );
	void setFramePrefix( vcl_string frame_prefix );
	void setFrameSuffix( vcl_string frame_suffix );

  void set_queue(ncm_video_frame_queue* queue);

	vcl_string makeSaveDir();

  //: Return true if there are more frames to save.
  bool is_saving() const;

  void set_properties_as_subdir( bool properties_as_subdir = true);


public slots:

	bool save_frame(int frame_num);


signals:

  void frame_saved(int);
	void saving_finished();


private: // Methods

	vcl_string output_dir();


private: // Variables

  // Properties common to all data streams
	static vcl_string root_dir_;
	static vcl_string patient_name_;
	static vcl_string sequence_name_;
  static vcl_string hand_;
	static vcl_string digit_;
	static vcl_string mag_;

  // Properties specific to this data stream
  vcl_string sub_dir_;
	vcl_string frame_prefix_;
	vcl_string frame_suffix_;

  //: Queue from which to pop frames for saving.
  ncm_video_frame_queue* queue_;

  //: Flag that indicates whether there are more frames to save.
  bool is_saving_;

  //: Flag to indicate whether the hand_digit_mag should be a subdirectory
  //  or part of every filename.
  bool properties_as_subdir_;
};

#endif // __ncm_qdata_saver_h__