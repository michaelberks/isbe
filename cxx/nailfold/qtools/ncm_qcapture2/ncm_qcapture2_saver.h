#ifndef NCM_QCAPTURE2_SAVER
#define NCM_QCAPTURE2_SAVER

#include <QtCore>
#include <vcl_iostream.h>
#include <vcl_iomanip.h>
#include <vcl_string.h>
#include <vcl_sstream.h>
#include <vcl_iosfwd.h>
#include <vil/vil_save.h>
#include <nailfold/qtools/ncm_video_frame.h>
#include "ncm_qcapture2_data_manager.h"

class ncm_qcapture2_saver : public QObject
{
    Q_OBJECT

public: // Methods

	ncm_qcapture2_saver();
	~ncm_qcapture2_saver();
	
	void setRootDir( vcl_string root_dir );
	void setSequenceName( vcl_string sequence_name );
	void setHand( vcl_string hand );
	void setDigit( vcl_string digit );
	void setMag( vcl_string mag );
	void setSubDir( vcl_string sub_dir );
	void setPatientName( vcl_string patient_name );
	void setCameraSuffix( vcl_string camera_suffix );

	void update_frame_name();

	vcl_string makeSaveDir();

  //: Return true if there are more frames to save.
  bool is_saving() const;

public slots:

	void save_frame(int frame_num, int cam_num);

signals:

  void frame_saved(int frame_num, int cam_num);
	void saving_finished(int cam_num);

private: // Methods

private: // Variables

  // All of these can be set by the user
	vcl_string root_dir_;
	vcl_string patient_name_;	
	vcl_string sequence_name_;
  vcl_string hand_;
	vcl_string digit_;
	vcl_string mag_;
  vcl_string sub_dir_;
	vcl_string frame_prefix_;
	vcl_string camera_suffix_;

	vcl_string frame_name_;
	vcl_string frame_dir_;

  //: Flag that indicates whether there are more frames to save.
  bool is_saving_;

};

#endif //NCM_QCAPTURE2_SAVER