#ifndef NCM_QCAPTURE_SAVER
#define NCM_QCAPTURE_SAVER

#include <QtCore>
#include <vcl_string.h>


class ncm_qcapture_saver : public QObject
{
    Q_OBJECT

public: // Methods

	ncm_qcapture_saver();
	
	void setRootDir( vcl_string root_dir );
	void setStudyName( vcl_string study_name );
	void setSubjectName( vcl_string patient_name );
	void setSessionName( vcl_string session_name );
	void setSequenceName( vcl_string sequence_name );
	void setHand( vcl_string hand );
	void setDigit( vcl_string digit );
	void setMag( vcl_string mag );
	
	vcl_string makeSaveDir();

  //: Return true if there are more frames to save.
  bool is_saving() const;

public slots:

	void save_frame(int frame_num);

signals:

  void frame_saved(int);
	void saving_finished();

private: // Methods

	void update_frame_name();

private: // Variables

  // All of these can be set by the user
	vcl_string root_dir_;
	vcl_string study_name_;
	vcl_string subject_name_;	
	vcl_string session_name_;
	vcl_string sequence_name_;
  vcl_string hand_;
	vcl_string digit_;
	vcl_string mag_;
  vcl_string frame_prefix_;

	vcl_string frame_name_;
	vcl_string frame_dir_;

  //: Flag that indicates whether there are more frames to save.
  bool is_saving_;

};

#endif //NCM_QCAPTURE_SAVER