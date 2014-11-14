#ifndef NAILFOLD_QCAPTURE_PREFERENCES_H
#define NAILFOLD_QCAPTURE_PREFERENCES_H

//Nailfold includes:

//VXL includes:
#include <vcl_string.h>

//QT includes:
#include <QString>

//Forward class defs
class QSqlRecord;
class QSqlQuery;
class QDateTime;

class ncm_qcapture_preferences
{

public:
  ncm_qcapture_preferences();
  ~ncm_qcapture_preferences();

	//Maximum frames to capture and align
	int max_frames_to_save_;
	int max_frames_to_align_;
	int max_mosaic_size_;
	int align_nth_frame_;

	//General flags
	bool auto_load_camera_;
	bool delete_empty_sessions_;
	bool process_unacceptable_;

	//Pixel resolutions: both motor and real world
	double pixels_per_mm_;
	double pixels_per_motor_mm_;

	//Use SQL database?
	bool use_sql_database_;

	//Home the motors on startup
	bool home_motors_on_startup_;

	//Autofocus settings
	bool use_autofocus_;
	bool compute_sharpness_;
	double modulate_weight_;
	unsigned n_sharpness_; //Num frames to use in running average of sharpness

	//Speed values for X/Y movement
	double joystick_max_;
	double dblclick_max_;

	//Speed/distance values for Z movement
	double z_scroll_dist_;
	double z_arrow_dist_;

	//Save defaults properties name
	vcl_string root_dir_;
	vcl_string sequence_properties_name_;
	vcl_string frame_prefix_;

	//Display defaults
	bool rotate_frames_;

	bool sql_write(QSqlQuery& query, int user_id) const;
	bool sql_update(QSqlQuery& query, int user_id) const;
	bool sql_read(QSqlRecord& record);

protected:

private:

};

#endif // NAILFOLD_QCAPTURE_PREFERENCES_H
