#include "ncm_qcapture_preferences.h"

// VXL libraries
#include <vcl_iostream.h>

// Qt libraries
#include <QtSql/QSqlQuery>
#include <QtSql/QSqlError>
#include <QtSql/QSqlRecord>
#include <QVariant>
#include <Qstring>

//
//: Constructor
ncm_qcapture_preferences::ncm_qcapture_preferences()
	: max_frames_to_save_(1e5),
		max_frames_to_align_(1e4),
		max_mosaic_size_(10),
		align_nth_frame_(10),
		auto_load_camera_(true),
		delete_empty_sessions_(true),
		process_unacceptable_(false),
		pixels_per_mm_(1160),
		pixels_per_motor_mm_(1160),
		use_sql_database_(true),
		home_motors_on_startup_(true),
		joystick_max_(0.15),
		dblclick_max_(1.5),
		z_scroll_dist_(0.1),
		z_arrow_dist_(0.1),
		compute_sharpness_(false),
		use_autofocus_(false),
		modulate_weight_(10.0),
		n_sharpness_(15),
		sequence_properties_name_("sequence_properties.txt"),
		frame_prefix_("frame"),
		root_dir_("C:/isbe/nailfold/camera_capture"),
		rotate_frames_(true)
{
}

//
//: Destructor
ncm_qcapture_preferences::~ncm_qcapture_preferences()
{
}


bool ncm_qcapture_preferences::sql_write(QSqlQuery& query, int user_id) const
{
	query.prepare("INSERT INTO user_preferences (user_id) VALUES (:user_id)");
	query.bindValue(":user_id", user_id);
	bool success = query.exec();

	//If successful write, get back the auto-incremented ID (this is why func can't be const)
	if (!success)
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;
	
	success = sql_update(query, user_id);

	return success;
}

bool ncm_qcapture_preferences::sql_update(QSqlQuery& query, int user_id) const
{
	query.prepare(
		"UPDATE user_preferences SET "
			"max_frames_to_save = :max_frames_to_save, "
			"max_frames_to_align = :max_frames_to_align, "
			"align_nth_frame = :align_nth_frame, "
			"auto_load_camera = :auto_load_camera, "
			"delete_empty_sessions = :delete_empty_sessions, "
			"process_unacceptable = :process_unacceptable, "
			"pixels_per_mm = :pixels_per_mm, "
			"pixels_per_motor_mm = :pixels_per_motor_mm, "
			"home_motors_on_startup = :home_motors_on_startup, "
			"use_autofocus = :use_autofocus, "
			"compute_sharpness = :compute_sharpness, "
			"modulate_weight = :modulate_weight, "
			"n_sharpness = :n_sharpness, "
			"joystick_max = :joystick_max, "
			"dblclick_max = :dblclick_max, "
			"z_scroll_dist = :z_scroll_dist, "
			"z_arrow_dist = :z_arrow_dist, "
			"root_dir = :root_dir, "
			"sequence_properties_name = :sequence_properties_name, "
			"frame_prefix = :frame_prefix "
			"WHERE user_id = " + QString::number(user_id));
	
	query.bindValue(":max_frames_to_save", max_frames_to_save_);
	query.bindValue(":max_frames_to_align", max_frames_to_align_);
	query.bindValue(":align_nth_frame", align_nth_frame_);

	//General flags
	query.bindValue(":auto_load_camera", auto_load_camera_);
	query.bindValue(":delete_empty_sessions", delete_empty_sessions_);
	query.bindValue(":process_unacceptable", process_unacceptable_);

	//Pixel resolutions: both motor and real world
	query.bindValue(":pixels_per_mm", pixels_per_mm_);
	query.bindValue(":pixels_per_motor_mm", pixels_per_motor_mm_);

	//Home the motors on startup
	query.bindValue(":home_motors_on_startup", home_motors_on_startup_);

	//Autofocus settings
	query.bindValue(":use_autofocus", use_autofocus_);
	query.bindValue(":compute_sharpness", compute_sharpness_);
	query.bindValue(":modulate_weight", modulate_weight_);
	query.bindValue(":n_sharpness", n_sharpness_); //Num frames to use in running average of sharpness

	//Speed values for X/Y movement
	query.bindValue(":joystick_max", joystick_max_);
	query.bindValue(":dblclick_max", dblclick_max_);

	//Speed/distance values for Z movement
	query.bindValue(":z_scroll_dist", z_scroll_dist_);
	query.bindValue(":z_arrow_dist", z_arrow_dist_);

	//Save defaults properties name
	query.bindValue(":root_dir", QString(root_dir_.c_str()));
	query.bindValue(":sequence_properties_name", QString(sequence_properties_name_.c_str()));
	query.bindValue(":frame_prefix", QString(frame_prefix_.c_str()));
  
  bool success = query.exec();

	//If successful write, get back the auto-incremented ID (this is why func can't be const)
	if (!success)
	{
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;
	}

	return success;
}

bool ncm_qcapture_preferences::sql_read(QSqlRecord& record)
{

	max_frames_to_save_ = record.value("max_frames_to_save").toInt();
	max_frames_to_align_ = record.value("max_frames_to_align").toInt();
	align_nth_frame_ = record.value("align_nth_frame").toInt();

	//General flags
	auto_load_camera_ = record.value("auto_load_camera").toBool();
	delete_empty_sessions_ = record.value("delete_empty_sessions").toBool();
	process_unacceptable_ = record.value("process_unacceptable").toBool();

	//Pixel resolutions: both motor and real world
	pixels_per_mm_ = record.value("pixels_per_mm").toDouble();
	pixels_per_motor_mm_ = record.value("pixels_per_motor_mm").toDouble();

	//Home the motors on startup
	home_motors_on_startup_ = record.value("home_motors_on_startup").toBool();

	//Autofocus settings
	use_autofocus_ = record.value("use_autofocus").toBool();
	compute_sharpness_ = record.value("compute_sharpness").toBool();
	modulate_weight_ = record.value("modulate_weight").toDouble();
	n_sharpness_ = record.value("n_sharpness").toUInt(); //Num frames to use in running average of sharpness

	//Speed values for X/Y movement
	joystick_max_ = record.value("joystick_max").toDouble();
	dblclick_max_ = record.value("dblclick_max").toDouble();

	//Speed/distance values for Z movement
	z_scroll_dist_  = record.value("z_scroll_dist").toDouble();
	z_arrow_dist_ = record.value("z_arrow_dist").toDouble();

	//Save defaults properties name
	root_dir_ = record.value("root_dir").toString().toStdString();
	sequence_properties_name_ = record.value("sequence_properties_name").toString().toStdString();
	frame_prefix_ = record.value("frame_prefix").toString().toStdString();

	return true;
}