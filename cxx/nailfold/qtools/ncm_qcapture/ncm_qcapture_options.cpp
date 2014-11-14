#include "ncm_qcapture_options.h"

// VXL libraries
#include <vcl_iostream.h>

//Nailfold libraries
#include "ncm_qcapture_preferences.h"

// Qt libraries
#include <QColorDialog>
#include <QtSql/QSqlRecord>
#include <QFileDialog>
#include <QMessageBox>
//
//: Constructor
ncm_qcapture_options::ncm_qcapture_options(ncm_qcapture_preferences & preferences,
																					QWidget *parent, Qt::WFlags flags)
: QDialog(parent, flags),
	current_preferences_(preferences)
{
  // setup the UI
  ui.setupUi(this);

	//General controls
	ui.autoLoadCamera->setChecked(current_preferences_.auto_load_camera_);
	ui.rotateFrames->setChecked(current_preferences_.rotate_frames_);
	ui.deleteEmptySessions->setChecked(current_preferences_.delete_empty_sessions_);
	ui.processUnacceptable->setChecked(current_preferences_.process_unacceptable_);
	
	//Save output settings
	ui.saveDirTextEdit->setText( current_preferences_.root_dir_.c_str() );
	ui.saveDirTextEdit->setReadOnly( true );
	ui.sequenceDataLineEdit->setText( current_preferences_.sequence_properties_name_.c_str() );
	ui.framePrefixLineEdit->setText( current_preferences_.frame_prefix_.c_str());
	
	ui.maxSaveFramesSpinBox->setValue( current_preferences_.max_frames_to_save_);
	ui.maxMosaicFramesSpinBox->setValue( current_preferences_.max_frames_to_align_);
	ui.maxMosaicSizeSpinBox->setValue( current_preferences_.max_mosaic_size_);
	ui.nthAlignFrameSpinBox->setValue( current_preferences_.align_nth_frame_ );

	//Autofocus controls
	ui.computeSharpness->setChecked(current_preferences_.compute_sharpness_);
	ui.sharpnessFrame->setEnabled(current_preferences_.compute_sharpness_);
	ui.enableAutofocus->setChecked(current_preferences_.use_autofocus_);
	ui.modulateFrame->setEnabled(current_preferences_.use_autofocus_);
	ui.modulateWeight->setValue(current_preferences_.modulate_weight_);

	//Motor controls
	ui.homeMotorsOnStartup->setChecked(current_preferences_.home_motors_on_startup_);
	ui.pixelResImageSpinBox->setValue(current_preferences_.pixels_per_mm_);
	ui.pixelResMotorSpinBox->setValue(current_preferences_.pixels_per_motor_mm_);

	ui.joystickSpeedMax->setValue(current_preferences_.joystick_max_);
	ui.dblclickSpeedMax->setValue(current_preferences_.dblclick_max_);
	ui.zarrowDist->setValue(current_preferences_.z_arrow_dist_);
	ui.zscrollDist->setValue(current_preferences_.z_scroll_dist_);
	
	//Now set the buttons false and unchanged flag (if you do this after any of the above
	//the widget's slot will enable the buttons and set the changed flag) 
	ui.cancelButton->setEnabled(false);
	ui.applyButton->setEnabled(false);
	preferences_changed_ = PreferencesUnchanged;

	ui.tabWidget->setCurrentWidget(ui.tabGeneral);
}

//
//: Destructor
ncm_qcapture_options::~ncm_qcapture_options()
{
}

//
//  Public slots
//

//
//: What to do when the user hits 'OK'

//
//  Private slots
//

//3 main controls

void ncm_qcapture_options::on_okButton_clicked()
{
	//Write values and close
	write_values();
	done(preferences_changed_);
}

void ncm_qcapture_options::on_cancelButton_clicked()
{
	//Close without writing values
	done(PreferencesUnchanged);
}
void ncm_qcapture_options::on_applyButton_clicked()
{
	//Write values without closing
	write_values();
	ui.cancelButton->setEnabled(false);
	ui.applyButton->setEnabled(false);
}

//---------------------------------------------------------
//General tab

void ncm_qcapture_options::on_autoLoadCamera_toggled(bool checked)
{
	new_preferences_.auto_load_camera_ = checked;
	preferences_modified();
}

void ncm_qcapture_options::on_rotateFrames_toggled(bool checked)
{
	new_preferences_.rotate_frames_ = checked;
	preferences_modified();
}

void ncm_qcapture_options::on_deleteEmptySessions_toggled(bool checked)
{
	new_preferences_.delete_empty_sessions_ = checked;
	preferences_modified();
}

void ncm_qcapture_options::on_processUnacceptable_toggled(bool checked)
{
	new_preferences_.process_unacceptable_ = checked;
	preferences_modified();
}

//Save options
void ncm_qcapture_options::on_saveDirSelect_clicked()
{
	// get filename from dialog box
	QString dirName = QFileDialog::getExistingDirectory(this, 
                       tr("Select folder to save images in"), "", 
                       QFileDialog::ShowDirsOnly);

  // Return if cancelled.
  if (dirName.isEmpty())
    return;

	dirName = dirName.replace("\\", "/");

	// Update the save_dir text edit box
	ui.saveDirTextEdit->setText( dirName );
	new_preferences_.root_dir_ = dirName.toStdString();

  preferences_modified();
}

void ncm_qcapture_options::on_sequenceDataLineEdit_textChanged(QString text)
{
	new_preferences_.sequence_properties_name_ = text.toStdString();
	preferences_modified();
}

void ncm_qcapture_options::on_framePrefixLineEdit_textChanged(QString text)
{
	new_preferences_.frame_prefix_ = text.toStdString();
	preferences_modified();
}

void ncm_qcapture_options::on_maxSaveFramesSpinBox_valueChanged(int value)
{
	new_preferences_.max_frames_to_save_ = value;
	preferences_modified();
}
void ncm_qcapture_options::on_maxMosaicFramesSpinBox_valueChanged(int value)
{
	new_preferences_.max_frames_to_align_ = value;
	preferences_modified();
}
void ncm_qcapture_options::on_maxMosaicSizeSpinBox_valueChanged(int value)
{
	new_preferences_.max_mosaic_size_ = value;
	preferences_modified();
}
void ncm_qcapture_options::on_nthAlignFrameSpinBox_valueChanged(int value)
{
	new_preferences_.align_nth_frame_ = value;
	preferences_modified();
}

//----------------------------------------------------------
//Auto focus tab
void ncm_qcapture_options::on_computeSharpness_toggled(bool checked)
{
	new_preferences_.compute_sharpness_ = checked;
	if (!checked)
	{
		ui.enableAutofocus->setChecked(false);
		ui.enableAutofocus->setEnabled(false);
	}
	else
		ui.enableAutofocus->setEnabled(true);

	preferences_modified();

}

void ncm_qcapture_options::on_nSharpnessSpinBox_valueChanged(int value)
{
	new_preferences_.n_sharpness_ = value;
	preferences_modified();
}

void ncm_qcapture_options::on_enableAutofocus_toggled(bool checked)
{
	new_preferences_.use_autofocus_ = checked;
	preferences_modified();
}

void ncm_qcapture_options::on_modulateWeight_valueChanged(int value)
{
	new_preferences_.modulate_weight_ = value;
	preferences_modified();
}

void ncm_qcapture_options::on_pixelResImageSpinBox_valueChanged(int value)
{
	new_preferences_.pixels_per_mm_ = value;
	preferences_modified();
}
void ncm_qcapture_options::on_pixelResMotorSpinBox_valueChanged(int value)
{
	new_preferences_.pixels_per_motor_mm_ = value;
	preferences_modified();
}

//home motors on startup
void ncm_qcapture_options::on_homeMotorsOnStartup_toggled(bool checked)
{
	new_preferences_.home_motors_on_startup_ = checked;
	preferences_modified();
}

//Set speeds of motors for X/Y panning
void ncm_qcapture_options::on_dblclickSpeedMax_valueChanged(double value)
{
	new_preferences_.dblclick_max_ = value;
	preferences_modified();
}

void ncm_qcapture_options::on_joystickSpeedMax_valueChanged(double value)
{
	new_preferences_.joystick_max_ = value;
	preferences_modified();
}

//Set speed/dist for Z motor
void ncm_qcapture_options::on_zarrowDist_valueChanged(double value)
{
	new_preferences_.z_arrow_dist_ = value;
	preferences_modified();
}

void ncm_qcapture_options::on_zscrollDist_valueChanged(double value)
{
	new_preferences_.z_scroll_dist_ = value;
	preferences_modified();
}

//
//  Private functions
//
void ncm_qcapture_options::write_values()
{
	current_preferences_ = new_preferences_;

}

void ncm_qcapture_options::preferences_modified()
{
	ui.cancelButton->setEnabled(true);
	ui.applyButton->setEnabled(true);
	preferences_changed_ = PreferencesChanged;
}