#include "ncm_qcapture_gui.h"

#include <dbt.h>

#include <vcl_iomanip.h>
#include <vcl_sstream.h>

#include <vil/vil_load.h>
#include <vil/vil_convert.h>
#include <vil/vil_save.h>

#include <vsl/vsl_quick_file.h>
#include <vnl/vnl_vector.h>
#include <vul/vul_file.h>

#include <qcore/qcore_convert_image.h>

#include <nailfold/ncm_apt_controller.h>
#include "ncm_qcapture_search.h"
#include "ncm_qcapture_end_sequence.h"
#include "ncm_qcapture_end_session.h"
#include "ncm_qcapture_options.h"

#include <QtSql/QSQlQuery>
#include <QtSql/QSQlRecord>
#include <QtSql/QSQlError>
#include <QVariant>
//
//: Constructor
ncm_qcapture_gui::ncm_qcapture_gui(ncm_qcapture_preferences& preferences,
                                    QWidget *parent, Qt::WFlags flags)
: QMainWindow(parent, flags),
	preferences_(preferences),
  velocityX_(0.0f),
  velocityY_(0.0f),
	velocityZ_(0.0f),
  // Constants
  motors_are_live_(true),
  autoFocus_state_(autoFocus_Inactive),
	calibrate_state_(calibrate_Inactive),
	activity_state_(activity_Login),
	motor_positions_(0, vcl_vector<double>(3)),
	motor_times_(0),
	mosaic_frame_positions_(0),
	processor_(preferences)
{
  // setup the UI
  ui.setupUi(this);

	//Load in the start/pause icons
	start_icon_.addFile(QString::fromUtf8(":/images/play.ico"), QSize(), QIcon::Normal, QIcon::Off);
	pause_icon_.addFile(QString::fromUtf8(":/images/pause.ico"), QSize(), QIcon::Normal, QIcon::Off);
	update_session_controls();

	//Now we don't do anything else until we've successfully loaded a user	
	
}

//
//: Destructor
ncm_qcapture_gui::~ncm_qcapture_gui()
{
}

 
// Events
void ncm_qcapture_gui::closeEvent(QCloseEvent *ev)
{
	if (activity_state_ == activity_Recording || activity_state_ == activity_RecordingPaused)
  {
    QMessageBox msgBox;
    msgBox.setText("Cannot close during live recording.");
    msgBox.setIcon(QMessageBox::Information);
    msgBox.exec();
    ev->ignore();
    return;
  }

	else if (activity_state_ == activity_Session)
  {
    QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.setText("Imaging session is still open");
    msgBox.setInformativeText("Close session and exit software?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::Yes);

    const int response = msgBox.exec();

    if (response == QMessageBox::No)
		{
			ev->ignore();
			return;
		}
		else
			ui.endSessionPushButton->click();
  }

  if (saver_.is_saving())
  {
    QMessageBox msgBox;
    msgBox.setText("Still saving frames - please wait.");
    msgBox.setIcon(QMessageBox::Information);
    msgBox.exec();
    ev->ignore();
    return;
  }

	/* This tidying up is done in end session now
	disconnect_database();
	disconnect_camera();
  apt_timer_.stop();

	// Close the threads
	saver_thread_.quit();
	processor_thread_.quit();
  apt_thread_.quit();
	aligner_thread_.quit();*/


  ev->accept();
}

bool ncm_qcapture_gui::winEvent(MSG * message, long * result)
{
	if (message->message == WM_DEVICECHANGE)
	{
		switch (message->wParam)
		{
		case DBT_CONFIGCHANGECANCELED:
			vcl_cout << "Device change: DBT_CONFIGCHANGECANCELED" << vcl_endl;
			break;

		case DBT_CONFIGCHANGED:
			vcl_cout << "Device change: DBT_CONFIGCHANGED" << vcl_endl;
			break;

		case DBT_DEVICEARRIVAL:
			vcl_cout << "Device change: DBT_DEVICEARRIVAL" << vcl_endl;
			break;

		case DBT_DEVICEQUERYREMOVE:
			vcl_cout << "Device change: DBT_DEVICEQUERYREMOVE" << vcl_endl;
			break;

		case DBT_DEVICEQUERYREMOVEFAILED:
			vcl_cout << "Device change: DBT_DEVICEQUERYREMOVEFAILED" << vcl_endl;
			break;

		case DBT_DEVICEREMOVECOMPLETE:
			vcl_cout << "Device change: DBT_DEVICEREMOVECOMPLETE" << vcl_endl;
			break;

		case DBT_DEVICEREMOVEPENDING:
			vcl_cout << "Device change: DBT_DEVICEREMOVEPENDING" << vcl_endl;
			break;

		case DBT_DEVICETYPESPECIFIC:
			vcl_cout << "Device change: DBT_DEVICETYPESPECIFIC" << vcl_endl;
			break;

		case DBT_DEVNODES_CHANGED:
			vcl_cout << "Device change: DBT_DEVNODES_CHANGED" << vcl_endl;
			break;

		case DBT_QUERYCHANGECONFIG:
			vcl_cout << "Device change: DBT_QUERYCHANGECONFIG" << vcl_endl;
			break;

		case DBT_USERDEFINED:
			vcl_cout << "Device change: DBT_USERDEFINED" << vcl_endl;
			break;

		default:
			vcl_cout << "Device change not recognised" << vcl_endl;

		}
		result = 0;
		return true;
	}
	else
		return false;
}

//
//  Private slots
// 

void ncm_qcapture_gui::on_actionExit_triggered()
{
	this->close();
}

//Log in to the system - cant do owt until this is achieved
void ncm_qcapture_gui::on_logInPushButton_clicked()
{
	QString username =  ui.usernameEdit->text();
	QString password = ui.passwordEdit->text();
	ui.passwordEdit->setText("");

	bool db_connected = connect_to_database(username, password);

	//What do we do if database connection unsuccessful?
	if (db_connected)
	{
		//Load user (this will also load their preferences from SQL)
		if(!load_user(username))
		{
			//Tell the user what's up
			QMessageBox msgBox;
			msgBox.setIcon(QMessageBox::Warning);
			msgBox.setText(
				"Could not find user in database. \n"
				"You will be able to run an image session as the default user, \n"
				"but will not be able to save any changes you make to the software preferences\n"
				"");
			msgBox.exec();
		}

		log_in();
		
	}
	else
	{
		//Tell the user what's up
		QMessageBox msgBox;
		msgBox.setIcon(QMessageBox::Warning);
		msgBox.setText(
			"Could not connect to database. \n"
			"Please check your username and password are correct. \n"
			"If you are certain they are, there may be a \n"
			"problem with the database and you need to contact \n"
			"your software administrator \n");
		msgBox.exec();
	}
}

//-------------------------------------------------------------------
// File menu call backs

void ncm_qcapture_gui::on_actionLoad_camera_triggered()
{
	load_camera(true);
}
void ncm_qcapture_gui::on_actionChange_camera_triggered()
{
	load_camera(false);
}

//
//: User selected 'Preferences' so bring up dialog box
void ncm_qcapture_gui::on_actionPreferences_triggered()
{
  ncm_qcapture_options optionsWindow(preferences_, this);
  const int response = optionsWindow.exec();

	if (response == ncm_qcapture_options::PreferencesChanged)
	{
		//Ask the user if we want to save the changed preferences for next time
		QMessageBox msgBox;
		msgBox.setIcon(QMessageBox::Question);
		msgBox.setText(
			"Save preferences for next time?");
		msgBox.setInformativeText(
			"If you choose 'No' your new preferences will apply for the rest of \n"
			"this session but will revert to your previous saved preferences the \n"
			"next time you run the software.");
		msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
		msgBox.setDefaultButton(QMessageBox::Yes);

		const int response = msgBox.exec();

		if (response == QMessageBox::Yes)
		{
			QSqlQuery query(database_);
			preferences_.sql_update(query, user_.user_id());
		}

		//Update a couple of bits and bits that may have changed - could be a waste
		//if these specific fields haven't actually changed but it's hardly expensive
		saver_.setRootDir( preferences_.root_dir_ );
		processor_.resize_sharpness_array();
		update_tab_controls();
	}
}

//----------------------------------------------------------------------------
// Welcome page control callbacks
void ncm_qcapture_gui::on_newSessionPushButton_clicked()
{
	start_image_session();
}

void ncm_qcapture_gui::on_processPushButton_clicked()
{
	start_offline_processing();
}

void ncm_qcapture_gui::on_reviewPushButton_clicked()
{
	start_review_session();
}

void ncm_qcapture_gui::on_logOffPushButton_clicked()
{
	log_off();
}

//----------------------------------------------------------------------------
// Image session callbacks
void ncm_qcapture_gui::on_startRecordingPushButton_clicked()
{
	if (activity_state_ == activity_Session)
	{

		//Create a new sequence
		ncm_image_sequence sequence(subject_.live_session()->session_id());
		sequence.start();

		//Set camera properties
		sequence.set_sequence_name(sequence.time_started().toString("hh_mm_ss"));
		sequence.set_camera_type(camera_.get_camera_type().c_str());
		sequence.set_camera_id(camera_.get_camera_id().c_str());
		sequence.set_frame_rate(camera_.get_current_FPS());
		if ( camera_.is_auto_exposure() )
			sequence.set_exposure_time(-1);
		else
			sequence.set_exposure_time(camera_.get_current_exposure());

		if ( camera_.is_auto_gain() )
			sequence.set_gain(-1);
		else
			sequence.set_gain(camera_.get_current_gain());		

		//Set motor properties
		sequence.set_apt_x_id(apt_.X()->id());
		sequence.set_apt_y_id(apt_.Y()->id());
		sequence.set_apt_z_id(apt_.Z()->id());
		sequence.set_x_reversed(apt_.X()->reversed());
		sequence.set_y_reversed(apt_.Y()->reversed());
		sequence.set_z_reversed(apt_.Z()->reversed());

		//Add the sequence name to the saver and make a new directory for the images
		saver_.setRootDir( preferences_.root_dir_ );
		saver_.setSequenceName(sequence.sequence_name().toStdString());
		vcl_string save_dir = saver_.makeSaveDir();

		//Add the filename of the datafile sequence object and then add the sequence object to the session
		sequence.set_images_datafile(QString(preferences_.sequence_properties_name_.c_str()));
		QString qsave_dir = QString(save_dir.c_str());
		sequence.set_images_full_dir(qsave_dir);

		subject_.live_session()->add_image_sequence(sequence);
		subject_.live_session()->display_details();

		//Make sure the aligner has up to date pixel resolution
		aligner_.setPixelsPerMm(preferences_.pixels_per_motor_mm_);
		aligner_.resetDisplacements();
		current_mosaic_displacement_ = QPoint(0, 0);

		//Turn on recording flag
		activity_state_ = activity_Recording;
		
		// Set the number of frames to save to MAX_FRAMES (achieved by sending a negative int)
		//and tell the processor to start
		processor_.set_num_frames_to_save( -1 );
		processor_.start_processing();
		
	}
	else if (activity_state_ == activity_Recording)
	{
		//Tell the processor to pause
		processor_.set_paused(true);
		activity_state_ = activity_RecordingPaused;
	}
	else if (activity_state_ == activity_RecordingPaused)
	{
		//Tell the processor to start again
		processor_.set_paused(false);
		activity_state_ = activity_Recording;
	}
	update_session_controls();
}

void ncm_qcapture_gui::on_stopRecordingPushButton_clicked()
{
	if (!(activity_state_ == activity_Recording || activity_state_ == activity_RecordingPaused)) //The button should be disabled anyway but just in case...
		return;

	//Set the activity state to saving to disable GUI until saving complete
	activity_state_ = activity_Saving;
	update_session_controls();

	//Tell the processor to stop saving - this will trigger signals
	//that will eventually call onProcessingFinished() where we finalise the sequence
	//and write it to SQL
	processor_.stop_saving();
}

void ncm_qcapture_gui::on_endSessionPushButton_clicked()
{
	end_image_session();
}

//----------------------------------------------------------------------------
//Reviewing session callbacks
void ncm_qcapture_gui::on_finishedReviewButton_clicked()
{
	end_review_session();
}

//----------------------------------------------------------------------------
// Display control callbacks
void ncm_qcapture_gui::on_brightnessSlider_sliderMoved(int brightness)
{
  scene_.set_contrast_from(ui.contrastSlider->value(), brightness);
}

void ncm_qcapture_gui::on_contrastSlider_sliderMoved(int contrast)
{
  scene_.set_contrast_from(contrast, ui.brightnessSlider->value());
}

void ncm_qcapture_gui::on_stopLiveButton_clicked( )
{
	if ( !ui.stopLiveButton->text().compare( QString( "Start live" ) ) ) 
  {
		camera_.start_live();
		ui.stopLiveButton->setText( QString( "Stop live" ) );
	}
	else 
  {
		camera_.stop_live();
		ui.stopLiveButton->setText( QString( "Start live" ) );
	}
}

//Display control callbacks for reviewing
//
//: Brightness and contrast controls
void ncm_qcapture_gui::on_brightnessRSlider_sliderMoved(int value)
{
  const int contrast = ui.contrastRSlider->value();
  const int brightness = value;
  const int lower_limit = (255-brightness) - (255-contrast)/2;
  const int upper_limit = (255-brightness) + (255-contrast)/2;

  imageProcessor_.set_contrast(lower_limit, upper_limit);
}

void ncm_qcapture_gui::on_contrastRSlider_sliderMoved(int value)
{
  const int contrast = value;
  const int brightness = ui.brightnessRSlider->value();
  const int lower_limit = (255-brightness) - (255-contrast)/2;
  const int upper_limit = (255-brightness) + (255-contrast)/2;

  imageProcessor_.set_contrast(lower_limit, upper_limit);
}

//:
//Zoom controls for reviewing

//
//: Update zoom controls
void ncm_qcapture_gui::updateZoomControls()
{
  const double scale = ui.reviewView->transform().m11();
  ui.zoomSlider->setValue(static_cast<int>(100*scale + 1));

  QString zoom_string = QString::number(static_cast<int>(100*scale + 1));
  ui.zoomEdit->setText(zoom_string+"%");
}

//
//: Update limits on zoom controls
void ncm_qcapture_gui::updateZoomLimits()
{
  const double minimum_scale = ui.reviewView->fitHeightScale();
  ui.zoomSlider->setMinimum(static_cast<int>(100*minimum_scale + 1));

  const double maximum_scale = 10 * minimum_scale;
  ui.zoomSlider->setMaximum(static_cast<int>(100*maximum_scale - 1));

  ui.reviewView->setMaxScale(maximum_scale);
}

void ncm_qcapture_gui::on_zoomSlider_sliderMoved(int value)
{
  setZoomFrom(value);
}

void ncm_qcapture_gui::on_zoomEdit_editingFinished()
{
  // Remove last character (%) from string and convert
  QString value_string = ui.zoomEdit->text();
  value_string.chop(1);
  const int value = value_string.toInt();

  setZoomFrom(value);
}

void ncm_qcapture_gui::on_zoomOutButton_clicked()
{
  // Avoid calling functions that generate more signals, such as this:
  //ui.zoomSlider->setValue(ui.zoomSlider->value()-ui.zoomSlider->pageStep());

  // Instead, set the graphicsView zoom directly, and let the resulting signal
  // trigger an update of the GUI controls
  int new_zoom = ui.zoomSlider->value() - ui.zoomSlider->pageStep();
  setZoomFrom(new_zoom);
}

void ncm_qcapture_gui::on_zoomInButton_clicked()
{
  // Avoid calling functions that generate more signals, such as this:
  //ui.zoomSlider->setValue(ui.zoomSlider->value()+ui.zoomSlider->pageStep());

  // Instead, set the graphicsView zoom directly, and let the resulting signal
  // trigger an update of the GUI controls
  int new_zoom = ui.zoomSlider->value() + ui.zoomSlider->pageStep();
  setZoomFrom(new_zoom);
}

void ncm_qcapture_gui::applyUserZoom()
{
  /*switch (scene_.edit_mode())
  {
    case ncm_qscene::ModeClassifyImage:
      applyImageZoom(preferences_.label_image_zoom());
      break;
    case ncm_qscene::ModeAddVessels:
      applyImageZoom(preferences_.add_vessels_zoom());
      break;
    case ncm_qscene::ModeSetVesselSize: // fall through
    case ncm_qscene::ModeSetVesselShape:
      applyImageZoom(preferences_.set_vessel_props_zoom());
      break;
    case ncm_qscene::ModeLabelApices: // fall through
    case ncm_qscene::ModeDrawVesselPath:
      applyVesselZoom(ui.vesselCombo->currentIndex());
      break;
    case ncm_qscene::ModeAddHaemorrhages:
      applyImageZoom(preferences_.add_haemorrhages_zoom());
      break;
		case ncm_qscene::ModeDisplayAll:
			applyImageZoom(preferences_.label_image_zoom());
      break;
    default:
      assert(false);
  }*/
}
void ncm_qcapture_gui::applyVesselZoom(int vessel_index)
{
  /*ncm_vessel* vessel = markup_.vessel(vessel_index);
  if (vessel == NULL)
    return;

  // Get normal zoom first
  double zoom = preferences_.normal_vessel_zoom(); // default
  ncm_vessel_properties vessel_props = vessel->properties();
  
  // Zoom further out for larger vessels
  if (vessel_props.is_size_enlarged())
    zoom *= preferences_.enlarged_zoom_relative();
  else if (vessel_props.is_size_giant() ||
           vessel_props.is_size_irregular())
  {
    zoom *= preferences_.enlarged_zoom_relative();
    zoom *= preferences_.giant_zoom_relative();
  }

  // Apply a shift if we're drawing the vessel path
  double y_shift = 0.0;
  if (scene_.edit_mode() == ncm_qscene::ModeDrawVesselPath)
    y_shift = preferences_.vessel_path_yshift();

  ui.graphicsView->focus_on_vessel(vessel_index, zoom, y_shift);*/
}

void ncm_qcapture_gui::on_autoContrast_toggled(bool checked)
{
	ui.reviewView->setAutoContrast(checked);
	updateContrastControls();
}

//
//: Update brightness and contrast that may change as a result of autocontrast
void ncm_qcapture_gui::updateContrastControls()
{
  const int min_contrast = imageProcessor_.min_contrast();
  const int max_contrast = imageProcessor_.max_contrast();
  const int brightness = 255 - (min_contrast + max_contrast)/2;
  const int contrast = 255 - max_contrast + min_contrast;

  ui.brightnessRSlider->setValue(brightness);
  ui.contrastRSlider->setValue(contrast);

  ui.brightnessRSlider->setEnabled(!ui.autoContrast->isChecked());
  ui.contrastRSlider->setEnabled(!ui.autoContrast->isChecked());
}

//----------------------------------------------------------------------------
// Camera control callbacks
void ncm_qcapture_gui::on_frameRateSelection_currentIndexChanged( 
  const QString &FpsText ) 
{
	// Stop live mode, change the FPS, then restart live mode
	const double frames_per_sec = FpsText.toDouble();
	camera_.set_FPS( frames_per_sec );
}

void ncm_qcapture_gui::on_gainSpinBox_valueChanged( double value )
{
	ui.gainSlider->setValue( int( value ) );
	camera_.set_gain( value );

}
void ncm_qcapture_gui::on_gainSlider_sliderMoved( int value )
{
	ui.gainSpinBox->setValue( double( value ) );
	camera_.set_gain( double( value ) );
}
void ncm_qcapture_gui::on_autoGain_stateChanged( int state )
{
	if ( state )
	{
		//disable slider and spin box
		ui.gainSlider->setEnabled( false );
		ui.gainSpinBox->setEnabled( false );

		//switch on camera's auto gain
		camera_.switch_auto_gain( true );

		vcl_cout << "Auto gain switched on" << vcl_endl;
	}
	else
	{
		//switch off camera's auto gain and get current gain value
		double gain;
		camera_.switch_auto_gain( false );
		camera_.get_current_gain( gain );

		//enable gain slider and spin box and set to current auto value
		ui.gainSlider->setEnabled( true );
		ui.gainSpinBox->setEnabled( true );
		ui.gainSpinBox->setValue( gain );
		ui.gainSlider->setValue( int( gain ) );
		vcl_cout << "Auto gain switched off" << vcl_endl;
	}
}
void ncm_qcapture_gui::on_exposureSpinBox_valueChanged( double value )
{
	double exposure = ui.exposureSpinBox->exposureFromValue( value );
	ui.exposureSlider->setValue( static_cast<int>(value) );
	camera_.set_exposure( exposure );
}
void ncm_qcapture_gui::on_exposureSlider_sliderMoved( int value )
{
	double exposure = ui.exposureSpinBox->exposureFromValue( value );
	ui.exposureSpinBox->setValue(static_cast<double>(value));
	camera_.set_exposure( exposure );
}
void ncm_qcapture_gui::on_autoExposure_stateChanged( int state )
{
	if ( state )
	{
		//disable slider and spin box
		ui.exposureSlider->setEnabled( false );
		ui.exposureSpinBox->setEnabled( false );

		//switch on camera's auto exposure
		camera_.switch_auto_exposure( true );

		vcl_cout << "Auto exposure switched on" << vcl_endl;
	}
	else
	{
		//switch off camera's auto exposure and get current exposure value
		double exposure, value;
		camera_.switch_auto_exposure( false );
		camera_.get_current_exposure( exposure );
		value = ui.exposureSpinBox->valueFromExposure( exposure );

		//enable exposure slider and spin box and set to current auto value
		ui.exposureSlider->setEnabled( true );
		ui.exposureSpinBox->setEnabled( true );
		ui.exposureSpinBox->setValue( value );
		ui.exposureSlider->setValue( int( value ) );
		vcl_cout << "Auto exposure switched off" << vcl_endl;
		vcl_cout << "Current exposure: " << exposure << vcl_endl;
	}
}

 
//: Add a line to the properties file summarizing the current frame.
void ncm_qcapture_gui::add_properties_header()
{
	/*
  // Write the frame number, time stamp, motor position (x, y, z), sharpness,
  // and so on.
  property_fs_ << vcl_endl
               << vcl_setw(12) << "Name"
							 << vcl_setw(12) << "No."
               << vcl_setw(12) << "Time"
               << vcl_setw(12) << "MotorX"
               << vcl_setw(12) << "MotorY"
               << vcl_setw(12) << "MotorZ"
               << vcl_setw(12) << "Sharpness"
               << vcl_endl;*/
}
 
//: Add metadata to the most recently captured frame.
void ncm_qcapture_gui::tag_frame()
{
	QSharedPointer<ncm_video_frame>& frame = 
    ncm_qcapture_data_manager::Instance()->main_queue()->tail();

	double mx = apt_.X()->position();
	double my = apt_.Y()->position();
	double mz = apt_.Z()->position();
	QTime t = QTime::currentTime();

  frame->set_motor_xyz( mx, my, mz);

	if (activity_state_ == activity_Recording)
	{
		vcl_vector<double> mv(3);
		mv.push_back(mx);
		mv.push_back(my);
		mv.push_back(mz);
		motor_positions_.push_back(mv);
		motor_times_.push_back(t.msecsTo(motor_time0_));
		//vcl_cout << "Time lag in tagging frames: " << t.msecsTo(frame->get_frame_header().frame_time_) << vcl_endl;
	}
  emit frame_tagged();
}
 
//: Add a line to the properties file summarizing the current frame.
void ncm_qcapture_gui::log_frame_properties(ncm_video_frame_header header)
{
	//Add this frame to the sequence object
	ncm_image_sequence &sequence = subject_.live_session()->image_sequences()->back();
	sequence.add_frame(header);
	int frame_number = header.frame_num_;

	/*
  if (!property_fs_.is_open())
    return;

	
	vcl_string frame_name = header.frame_name_;

	float x = header.motor_position_[0];
	float y = header.motor_position_[1];
	float z = header.motor_position_[2];

  double sharpness = header.sharpness_;

  vcl_string time_str = "";
	time_str = header.frame_time_.toString("hhmmsszzz").toStdString();

	
  // Write the frame number, time stamp, motor position (x, y, z), sharpness,
  // and so on.
  property_fs_ << vcl_setw(12) << frame_name
							 << vcl_setw(12) << frame_number
               << vcl_setw(12) << time_str
               << vcl_setw(12) << x
               << vcl_setw(12) << y
               << vcl_setw(12) << z
               << vcl_setw(12) << sharpness
               << vcl_endl;*/

  emit frame_logged(frame_number);
}
 
//: Update the progress bar when a frame is saved.
void ncm_qcapture_gui::onFrameSaved(int frame_number)
{
  // Don't update the progress bar if the processor is still processing.
  if (processor_.is_processing())
    return;

  double progress = static_cast<double>(frame_number) /
                    processor_.num_frames_to_save();

  updateProgress(progress);
}
 
//: Respond when all frames have been saved.
void ncm_qcapture_gui::onSavingFinished()
{
	emit finalize_sequence();
}


void ncm_qcapture_gui::onFinalizeSequence()
{
	//Get alias to current sequence
	ncm_image_sequence &sequence = subject_.live_session()->image_sequences()->back();

	//Only finish the sequence when we've labelled (and thus saved) the sequence and saved all the frames
	if (saver_.is_saving() || sequence.sequence_id() < 0)
		return;

	//try and sync the sequence name (if it's alreday synced or the name is empty nothing will happen)
	sync_sequence_name(sequence);

	//try and recompute motor positions
	interpolate_motor_positions(sequence);
	
	vcl_ofstream ofs;
	ofs.open(sequence.images_full_dir().toStdString() + "/" + sequence.images_datafile().toStdString());
	sequence.t_write(ofs);
	ofs.close();

	//Clear the motor paths
	if (motors_are_live_)
	{
		for ( vcl_vector<QGraphicsRectItem*>::iterator itr = mosaic_frame_positions_.begin(), 
			end = mosaic_frame_positions_.end(); itr != end; ++itr )
	
		 motors_scene_.removeItem(*itr);
		
		mosaic_frame_positions_.clear();
	}

	//Clear the mosaic scene
	mosaic_scene_.clear();

	//Hide the saving frames notice
	scene_.hide_saving_msg();

	//Now we can release the GUI
	activity_state_ = activity_Session;
	update_session_controls();
	updateProgress(0);
}
 

//void ncm_qcapture_gui::correct_artefacts()
//{
  //  vil_image_view<int> vxl_image_int;
  //  vil_convert_cast(vil_plane(vxl_image, 0), vxl_image_int);

  //  vil_math_image_difference(vxl_image_int, artefact_image_, vxl_image_int);

  //  vil_math_truncate_range(vxl_image_int, 0, 255);
  //  vil_convert_cast(vxl_image_int, vil_plane(vxl_image, 0));

  //  qcore_convert_image(*qimage, vil_plane(vxl_image, 0));
//}

//void ncm_qcapture_gui::stabilize_image()
//{
 //   if (!aligner_.is_ready())
 //     aligner_.set_destination(vxl_image);
 //   aligner_.set_source(vxl_image);

 //   int di = 0, dj = 0;
 //   aligner_.align_src_to_dest(di, dj);
 //   aligner_.swap_src_with_dest();

 //   di_cumul += di;
 //   dj_cumul += dj;
 //   scene_.move_pixmap_to(di_cumul, dj_cumul);
//}

//: Redraw items on canvas
void ncm_qcapture_gui::redraw_scene()
{
  QImage* qimage = ncm_qcapture_data_manager::Instance()->get_qimage();

	/*
  const bool is_correcting_artefacts = 
      (ui.correctImageCheckbox->isChecked()) &&
      (artefact_image_.ni() != 0);

  const bool is_stabilizing = (ui.stabilizeCheckbox->isChecked());

  //vil_image_view<vxl_byte> vxl_image;
  //qcore_convert_image(vxl_image, *qimage);

  //if (is_correcting_artefacts)
  //  correct_artefacts();*/
  
  scene_.set_image( qimage );

 // static int di_cumul = 0, dj_cumul = 0;

 // if (is_stabilizing)
 //   stabilize_image();
 // else
 // {
 //   di_cumul = 0;
 //   dj_cumul = 0;
 //   scene_.move_pixmap_to(0, 0);
 // }

  // Don't refit the view every time - it causes minute differences in the 
  // size of the bounding rectangle, which subsequently induce mouseMove events
  // (There has to be a better way than this...)
  static bool first_time = true;
  if (first_time && qimage != NULL && !qimage->isNull())
  {
    ui.graphicsView->fitInView(scene_.sceneRect(), 
                               Qt::KeepAspectRatio);
    first_time = false;
  }

	if (motors_are_live_)
	{
		//Update motor positons view (we'll prob move this to a new function)
		float x = apt_.X()->position();
		float y = apt_.Y()->position();

		if (apt_.X()->reversed())
			x = -x;

		if (apt_.Y()->reversed())
			y = -y;

		current_frame_position_->setRect(x-1.0, y-0.75, 2.0, 1.5);

		if (x < -12.3)
			scene_.show_x_left_limit();
		else
			scene_.hide_x_left_limit();

		if (x > 12.3)
			scene_.show_x_right_limit();
		else
			scene_.hide_x_right_limit();

		if (y < -12.3)
			scene_.show_y_top_limit();
		else
			scene_.hide_y_top_limit();

		if (y > 12.3)
			scene_.show_y_bottom_limit();
		else
			scene_.hide_y_bottom_limit();

		motors_scene_.update();

		float z = apt_.Z()->position();
		if (apt_.Z()->reversed())
			z = -z;

		ui.zposSlider->setValue(10*z);

		if (z < -12.3)
			scene_.show_z_lower_limit();
		else
			scene_.hide_z_lower_limit();

		if (z > 12.3)
			scene_.show_z_upper_limit();
		else
			scene_.hide_z_upper_limit();
	}
	scene_.update();

  //if (autoFocus_state_ != autoFocus_Inactive)
  //{
  //  const QString s = QString::number(processor_.image_sharpness());
  //  ui.autofocusTextEdit->appendPlainText(s);
  //}

  // Compute, display and store the sharpness of the image.
  QString msg;
  if (apt_.n_controllers() > 0)
  {
    msg += "X: " + QString::number(apt_.X()->absolute_position()) + " ";
    msg += "Y: " + QString::number(apt_.Y()->absolute_position()) + " ";
    msg += "Zabs: " + QString::number(apt_.Z()->absolute_position()) + " ";
    msg += "Zrel: " + QString::number(apt_.Z()->position()) + " ";
  }
  msg += " Sh: " + QString::number(processor_.image_sharpness());

	msg += "              AutoFocus: ";
  switch (autoFocus_state_)
  {
    case autoFocus_Inactive: msg += " Inactive"; break;
    case autoFocus_GetDirectionForward: msg += " GetDirectionForward"; break;
    case autoFocus_GetDirectionReverse: msg += " GetDirectionReverse"; break;
    case autoFocus_FindPeak: msg += " FindPeak"; break;
    case autoFocus_GoToPeak: msg += " GoToPeak"; break;
    case autoFocus_Interrupt: msg += " Interrupt"; break;
    case autoFocus_Error: msg += " Error"; break;
    default: msg += " ???"; break;
  }

	msg += "              Session: ";
	switch (activity_state_)
	{
		case activity_Login: break;
		case activity_Home: break;
		case activity_Session: msg +=  "Ready to record"; break;
		case activity_Recording: msg +=  "Recording"; break;
		case activity_RecordingPaused: msg +=  "Recording paused"; break;
		case activity_Saving: msg +=  "Saving frames"; break;
		case activity_Reviewing: break;
		default: assert(false); //should never get here
	}

  updateStatus(msg);
}
 
//:
void ncm_qcapture_gui::draw_mosaic()
{
  //: Return a pixmap object from the input image, scaled in size by scale, and 
//  with grey levels adjusted by gain and offset.

	//Clear any existing stuff from the mosaic view
	mosaic_scene_.clear();

  double scale = 1.0;
  double gain = 1.0;
  double offset = 1.0;

  vil_image_view<vxl_byte> vxl_image;

  if ( false) //(vcl_abs(gain - 1.0) > 1e-6) || (vcl_abs(offset - 0.0) > 1e-6) )
  {
    // Scale and offset a copy of the input image, since a vil_image_view acts
    // as a pointer and will change the aligner's copy of the image.
    vxl_image.deep_copy(stitcher_->mosaic());
    vil_math_scale_and_offset_values(vxl_image, gain, offset);
  }
  else
    vxl_image = stitcher_->mosaic();

  QImage qt_image;
  qcore_convert_image(qt_image, vxl_image);
  
  QPixmap pixmap;
  int scaled_width = scale * qt_image.width();
  if (scaled_width == qt_image.width())
    pixmap.convertFromImage(qt_image);
  else
    pixmap.convertFromImage(qt_image.scaledToWidth(scaled_width));

	QGraphicsPixmapItem* p = mosaic_scene_.addPixmap(pixmap);
  mosaic_scene_.setSceneRect(p->boundingRect());

  ui.mosaicView->fitInView(mosaic_scene_.sceneRect(), Qt::KeepAspectRatio	);

}

void ncm_qcapture_gui::onProcessingFinished(int n_processed)
{	
	
	//Get reference to sequence		
	ncm_image_sequence &sequence = subject_.live_session()->image_sequences()->back();

	//Finalise the current sequence
	sequence.finalise();
	sequence.set_num_frames(n_processed);

	scene_.show_mosaic_msg();

	//Update the session display but don't change the activity state (and thus release the
	//GUI until svaing is finsihed)
	update_session_display();

	//We're now waiting for onMosaicFinished() before we get the user to tag this session
}
 
//: Update the scene at regular intervals determined by a timer.
//  (TODO: This is unnecessarily messy and should be tidied up soon.)
void ncm_qcapture_gui::onAptTimer()
{
	/*
  ui.aptXPosition->setValue(apt_.X()->position());
  ui.aptYPosition->setValue(apt_.Y()->position());
  ui.aptZPosition->setValue(apt_.Z()->position());
	*/
	
  // It's ok to call the move_at_velocity() function directly here because it
  // always returns immediately, unlike move_by().

  // Set up the velocity look-up table
  const unsigned array_width = 16; // multiples of two are good
  const unsigned array_size = 2*array_width + 1;
  static int velocity_index = -1;
  static vnl_vector<float> velocity_array(array_size, 0.0);

  const double sinusoid_amplitude = 1.0;

  // There is an interdependence here between the frequency of the wave 
  // determined by the frequency of the timer and array_width) and the 
  // amplitude of the sinusoid (sinusoid_amplitude).
  // These determine the maximum rate of velocity change, which must be 
  // within the acceleration limits of the motorized platform.

  if (velocity_index == -1)
  {
    velocity_index = 0;

    for (unsigned i = 0; i < array_size; ++i)
    {
      const float x = (static_cast<float>(i) / array_size) * (2 * M_PI);
      velocity_array[i] = sinusoid_amplitude * vcl_sin(x);
    }
  }

  static vnl_vector<float> sharpness_array(array_size, 0.0);
  static vnl_vector<float> position_array(array_size, 0.0);
  static vnl_vector<float> weighted_dz_array(array_size, 0.0);
  static bool array_is_filled = false;

  if (motors_are_live_)
  {
    if (velocityX_ != 0)
      apt_.X()->move_at_relative_velocity(velocityX_);


    if (velocityY_ != 0)
      apt_.Y()->move_at_relative_velocity(velocityY_);


    const bool is_panning = (velocityX_ != 0) ||
                            (velocityY_ != 0);

    if (velocityZ_ != 0)
    {
      // Interrupt any autofocus process that's underway.
      autoFocus(autoFocus_Interrupt);

      apt_.Z()->move_at_velocity(velocityZ_);
    }
    else if ( (autoFocus_state_ == autoFocus_Inactive) && 
              is_panning )
    {
      // Update the sharpness array
      sharpness_array[velocity_index] = processor_.image_sharpness();
      position_array[velocity_index] = apt_.Z()->position();

      // Estimate the direction in which the z-axis should move.
      double weighted_dz = 0.0;
      if (array_is_filled)
      {
        double sum_sharpness = 0.0;
        double sum_position = 0.0;
        for (unsigned i = 0; i < array_size; ++i)
        {
          sum_sharpness += sharpness_array[i];
          sum_position += position_array[i];

          weighted_dz += position_array[i] * sharpness_array[i];
        }
        // Weighted position
        weighted_dz /= sum_sharpness;

        // Weighted displacement with respect to mean.
        weighted_dz -= (sum_position / array_size);
        weighted_dz_array[velocity_index] = weighted_dz;
      }

      double mean_mag_dz = 0.0;
      for (unsigned i = 0; i < array_size; ++i)
        mean_mag_dz += vcl_abs(weighted_dz_array[i]);
      mean_mag_dz /= array_size;

			if ( preferences_.use_autofocus_ && ui.focusDuringPan->isChecked() )
      {
        float velocityZ = velocity_array[velocity_index];
				double modulated_weight = preferences_.modulate_weight_ * weighted_dz;

        if (array_is_filled)
          velocityZ += modulated_weight;

        // Scale the sinusoid speed by the panning speed.
        const double speed = vcl_sqrt(velocityX_*velocityX_ + 
                                      velocityY_*velocityY_);
        velocityZ *= speed * 5;
        //velocityZ *= mean_mag_dz;

        if (velocityZ != 0)
          apt_.Z()->move_at_velocity(velocityZ);
        else
          apt_.Z()->stop();

        if ( focus_log_ && array_is_filled )
        {
          focus_log_ << apt_.Z()->position() << '\t'
                     << sharpness_array[velocity_index] << '\t'
                     << velocity_array[velocity_index] << '\t'
                     << modulated_weight << '\t'
                     << speed << '\t'
                     << velocityZ << '\t'
                     << mean_mag_dz << '\n';
        }
      }
      else // don't modulate Z axis (zoom)
        apt_.Z()->stop();

      // Increment the velocity counter.
      ++velocity_index;
      if (velocity_index == array_size)
      {
        velocity_index = 0;
        array_is_filled = true;
      }
		}
  }
}
 
//: Change the controller for each axis
void ncm_qcapture_gui::on_aptXCombobox_activated(int index)
{
  // If combobox is empty (i.e. the list has just been cleared) then do nothing.
  switch (index)
  {
    case -1: // List is empty
      ui.aptXHome->setEnabled(false);
      return;
    case 0: // 'None' chosen
      ui.aptXHome->setEnabled(false);
      apt_.set_X_id(-1);
      break;
    default: // 
      ui.aptXHome->setEnabled(true);
      apt_.set_X_id(ui.aptXCombobox->currentText().toLong());
			QSqlQuery query("UPDATE motor_ids SET AptXId = " + QString::number(apt_.X()->id()) + " WHERE motor_id=1", database_);
  }

  update_apt_comboboxes();
}
void ncm_qcapture_gui::on_aptYCombobox_activated(int index)
{
  // If combobox is empty (i.e. the list has just been cleared) then do nothing.
  switch (index)
  {
    case -1:
      ui.aptYHome->setEnabled(false);
      return;
    case 0:
      ui.aptYHome->setEnabled(false);
      apt_.set_Y_id(-1);
      break;
    default:
      ui.aptYHome->setEnabled(true);
      apt_.set_Y_id(ui.aptYCombobox->currentText().toLong());
      QSqlQuery query("UPDATE motor_ids SET AptYId = " + QString::number(apt_.Y()->id()) + " WHERE motor_id=1", database_);
  }

  update_apt_comboboxes();
}
void ncm_qcapture_gui::on_aptZCombobox_activated(int index)
{
  // If combobox is empty (i.e. the list has just been cleared) then do nothing.
  switch (index)
  {
    case -1:
      ui.aptZHome->setEnabled(false);
      return;
    case 0:
      ui.aptZHome->setEnabled(false);
      apt_.set_Z_id(-1);
      break;
    default:
      ui.aptZHome->setEnabled(true);
      apt_.set_Z_id(ui.aptZCombobox->currentText().toLong());
      QSqlQuery query("UPDATE motor_ids SET AptZId = " + QString::number(apt_.Z()->id()) + " WHERE motor_id=1", database_);
  }

  update_apt_comboboxes();
}

//
//: Reverse directions of relative movements.
void ncm_qcapture_gui::on_aptXReverse_toggled(bool checked)
{
  apt_.X()->set_reverse(checked);
	QSqlQuery query("UPDATE motor_ids SET XReversed = " + QString::number(checked) + " WHERE motor_id=1", database_);
}
void ncm_qcapture_gui::on_aptYReverse_toggled(bool checked)
{
  apt_.Y()->set_reverse(checked);
  QSqlQuery query("UPDATE motor_ids SET YReversed = " + QString::number(checked) + " WHERE motor_id=1", database_);
}
void ncm_qcapture_gui::on_aptZReverse_toggled(bool checked)
{
  apt_.Z()->set_reverse(checked);
  QSqlQuery query("UPDATE motor_ids SET ZReversed = " + QString::number(checked) + " WHERE motor_id=1", database_);
}

//
//: Home the motors
void ncm_qcapture_gui::on_aptXHome_clicked()
{
  home_x();
}
void ncm_qcapture_gui::on_aptYHome_clicked()
{
  home_y();
}
void ncm_qcapture_gui::on_aptZHome_clicked()
{
  home_z();
}
 
void ncm_qcapture_gui::onSceneJoystickMoved(double x, double y)
{
  // 1 normalized unit half the frame height (= 240px on DMK cameras)
  //
	float scale = preferences_.joystick_max_;

  velocityX_ = scale * static_cast<float>(x);
  velocityY_ = scale * static_cast<float>(y);

}

void ncm_qcapture_gui::onSceneJoystickReleased()
{
  velocityX_ = 0.0f;
  apt_.X()->stop();

  velocityY_ = 0.0f;
  apt_.Y()->stop();

	if ( preferences_.use_autofocus_ && 
		ui.focusOnRelease->isChecked() )
    autoFocus(autoFocus_GetDirectionForward);
}

void ncm_qcapture_gui::onSceneDoubleClicked(
  double x, 
  double y, 
  Qt::MouseButton button)
{
  if (button == Qt::LeftButton)
  {
    set_apt_targetX(x);
    set_apt_targetY(y);
  }
}

//: React to mouse wheel event on the graphics view
void ncm_qcapture_gui::on_graphicsView_zoomed(int delta)
{
  // Apply a sensible relative movement in Z.
  // Typically, one 'click' of the mouse wheel will generate a delta of 120.
  // Scale down such that each wheel 'click' moves by value selected in preferences
	//vcl_cout << "Set APT z target " << static_cast<double>(delta)/12000 << vcl_endl;
	double one_mm = static_cast<double>(delta) / 120;
	set_apt_targetZ( one_mm * preferences_.z_scroll_dist_);
}

//:Event filter so we can process arrow keys to change the Z-axis of the motor
//: these will only be processed if the camera scene is in focus
//: (this may require clicking on the scene if X/Y motor controls haven't be used)
bool ncm_qcapture_gui::eventFilter(QObject *object, QEvent *event)
{
 if (event->type() == QEvent::KeyPress)
 {

	if (ui.graphicsView->hasFocus())
	{
		QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
		if (keyEvent->key() == Qt::Key_Up)
			set_apt_targetZ(preferences_.z_arrow_dist_);
		else if (keyEvent->key() == Qt::Key_Down)
			set_apt_targetZ(-preferences_.z_arrow_dist_);

		return true;
	}
 }
 return false;
 
}

//Review slots
void ncm_qcapture_gui::onSessionSelectionChanged()
{
	//Get the subject record from the sql model, then load this record into the current subject object
	QModelIndex filterIndex = ui.sessionListView->selectionModel()->selectedIndexes().first();
	QModelIndex srcIndex = sessionFilterModel_->mapToSource(filterIndex);
	
	QRect ri = ui.sessionListView->visualRect(filterIndex);

	//QSqlQuery query(database_);
	ncm_image_session session;
	session.sql_read( sessionModel_->record(srcIndex.row()));//, query

	//Find the index of this session in the subject's sessions list
	//int idx = subject_.find_session(session.se

	//populate the session list view with sessions from this subject
	QString sequence_query_str = "SELECT * FROM image_sequence WHERE image_session_id = " + QString::number(session.session_id());

	//Query the database and get the column index of the subject name
	sequenceModel_->setQuery(sequence_query_str, database_);
	int fieldNo = sequenceModel_->record().indexOf("sequence_name");
	
	//Link the SQL model to the filter model so we can filter on subject name
	sequenceFilterModel_->setSourceModel( sequenceModel_ );
	sequenceFilterModel_->setFilterKeyColumn(fieldNo);

	//Link the filter model to the list view to display the filtered subjects
	ui.sequenceListView->setModel(sequenceFilterModel_);
	ui.sequenceListView->setModelColumn(fieldNo);
	ui.sequenceListView->setSelectionMode(QAbstractItemView::SingleSelection);
	ui.sequenceListView->setSelectionBehavior(QAbstractItemView::SelectRows);
  ui.sequenceListView->show();

	//Connect up the sequence view selection model so we process selection changes
	connect(ui.sequenceListView->selectionModel(), SIGNAL( selectionChanged (QItemSelection, QItemSelection)),
				this, SLOT(onSequenceSelectionChanged()),
				Qt::UniqueConnection);

	connect(this, SIGNAL( review_sequence_selected(QModelIndex)),
				ui.sequenceListView, SIGNAL(pressed(QModelIndex) ),
				Qt::UniqueConnection);

	//Set with the previous or next sequence buttons are enabled
	ui.prevSessionButton->setEnabled(filterIndex.row() > 0);
	ui.nextSessionButton->setEnabled(filterIndex.row() < (sessionFilterModel_->rowCount()-1));

	if (sequenceFilterModel_->rowCount() > 0)
	{
		const QModelIndex i = sequenceFilterModel_->index(0,0);
		ui.sequenceListView->setCurrentIndex(i);
	}
}

void ncm_qcapture_gui::onSequenceSelectionChanged()
{
	if (ui.sequenceListView->selectionModel()->hasSelection())
	{
		//Get the subject record from the sql model, then load this record into the current subject object
		QModelIndex filterIndex = ui.sequenceListView->selectionModel()->selectedIndexes().first();
		QModelIndex srcIndex = sequenceFilterModel_->mapToSource(filterIndex);

		ncm_image_sequence sequence;
		sequence.sql_read( sequenceModel_->record(srcIndex.row()));

		imageProcessor_.clear();
		review_scene_.clear();

		vcl_string mosaic_name;
		if (!sequence.full_mosaic().isEmpty())
		{
			mosaic_name = sequence.images_full_dir().toStdString() + "/" + sequence.full_mosaic().toStdString();
		}
		else if (!sequence.preview_mosaic().isEmpty())
		{
			mosaic_name = sequence.images_full_dir().toStdString() + "/" + sequence.preview_mosaic().toStdString();
		}
		if (!mosaic_name.empty())
		{	
			// Tell the image processor which file to load then load it.
			imageProcessor_.set_filename(mosaic_name);
			imageProcessor_.load_image();
			review_scene_.fit_to_image();
			ui.reviewView->fit_both();
			ui.reviewView->setEnabled(true);
			ui.zoomGroup->setEnabled(true);
			ui.imageGroup->setEnabled(true);
			no_mosaic_msg_->hide();
		}
		else
		{
			QRectF r = no_mosaic_msg_->boundingRect();
			QRectF sr = review_scene_.sceneRect();
			if (preferences_.rotate_frames_)
				no_mosaic_msg_->setPos(sr.center().x()+r.center().x(),sr.center().y()+r.center().y());
			else
				no_mosaic_msg_->setPos(sr.center().x()-r.center().x(),sr.center().y()-r.center().y());
			no_mosaic_msg_->show();
			ui.reviewView->setEnabled(false);
			ui.zoomGroup->setEnabled(false);
			ui.imageGroup->setEnabled(false);
		}
		QString hand = QString(ncm_hand::toLetter(sequence.hand()).c_str());
		QString digit = QString::number(sequence.digit());
		QString nframes = QString::number(sequence.num_frames());
		QString time = sequence.time_finished().toString("hh:mm:ss");
		QString sequence_string;
		sequence_string +=  ("Hand: " + hand + ", digit: " + digit + "\n");
		sequence_string += (nframes + " frames, " + time);
		ui.sequenceDetailsText->setText(sequence_string);

		//Set with the previous or next sequence buttons are enabled
		ui.prevSequenceButton->setEnabled(filterIndex.row() > 0);
		ui.nextSequenceButton->setEnabled(filterIndex.row() < (sequenceFilterModel_->rowCount()-1));

		vcl_cout << "Sequence selection changed, new sequence = " << sequence.sequence_name().toStdString() << vcl_endl;
	}
}

//Move between sessions in review mode
void ncm_qcapture_gui::on_prevSessionButton_clicked()
{
	const QModelIndex curr = ui.sessionListView->currentIndex();
	const QModelIndex prev = sessionFilterModel_->index(curr.row()-1,0);
	ui.sessionListView->setCurrentIndex(prev);
}
void ncm_qcapture_gui::on_nextSessionButton_clicked()
{
	const QModelIndex curr = ui.sessionListView->currentIndex();
	const QModelIndex next = sessionFilterModel_->index(curr.row()+1,0);
	ui.sessionListView->setCurrentIndex(next);
}

//Move between sequences in review mode
void ncm_qcapture_gui::on_prevSequenceButton_clicked()
{
	const QModelIndex curr = ui.sequenceListView->currentIndex();
	const QModelIndex prev = sequenceFilterModel_->index(curr.row()-1,0);
	ui.sequenceListView->setCurrentIndex(prev);
}
void ncm_qcapture_gui::on_nextSequenceButton_clicked()
{
	const QModelIndex curr = ui.sequenceListView->currentIndex();
	const QModelIndex next = sequenceFilterModel_->index(curr.row()+1,0);
	ui.sequenceListView->setCurrentIndex(next);
}

// Connect to motors.
void ncm_qcapture_gui::on_connectToMotors_clicked()
{
  if (apt_timer_.isActive())
  {
    vcl_cout << "Disconnecting motors" << vcl_endl;

    apt_.X()->stop();
    apt_.Y()->stop();
    apt_.Z()->stop();

    apt_timer_.stop();
    focus_log_.close();
  }
  else
  {
    vcl_cout << "Connecting motors" << vcl_endl;
    initialize_apt_timer();

    focus_log_.open("u:/projects/nailfold/tmp/focus_log.txt");
  }
}

// Initiate the calibration motion.
void ncm_qcapture_gui::on_calibrate_clicked()
{
  handleCalibrateState();
}

//Autofocus callback
void ncm_qcapture_gui::on_autoFocus_clicked()
{
  if (!motors_are_live_)
    return;

  //ui.plainTextEdit->setPlainText("");

  // Begin the autofocus procedure by searching in the forward direction for
  // a measurable increase in sharpness.
  autoFocus(autoFocus_GetDirectionForward);
}
 
//: React to motors completing a move.
void ncm_qcapture_gui::onAptMoveComplete(ncm_qapt_server::apt_axis)
{
  if (autoFocus_state_ != autoFocus_Inactive)
    handleAutoFocusState();
  else if (calibrate_state_ != calibrate_Inactive)
    handleCalibrateState();
}
 
//: Handle changes in autofocus_state_.
void ncm_qcapture_gui::handleAutoFocusState()
{
  switch (autoFocus_state_)
  {
    case autoFocus_Inactive:
      // What are we doing here?
      break;

    case autoFocus_GetDirectionForward:
      // Reverse the direction of search because we've gone far enough that we
      // would have found an improvement had we been going in the right direction.
      apt_.set_max_displacement(ncm_qapt_server::AptZ, 1e9);
      autoFocus(autoFocus_GetDirectionReverse);
      break;

    case autoFocus_GetDirectionReverse:
      // Go to error state since we can't find a direction that generates a 
      // measurable improvement in sharpness. This is likely to happen if 
      // autofocus is started too far from the true position (in the flat region
      // of the function).
      apt_.set_max_displacement(ncm_qapt_server::AptZ, 1e9);
      autoFocus(autoFocus_Error);
      break;

    case autoFocus_FindPeak:
      break;

    case autoFocus_GoToPeak:
      autoFocus(autoFocus_Inactive);
      break;

    case autoFocus_Error:
      break;

    default:
      break;
  }
}
 
//: Handle changes in calibrate_state_.
void ncm_qcapture_gui::handleCalibrateState()
{
  const float move_distance = 0.15f;

  switch (calibrate_state_)
  {
  case calibrate_InvalidFirst:
  case calibrate_Error:
    // Do nothing (for now)
    break;

  case calibrate_Inactive:
    // Go up.
    calibrate_state_ = calibrate_Up;
    emit apt_move_by(ncm_qapt_server::AptY, move_distance, false);
    break;

  case calibrate_Up:
    // Go right.
    calibrate_state_ = calibrate_Right;
    emit apt_move_by(ncm_qapt_server::AptX, move_distance, false);
    break;

  case calibrate_Right:
    // Go down.
    calibrate_state_ = calibrate_Down;
    emit apt_move_by(ncm_qapt_server::AptY, -move_distance, false);
    break;

  case calibrate_Down:
    // Go left.
    calibrate_state_ = calibrate_Left;
    emit apt_move_by(ncm_qapt_server::AptX, -move_distance, false);
    break;

  case calibrate_Left:
    // Reset to Inactive state.
    calibrate_state_ = calibrate_Inactive;
    break;

  default:
    // Huh?
    break;
  }
}
 
//:
void ncm_qcapture_gui::onLowerSharpnessThresholdCrossed()
{
  vcl_stringstream ss;
  ss << "LowerThresholdCrossed" << vcl_endl;
  ui.autofocusTextEdit->appendPlainText(ss.str().c_str());

  switch (autoFocus_state_)
  {
    case autoFocus_GetDirectionForward:
      // Reverse the direction of search because we're going the wrong way.
      autoFocus(autoFocus_GetDirectionReverse);
      break;

    case autoFocus_FindPeak:
      // Stop looking for the peak once the lower threshold has been crossed, as
      // we're now 'over the hump' and getting worse again.
      autoFocus(autoFocus_GoToPeak);
      break;

    case autoFocus_Inactive:
      // Fall through:
      //    Limits should not be active at this stage of the autofocus.
      //    Something's gone wrong.
    case autoFocus_GetDirectionReverse:
      // Fall through:
      //    Go to error state as this shouldn't, in theory, happen: we're in this
      //    state because going in the opposite direction decreased sharpness, so 
      //    this direction *must* increase sharpness.
      vcl_cout << "Error: Lower threshold triggered during ReverseSearch"
               << vcl_endl;
      autoFocus(autoFocus_Error);
      break;

    case autoFocus_GoToPeak:
      // Fall through:
      //    Limits should not be active at this stage of the autofocus.
      //    Something's gone wrong.
      vcl_cout << "Error: Lower threshold triggered during GoToPeak"
               << vcl_endl;
      autoFocus(autoFocus_Error);
      break;

    case autoFocus_Error:
      autoFocus(autoFocus_Error);
      break;

    default:
      break;
  }
}

void ncm_qcapture_gui::onUpperSharpnessThresholdCrossed()
{
  vcl_stringstream ss;
  ss << "UpperThresholdCrossed" << vcl_endl;
  ui.autofocusTextEdit->appendPlainText(ss.str().c_str());

  switch (autoFocus_state_)
  {
    case autoFocus_GetDirectionForward:
      // Fall through
    case autoFocus_GetDirectionReverse:
      // Find the peak, since we know we're now going in the right direction.
      autoFocus(autoFocus_FindPeak);
      break;

    case autoFocus_FindPeak:
    {
      // Update the thresholds used for peak detection.
      autoFocus(autoFocus_FindPeak);
      break;
    }

    case autoFocus_Inactive:
      // Limits should not be active at this stage of the autofocus.
      // Something's gone wrong.
      // (Fall through)
    case autoFocus_GoToPeak:
      // Limits should not be active at this stage of the autofocus.
      // Something's gone wrong.
      // (Fall through)
      vcl_cout << "Error: Upper threshold triggered during GoToPeak"
               << vcl_endl;
      autoFocus(autoFocus_Error);
      break;

    case autoFocus_Error:
      autoFocus(autoFocus_Error);
      break;

    default:
      break;
  }
}

//Align frames for mosaic
void ncm_qcapture_gui::onFrameToAlign(int frame_num)
{
	if (motors_are_live_)
	{
		float w = 640 / preferences_.pixels_per_motor_mm_;//scene_.width()
		float h = 480 / preferences_.pixels_per_motor_mm_;//scene_.height()
	
		//Update motor positons view (we'll prob move this to a new function)
		float x = apt_.X()->position();
		float y = apt_.Y()->position();

		if (apt_.X()->reversed())
			x = -x;

		if (apt_.Y()->reversed())
			y = -y;

		QGraphicsRectItem* r = motors_scene_.addRect(x - w/2, y - h/2, w, h,
			QPen(QColor(0,0,255)), QBrush());
		mosaic_frame_positions_.push_back(r);
	}
		
	//vcl_cout << "frame to align: " << frame_num << vcl_endl;
	emit frame_to_align();
}

void ncm_qcapture_gui::onFramesAligned(int di, int dj, double scale, double offset, double mse)
{
	//Take frame from head of align queue - this should always be the 'dest' frame
	//i.e. the displacements di and dj are where the next frame should be placed...
	QSharedPointer< ncm_video_frame > vid_frame = ncm_qcapture_data_manager::Instance()->remove_from_align_queue();
	
	//Take a deep copy of the image in the vid frame so we can transform it
	vil_image_view<vxl_byte> vxl_image;
	vxl_image.deep_copy(*(vid_frame->frame()));

	if ( (vcl_abs(scale - 1.0) > 1e-6) ||
       (vcl_abs(offset - 0.0) > 1e-6) )
  {
    // Scale and offset the image
    vil_math_scale_and_offset_values(vxl_image, scale, offset);
  }

  QImage qt_image;
  qcore_convert_image(qt_image, vxl_image);
  
  QPixmap pixmap;
  int scaled_width = scale * qt_image.width();
  if (scaled_width == qt_image.width())
    pixmap.convertFromImage(qt_image);
  else
    pixmap.convertFromImage(qt_image.scaledToWidth(scaled_width));

	QGraphicsPixmapItem* item = mosaic_scene_.addPixmap(pixmap);
  item->setPos(current_mosaic_displacement_);

	//Increment the displacement for next time
	current_mosaic_displacement_ += QPoint(di, dj);	
}

void ncm_qcapture_gui::onAlignmentFinished()
{
	vcl_cout << "alignment finished!" << vcl_endl;

	//In theory there should always be one more frame to pop from the aligner queue, so clear this now...
	//We could draw it if we want but it'll be replaced by the full mosaic in a second anyway so not much point
	assert(ncm_qcapture_data_manager::Instance()->align_queue()->count()==1);
	ncm_qcapture_data_manager::Instance()->remove_from_align_queue();

	//Instantiate the stitcher, connect up sticther signals and move to the aligner thread
	stitcher_ = new ncm_qcapture_mosaic_maker(aligner_.ni(), aligner_.nj());
	stitcher_->set_max_pixels((unsigned)preferences_.max_mosaic_size_*1e6);

	// this -> stitcher_
  QObject::connect( this,  SIGNAL(make_mosaic(vcl_vector<QPoint>)),
                    stitcher_, SLOT(makeMosaicFrom(vcl_vector<QPoint>)) ); 

  // stitcher_ -> this
  QObject::connect( stitcher_, SIGNAL(mosaicUpdated(int)),
                    this, SLOT(onMosaicUpdated(int)) );
  QObject::connect( stitcher_, SIGNAL(mosaicFinished(bool)),
                    this, SLOT(onMosaicFinished(bool)) );

	//Move to thread
	stitcher_->moveToThread(&aligner_thread_);

	//Get aligner displacements and start mosaicing...
	vcl_vector<QPoint> displacements = aligner_.displacements();
	emit make_mosaic(displacements);
}

void ncm_qcapture_gui::onMosaicUpdated(int frame_num)
{
	//vcl_cout << "Frame " << frame_num << " added to mosaic" << vcl_endl;
}

void ncm_qcapture_gui::onMosaicFinished(bool success)
{
	//Can disconnect the mosaic signals
	QObject::disconnect( stitcher_, SIGNAL(mosaicUpdated(int)),
                    this, SLOT(onMosaicUpdated(int)) );
  QObject::disconnect( stitcher_, SIGNAL(mosaicFinished(bool)),
                    this, SLOT(onMosaicFinished(bool)) );

	//Switch message on scene to saving frames (if saving already finished this
	//will instantly be removed in finalize sequence)
	scene_.hide_mosaic_msg();
	scene_.show_saving_msg();

	//Get reference to sequence		
	ncm_image_sequence &sequence = subject_.live_session()->image_sequences()->back();

	if (success)
	{
		vcl_cout << "Mosaic finished!" << vcl_endl;
		vcl_cout << "Mosaic size, w = " << stitcher_->mosaic().ni() << ", h = " << stitcher_->mosaic().nj() << vcl_endl;

		//Draw this mosaic
		draw_mosaic();

		//Save the mosaic
		sequence.set_preview_mosaic("preview_mosaic.png"); //TO DO make this a user preference...
		vcl_string mosaic_name = sequence.images_full_dir().toStdString() + "/" + sequence.preview_mosaic().toStdString();
		vil_save(stitcher_->mosaic(), mosaic_name.c_str() );
		
	}
	else
	{
		QMessageBox msgBox;
		msgBox.setIcon(QMessageBox::Warning);
		msgBox.setText(
			"Mosaic too large to build during live imaging. \n"
			"All image frames will be saved and the mosaic can \n"
			"be completed offline");
		msgBox.exec();
	}

	//Reset the stitcher to NULL
	delete stitcher_;
  stitcher_ = NULL;

	//Ask the user what hand and digit were imaged (and update the sequence name)
	ncm_qcapture_end_sequence sequence_ui(sequence, study_name_, subject_.name_in_study(), this);
	sequence_ui.exec();

	//Write the sequence to the dabatase
	QSqlQuery query(database_);
	sequence.sql_write(query);

	

	emit finalize_sequence();
}
 
//
//  Private methods
//
 
//: Create the status bar's labels and progress bar
void ncm_qcapture_gui::initializeStatusbar()
{
  statusBar_main_.setFrameStyle(QFrame::Panel & QFrame::Sunken);
  statusBar_main_.setLineWidth(1);
  statusBar_main_.setText("");
  statusBar_main_.setContentsMargins(4, 0, 4, 0);
  statusBar()->addWidget(&statusBar_main_, 1);

  statusBar_progress_.setFixedWidth(160);
  statusBar()->addWidget(&statusBar_progress_);
}

 
//: Update status bar.
//  A progress value in [0..1] sets the progress bar value.
//  progress < 0 makes no change.
void ncm_qcapture_gui::updateStatus(
  const QString& status, 
  double progress /* = -1 */)
{
  statusBar_main_.setText(status);
  updateProgress(progress);
  statusBar()->update();
}
 
//: Update progress bar only.
//  A progress value in [0..1] sets the progress bar value.
//  progress < 0 makes no change.
void ncm_qcapture_gui::updateProgress(
  double progress)
{
  if ((0.0 <= progress) && (progress <= 1.0))
  {
    const int progress_min = statusBar_progress_.minimum();
    const int progress_range = statusBar_progress_.maximum() - progress_min;
    statusBar_progress_.setValue(progress_min + progress*progress_range);
  }
}
 
//:
void ncm_qcapture_gui::initialize_aligner_thread()
{
  aligner_thread_.setObjectName("aligner_thread");
  
	//radius for displacing??
	int radius = 40;
	aligner_.set_filter(ncm_frame_aligner::filter_g1d);
  aligner_.set_n_levels(1);
  aligner_.set_di_radius(radius);
  aligner_.set_dj_radius(radius);
	aligner_.setNumToAlign(1e4);
  aligner_.moveToThread(&aligner_thread_);

	// Register the points vector so that it can be used in signals
  qRegisterMetaType< vcl_vector<int> >("vcl_vector<QPoint>");

  aligner_thread_.start();
	
}
 
//: Write sequence properties to a text file.
void ncm_qcapture_gui::write_properties( vcl_ostream& tfs )
{
	tfs << vsl_indent() << "NCM QCapture Sequence Properties" << vcl_endl;
  vsl_indent_inc(tfs);

	tfs << vsl_indent() << "Date: " << QDate::currentDate().toString( "ddd dd MMM yyyy" ).toStdString() << vcl_endl;
	tfs << vsl_indent() << "Time: " << QTime::currentTime().toString( "hh:mm:ss" ).toStdString() << vcl_endl;
	tfs << vcl_endl;
	tfs << vsl_indent() << "Camera Type: " << camera_.get_camera_type() << vcl_endl;
	tfs << vsl_indent() << "Camera ID: " << camera_.get_camera_id() << vcl_endl;
  vsl_indent_inc(tfs);
	tfs << vsl_indent() << "Frame rate: " << camera_.get_current_FPS() << " fps" << vcl_endl;
	tfs << vsl_indent() << "Exposure time: ";
	if ( camera_.is_auto_exposure() )
		tfs << "auto" << vcl_endl;
	else
		tfs << camera_.get_current_exposure() << vcl_endl;

	tfs << vsl_indent() << "Gain: ";
	if ( camera_.is_auto_gain() )
		tfs << "auto" << vcl_endl;
	else
		tfs << camera_.get_current_gain() << vcl_endl;

  tfs << vcl_flush;
}

//-------------------------------------------------------------------------
//Private methods
void ncm_qcapture_gui::log_in()
{
	//We now have a valid user...

	//Should initialize anything generic that wil be used by all of the
	//main functions (capturing, reviewing, reporting)

	
	//If there are tasks specific to a particular function (e.g. connecting
	// and homing the motors for an image session) they should live in the start 
	// function for that session

	//Move to the home screen and continue with the software
	activity_state_ = activity_Home;
	update_session_controls();

}

void ncm_qcapture_gui::log_off()
{
	disconnect_database();
  
	//reset the user to the default
	user_ = ncm_user(1);

	activity_state_ = activity_Login;
	update_session_controls();
}

void ncm_qcapture_gui::start_image_session()
{
	//Create new instance of search for subject GUI
	ncm_qcapture_search search_ui(subject_, study_name_, database_, ncm_qcapture_search::ImageSession, this);
	int result = search_ui.exec();

	if (subject_.ncm_id() == -1)
		//Not subject was selected, just return without doing anything
		return;

	//This all seems to take a little time, so show a busy cursor until the select subject
	//dialog box opens
	setCursor(QCursor(Qt::WaitCursor));

	//and only quit them when the program quits
	initialize_processor_thread();
  initialize_saver_thread();
  initialize_aligner_thread();
	initialize_apt_thread();
	connect_session_signals_to_slots();

	//What happens below is only applicable to running
	//a capture session - e.g. if we're going to reveiw images we don't need
	//to connect and home motors.
	initialize_motors();
	initialize_scene();

	//Link the data managers queue to the camera
	camera_.set_queue( ncm_qcapture_data_manager::Instance()->main_queue() );

  //initialize_saver(); Don't need to do this

	processor_.resize_sharpness_array(preferences_.n_sharpness_);

  //vcl_string artefacts_filename = 
  //    "U:/projects/nailfold/capture/2013_02_21/Left/Digit4/x300/15_01_31/_diff.png";
  //vil_image_view<vxl_byte> artefact_image_byte = 
  //    vil_load(artefacts_filename.c_str());
  //vil_convert_cast(artefact_image_byte, artefact_image_);
  //vil_math_scale_and_offset_values(artefact_image_, 1.0, -128);
  
  initializeStatusbar();
	update_tab_controls();
	
	//Create a new session
	ncm_image_session session(subject_.ncm_id(), user_.user_id());
	session.start();
	QString session_name = session.time_started().toString("yyyy_MM_dd");

	//Write the new session to the database, then add it as the current session of the subject
	QSqlQuery query(database_);
	session.sql_write(query); // This generates a unique session ID
	subject_.set_live_session(session);

	//Set the subject name in the saver object
	saver_.setStudyName(study_name_.toStdString() );
	saver_.setSubjectName( subject_.name_in_study().toStdString() );
	saver_.setSessionName( session_name.toStdString() );

	//Update the subject display and switch to the session tab
	activity_state_ = activity_Session;
	update_subject_display();
	update_session_display();
	update_session_controls();
	ui.stackedWidget->setCurrentWidget(ui.pageSession);
	ui.tabWidget->setCurrentWidget(ui.tabPosition);

	//Load up the new camera (try loading existing camera)
	load_camera(preferences_.auto_load_camera_);

	setCursor(QCursor(Qt::ArrowCursor));
}

void ncm_qcapture_gui::start_offline_processing()
{
	QSqlQuery query(database_);
	bool success = query.exec("SELECT sequence_id FROM sequences_to_process");

	if (success)
	{
		//Load sequence ids into a vector
		vcl_vector<int> sequences_to_process;
		while (query.next())
			sequences_to_process.push_back(query.value(0).toInt());

		int num_seqs_to_process = sequences_to_process.size();
		vcl_cout << "Num sequences to process offline = " << num_seqs_to_process << vcl_endl;

		//Now change the query to select sequences from the databse
		for (size_t i = 0; i < num_seqs_to_process; ++i)
		{
			query.exec("SELECT * FROM image_sequence WHERE image_sequence_id = " + QString::number(sequences_to_process[i]));
			assert(query.first());
			ncm_image_sequence sequence;
			sequence.sql_read(query.record());
			align_full_mosaic(sequence);
		}

		

	}
	else
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;

	 
}

void ncm_qcapture_gui::end_image_session()
{
	bool finalise_session = true;

	//If we're set to delete empty session, check if the session had any sequences
	if (preferences_.delete_empty_sessions_ &&
		subject_.live_session()->image_sequences()->empty())
	{
		QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Question);
    msgBox.setText("No sequences captured for this session.");
    msgBox.setInformativeText("Delete the session?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::Yes);

    const int response = msgBox.exec();

    if (response == QMessageBox::Yes)
		{
			//Delete session from database
			QSqlQuery query(database_);
			subject_.live_session()->sql_delete(query);
			finalise_session = false;
		}
	}

	if (finalise_session)
	{
		//Let the user review the session
		ncm_image_session &session = *(subject_.live_session());
		ncm_qcapture_end_session session_ui(session, study_name_, subject_.name_in_study(), this);
		session_ui.exec();

		//Finalise the session and update the database
		session.finalise();
		QSqlQuery query(database_);
		session.sql_update(query, "time_finished", session.time_finished());
		session.sql_update(query, "session_comments", session.notes());

		//Update the list of sequences to process offline (for now add them all - TO DO only add selected sequences)
		bool process_bad = true;
		session.sql_add_to_process_list(query, process_bad);

		//Tell the subject to move its current session to the sessions list
		subject_.add_live_to_previous();
	}

	//Clear up the motors scene
	if (motors_are_live_)
		motors_scene_.clear();

	//Reset the subject
	subject_.reset();

	//Reset the motors
	motors_are_live_ = false;
	apt_.destroy_controllers();	

	//Disconnect the camera
	disconnect_camera();

	//disconnect signals and end all the threads
	disconnect_session_signals_from_slots();

	//Empty the main queue
	ncm_qcapture_data_manager::Instance()->empty_main_queue();

	// Close the threads
	apt_timer_.stop();
	apt_thread_.quit();
	aligner_thread_.quit();
	saver_thread_.quit();
	processor_thread_.quit();

	//Disable session tab, go back to the home tab
	activity_state_ = activity_Home;
	update_session_controls();
}

void ncm_qcapture_gui::start_review_session()
{
	//Create new instance of search for subject GUI
	ncm_qcapture_search search_ui(subject_, study_name_, database_, ncm_qcapture_search::ReviewSession, this);
	search_ui.exec();

	if (subject_.ncm_id() == -1)
		//Not subject was selected, just return without doing anything
		return;


	//This all seems to take a little time, so show a busy cursor until the select subject
	//dialog box opens
	setCursor(QCursor(Qt::WaitCursor));

	//Initialise the review scene...
	// use Pan mode by default, with reverse zoom (wheel up = zoom in)
	ui.reviewView->usePanMode();
  ui.reviewView->reverse_zoom();
  ui.reviewView->setFocus();

  // tell graphics view to use this scene
  ui.reviewView->setScene(&review_scene_);
	review_scene_.set_image_processor(&imageProcessor_);
	review_scene_.set_annotation(&markup_);

	setCursor(QCursor(Qt::ArrowCursor));

	//Update the subject display and switch to the session tab
	activity_state_ = activity_Reviewing;
	update_review_display();
	update_session_controls();

	int i = ui.editWidget->indexOf(ui.markupTab);
	if (i >= 0)
		ui.editWidget->removeTab(i);
	i = ui.editWidget->indexOf(ui.editTab);
	if (i >= 0)
		ui.editWidget->removeTab(i);
	
	//Initialise subjects view model
	sessionModel_ = new QSqlQueryModel;
	sessionFilterModel_ = new QSortFilterProxyModel();
	sequenceModel_ = new QSqlQueryModel;
	sequenceFilterModel_ = new QSortFilterProxyModel();

	//populate the session list view with sessions from this subject
	QString session_query_str = "SELECT * FROM session WHERE subject_id = " + QString::number(subject_.ncm_id());

	//Query the database and get the column index of the subject name
	sessionModel_->setQuery(session_query_str, database_);


	int fieldNo = sessionModel_->record().indexOf("time_started");
	
	//Link the SQL model to the filter model so we can filter on subject name
	sessionFilterModel_->setSourceModel( sessionModel_ );
	sessionFilterModel_->setFilterKeyColumn(fieldNo);

	//Link the filter model to the list view to display the filtered subjects
	ui.sessionListView->setModel(sessionFilterModel_);
	ui.sessionListView->setModelColumn(fieldNo);
	ui.sessionListView->setSelectionBehavior(QAbstractItemView::SelectRows);
	ui.sessionListView->setSelectionMode(QAbstractItemView::SingleSelection);
  ui.sessionListView->show();

	//Connect up the session view selection model so we process selection changes
	connect(ui.sessionListView->selectionModel(), SIGNAL( selectionChanged (QItemSelection, QItemSelection)),
				this, SLOT(onSessionSelectionChanged()),
				Qt::UniqueConnection);

	QBrush b(QColor(255,0,0,255)); //Semi-opaque red lettering
	no_mosaic_msg_ = review_scene_.addSimpleText("No mosaic available to review", QFont("Arial", 48));
	no_mosaic_msg_->setBrush(b);
	no_mosaic_msg_->hide();

	//rotate the review scene if we need to
	if (preferences_.rotate_frames_)
	{
		QTransform t1 = ui.reviewView->transform();
		if (t1.m11()==1 && t1.m22()==1)
		{
			ui.reviewView->rotate(180);
		}
		//QPointF p = no_mosaic_msg_->pos();
		//no_mosaic_msg_->setPos(0,0);
		no_mosaic_msg_->rotate(180);
		//no_mosaic_msg_->setPos(-p);		
	}

	connect_review_signals_to_slots();

	//Quite frankly, if the user has no sessions we shouldn't be here!
	if (sessionFilterModel_->rowCount() > 0)
	{
		const QModelIndex i = sessionFilterModel_->index(sessionFilterModel_->rowCount()-1,0);
		ui.sessionListView->setCurrentIndex(i);
	}
}

void ncm_qcapture_gui::end_review_session()
{
	disconnect_review_signals_from_slots();

	//Reset the scene
	imageProcessor_.clear();
	review_scene_.clear();
	review_scene_.removeItem(no_mosaic_msg_);

	//Reset the subject
	subject_.reset();

	//Delete the pointers to the view models
	delete sessionModel_;
	delete sessionFilterModel_;
	delete sequenceModel_;
	delete sequenceFilterModel_;

	//Disable session tab, go back to the home tab
	activity_state_ = activity_Home;
	update_session_controls();
}

bool ncm_qcapture_gui::load_camera(bool use_existing)
{
	//Try and connect the device
	if ( !camera_.connect_device(use_existing) ) {
		vcl_cout << "Camera failed to load" << vcl_endl;
		return false;
	}

  update_framerate_controls();
  update_exposure_controls();
  update_gain_controls();

  // Use highest frame rate by default
  // This is done already by the Load Camera dialog
  // And is irritating if you select something other than the highest from the 
  // dialog, only for your selection to be reset.

  //ui.frameRateSelection->setCurrentIndex(0);

	//Start the device's live mode
	camera_.start_live();

	//Set exposure and gain to auto
	if (!camera_.is_auto_exposure())
		camera_.switch_auto_exposure(true);
	if (!camera_.is_auto_gain())
		camera_.switch_auto_gain(true);

	ui.stopLiveButton->setText( QString( "Stop live" ) );

	// draw scene and fit image in view
  redraw_scene();

	//Enable the camera and save sequence controls
	ui.tabDisplay->setEnabled( true );
	ui.tabCamera->setEnabled( true );

	ui.gainSlider->setEnabled( false ); //Start in auto mode
	ui.exposureSlider->setEnabled( false ); //Start in auto mode
	ui.gainSpinBox->setEnabled( false ); //Start in auto mode
	ui.exposureSpinBox->setEnabled( false ); //Start in auto mode
	ui.autoExposure->setChecked(true);
	ui.autoGain->setChecked(true);
	ui.actionLoad_camera->setEnabled( false ); //Can change camera but no point trying to load the same one
	ui.graphicsView->setVisible(true);
	return true;
}

void ncm_qcapture_gui::disconnect_camera()
{
	// Disconnect the camera
	if ( camera_.is_connected( ) )
		camera_.disconnect_device();

	ui.graphicsView->setVisible(false);
}

bool ncm_qcapture_gui::connect_to_database(const QString & username, const QString & password) 
{
	//Just tryopening a connecting to the nailfold database here... we'll play around later
	//with exactly what and where this should be done 
	//database_.database("QMYSQL", "nailfold_connection");
	database_ = QSqlDatabase::addDatabase("QMYSQL", "nailfold_connection");
  database_.setHostName("");
  database_.setDatabaseName("nailfold_initial_test");
  database_.setUserName(username);
  database_.setPassword(password);

	database_connected_ = database_.open();
	if (database_connected_)
		vcl_cout << "Connected successfully to the database" << vcl_endl;
	else
		vcl_cout << "Couldn't connect to database - what now?" << vcl_endl;

	return database_connected_;
}

void ncm_qcapture_gui::disconnect_database()
{
	database_.removeDatabase("nailfold_connection");
	database_.close();
}

bool ncm_qcapture_gui::load_user(const QString & username)
{
	//Also load in this user - are we just going to assume this will work given a successful connection
	//to the database??
	QSqlQuery query(database_);
	query.prepare("SELECT user_id FROM user WHERE username = '" + username + "'");
	if(!query.exec())
	{
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;
		return false; //Problem with database?
	}

	if (!query.next())
		return false; //Username not found?

	int user_id = query.value(0).toInt();

	//Try and load preferences for this user
	query.prepare("SELECT * FROM user_preferences WHERE user_id = " + QString::number(user_id));
	if(!query.exec())
		vcl_cout << query.lastError().databaseText().toStdString() << vcl_endl;

	if (!query.next())
		preferences_.sql_write(query, user_id);
	
	else
		preferences_.sql_read(query.record());

	//Also create this user object
	user_ = ncm_user(user_id);
	user_.set_name(username);
	//user_.set_host(host);

	//Success!
	return true;
}

void ncm_qcapture_gui::update_session_controls()
{
	switch (activity_state_)
	{
		case (activity_Login):
		{
			//Disable everything except the login page
			ui.stackedWidget->setCurrentWidget(ui.pageLogin);

			ui.pageLogin->setEnabled(true);
			ui.pageLogin->setVisible(true);

			ui.pageHome->setEnabled(false);
			ui.pageHome->setVisible(false);

			ui.pageSession->setEnabled(false);
			ui.pageSession->setVisible(false);

			ui.pageReview->setEnabled(false);
			ui.pageReview->setVisible(false);

			ui.menubar->setEnabled(false);
			ui.menubar->setVisible(false);

			ui.statusbar->setVisible(false);

			break;
		}
		case (activity_Home):
		{
			//Current page is home page
			ui.stackedWidget->setCurrentWidget(ui.pageHome);

			ui.pageLogin->setEnabled(false);
			ui.pageLogin->setVisible(false);

			ui.pageHome->setEnabled(true);
			ui.pageHome->setVisible(true);

			ui.pageSession->setEnabled(false);
			ui.pageSession->setVisible(false);

			ui.pageReview->setEnabled(false);
			ui.pageReview->setVisible(false);

			ui.menubar->setEnabled(true);
			ui.menubar->setVisible(true);
			ui.statusbar->setVisible(false);
			break;
		}

		case (activity_Session):
		{
			//Current page is session page, current tab is session tab, record but showing "start"
			ui.stackedWidget->setCurrentWidget(ui.pageSession);
			ui.pageLogin->setEnabled(false);
			ui.pageLogin->setVisible(false);

			ui.pageHome->setEnabled(false);
			ui.pageHome->setVisible(false);

			ui.pageSession->setEnabled(true);
			ui.pageSession->setVisible(true);

			ui.pageReview->setEnabled(false);
			ui.pageReview->setVisible(false);
			ui.statusbar->setVisible(true);

			ui.tabWidget->setCurrentWidget(ui.tabPosition);
			ui.startRecordingPushButton->setIcon(start_icon_);
			ui.startRecordingPushButton->setText("Start recording");
			ui.startRecordingPushButton->setToolTip("Start recording new sequence");

			//Start is enabled, stop is disabled
			ui.stopRecordingPushButton->setEnabled(false);
			ui.startRecordingPushButton->setEnabled(true);
			ui.endSessionPushButton->setEnabled(true);

			// Enable the camera controls
			ui.tabCamera->setEnabled( true );
			break;
		}

		case (activity_Recording):
		{
			// record button showing "pause", stop button enabled
			ui.startRecordingPushButton->setIcon(pause_icon_);
			ui.startRecordingPushButton->setToolTip("Pause recording");
			ui.startRecordingPushButton->setText("Pause recording");
			ui.stopRecordingPushButton->setEnabled(true);

			// Disable the camera controls
			ui.tabCamera->setEnabled( false );
			ui.endSessionPushButton->setEnabled( false );
			break;
		}

		case (activity_RecordingPaused):
		{
			// record button showing "start"
			ui.startRecordingPushButton->setIcon(start_icon_);
			ui.startRecordingPushButton->setToolTip("Resume recording");
			ui.startRecordingPushButton->setText("Resume recording");
			//As this mode can only be entered from activity_Recording
			//we don't need disable everything again
			break;
		}

		case (activity_Saving):
		{
			//Disable start/stop buttons (everything else that was disable should still be disabled)
			ui.startRecordingPushButton->setEnabled(false);
			ui.stopRecordingPushButton->setEnabled(false);

			//Display something in the status to say we're saving?
			break;
		}

		case (activity_Reviewing):
		{
			//Current page is reviewing page
			ui.stackedWidget->setCurrentWidget(ui.pageReview);
			ui.pageLogin->setEnabled(false);
			ui.pageLogin->setVisible(false);

			ui.pageHome->setEnabled(false);
			ui.pageHome->setVisible(false);

			ui.pageSession->setEnabled(false);
			ui.pageSession->setVisible(false);

			ui.pageReview->setEnabled(true);
			ui.pageReview->setVisible(true);

			break;
		}

		default:
			assert(false); //should never get here
	}
}

void ncm_qcapture_gui::update_tab_controls()
{
	if (!motors_are_live_)
	{
		//disable motor controls
		const int positionIndex = ui.tabWidget->indexOf(ui.tabPosition);
		if (positionIndex >= 0 )
			ui.tabWidget->removeTab(positionIndex);

		const int motorsIndex = ui.tabWidget->indexOf(ui.tabMotors);
		if (motorsIndex >= 0 )
			ui.tabWidget->removeTab(motorsIndex);

		const int autofocusIndex = ui.tabWidget->indexOf(ui.tabAutofocus);
		if (autofocusIndex >= 0 )
			ui.tabWidget->removeTab(autofocusIndex);

	}
	else if (!preferences_.use_autofocus_)
	{
		//disable autofocus
		const int autofocusIndex = ui.tabWidget->indexOf(ui.tabAutofocus);
		if (autofocusIndex >= 0 )
			ui.tabWidget->removeTab(autofocusIndex);
	}
	else
	{
		//Motors and autofocus enabled, add all the tabs (if they're not already there)
		const int positionIndex = ui.tabWidget->indexOf(ui.tabPosition);
		if (positionIndex < 0 )
			ui.tabWidget->insertTab(0, ui.tabPosition, "Frame position");

		const int motorsIndex = ui.tabWidget->indexOf(ui.tabMotors);
		if (motorsIndex < 0 )
			ui.tabWidget->insertTab(-1, ui.tabMotors, "Motors");

		
		const int autofocusIndex = ui.tabWidget->indexOf(ui.tabAutofocus);
		if (autofocusIndex < 0 )
			ui.tabWidget->insertTab(-1, ui.tabAutofocus, "Autofocus");

		//Make the motor position tab the current one
		ui.tabWidget->setCurrentWidget(ui.tabPosition);
	}
}

void ncm_qcapture_gui::update_framerate_controls()
{	
  //Get available frame rates and update frame selection combo box
	vcl_vector<double> available_FPS;
	double fps;

  // Disconnect signal to slot temporarily to avoid unnecessary event
  // handling.
  QObject::disconnect( ui.frameRateSelection, 
                       SIGNAL(currentIndexChanged(const QString)),
                       this, 
                       SLOT(on_frameRateSelection_currentIndexChanged(const QString)) );

	if ( camera_.get_available_FPS( available_FPS ) &&
       camera_.get_current_FPS( fps ) )
	{
		QString text_double;
		for (unsigned i = 0; i < available_FPS.size(); i++)
		{
			double fps_i = available_FPS[ i ];
			text_double.setNum( fps_i );
			ui.frameRateSelection->addItem( text_double );
		
			if ( fps_i == fps )
				ui.frameRateSelection->setCurrentIndex( i );
		}
	}

  // Reconnect signal to slot
  QObject::connect( ui.frameRateSelection, 
                    SIGNAL(currentIndexChanged(const QString)),
                    this, 
                    SLOT(on_frameRateSelection_currentIndexChanged(const QString)) );
}
 
void ncm_qcapture_gui::update_exposure_controls()
{
	//Get exposure range of camera and set spin box and slider
	double min_exposure, max_exposure;
	if ( camera_.get_exposure_range( min_exposure, max_exposure ) )
	{
		ui.exposureSpinBox->setRange( 0, 100 );
		ui.exposureSlider->setRange( 0, 100 );
		ui.exposureSpinBox->setExposureMin( min_exposure );
		ui.exposureSpinBox->setExposureMax( max_exposure );
	}
}
 
void ncm_qcapture_gui::update_gain_controls()
{
	//Get gain range of camera and set spin box and slider	
	double min_gain, max_gain;
	if ( camera_.get_gain_range( min_gain, max_gain ) )
	{
		ui.gainSpinBox->setRange( min_gain, max_gain );
		ui.gainSlider->setRange( int( min_gain ), int( max_gain ) );
	}
}
 
//: Disable the motor controls
void ncm_qcapture_gui::disable_motor_controls()
{
  const int motorsIndex = ui.tabWidget->indexOf(ui.tabMotors);
  //ui.tabWidget->setTabEnabled(motorsIndex, false);
	ui.tabWidget->removeTab(motorsIndex);

	const int positionIndex = ui.tabWidget->indexOf(ui.tabPosition);
  //ui.tabWidget->setTabEnabled(positionIndex, false);
	ui.tabWidget->removeTab(positionIndex);

	const int autofocusIndex = ui.tabWidget->indexOf(ui.tabAutofocus);
  //ui.tabWidget->setTabEnabled(autofocusIndex, false);
	ui.tabWidget->removeTab(autofocusIndex);

}
 
//: Update controller comboboxes that enable you to assign controllers to axes.
void ncm_qcapture_gui::update_apt_comboboxes()
{
  // Get list of IDs that are connected
  vcl_vector<long> id_vector = apt_.ids();

  // For each axis, remove the IDs of any controller currently assigned to 
  // another axis and set the current item index to that matching the currently
  // assigned controller (if any).

  ui.aptXCombobox->clear();
  ui.aptXCombobox->addItem("None");
  for (unsigned i = 0; i < id_vector.size(); ++i)
  {
    // Add controller to list if not already in use by another axis.
    if ((apt_.Y()->id() != id_vector[i]) && 
        (apt_.Z()->id() != id_vector[i]))
      ui.aptXCombobox->addItem(QString::number(id_vector[i]));

    // Store the index of the currently assigned controller.
    if (apt_.X()->id() == id_vector[i])
      ui.aptXCombobox->setCurrentIndex(ui.aptXCombobox->count()-1);
  }
    
  ui.aptYCombobox->clear();
  ui.aptYCombobox->addItem("None");
  for (unsigned i = 0; i < id_vector.size(); ++i)
  {
    // Add controller to list if not already in use by another axis.
    if ((apt_.X()->id() != id_vector[i]) && 
        (apt_.Z()->id() != id_vector[i]))
      ui.aptYCombobox->addItem(QString::number(id_vector[i]));

    // Store the index of the currently assigned controller.
    if (apt_.Y()->id() == id_vector[i])
      ui.aptYCombobox->setCurrentIndex(ui.aptYCombobox->count()-1);
  }

  ui.aptZCombobox->clear();
  ui.aptZCombobox->addItem("None");
  for (unsigned i = 0; i < id_vector.size(); ++i)
  {
    // Add controller to list if not already in use by another axis.
    if ((apt_.X()->id() != id_vector[i]) && 
        (apt_.Y()->id() != id_vector[i]))
      ui.aptZCombobox->addItem(QString::number(id_vector[i]));

    // Store the index of the currently assigned controller.
    if (apt_.Z()->id() == id_vector[i])
      ui.aptZCombobox->setCurrentIndex(ui.aptZCombobox->count()-1);
  }
}
 
//: Set the target position for X axis.
void ncm_qcapture_gui::set_apt_targetX(double distance_pixels)
{
  const double mm_per_pixel = 1.0 / preferences_.pixels_per_motor_mm_;
  double distance_mm = mm_per_pixel * distance_pixels;

  if (motors_are_live_)
  {
		apt_.X()->set_velocity(preferences_.dblclick_max_);
    apt_.X()->move_by(distance_mm);
  }
}

void ncm_qcapture_gui::set_apt_targetY(double distance_pixels)
{
  const double mm_per_pixel = 1.0 / preferences_.pixels_per_motor_mm_;
  double distance_mm = mm_per_pixel * distance_pixels;
  
  if (motors_are_live_)
  {
    apt_.Y()->set_velocity(preferences_.dblclick_max_);
    apt_.Y()->move_by(distance_mm);
  }
}

void ncm_qcapture_gui::set_apt_targetZ(double z)
{
  if (motors_are_live_)
  {
    apt_.Z()->set_velocity(1.0f);
    apt_.Z()->move_by(z);
  }
}

void ncm_qcapture_gui::home_x()
{
	if (!motors_are_live_)
    return;

  // Emit a signal telling apt_server to go home.
  const ncm_qapt_server::apt_axis axis = ncm_qapt_server::AptX;
  apt_.set_max_displacement(axis, 1e9);
  emit apt_move_home(axis, /* return_immediately = */ true);
}

void ncm_qcapture_gui::home_y()
{
	if (!motors_are_live_)
    return;

  // Emit a signal telling apt_server to go home.
  const ncm_qapt_server::apt_axis axis = ncm_qapt_server::AptY;
  apt_.set_max_displacement(axis, 1e9);
  emit apt_move_home(axis, /* return_immediately = */ true);
}

void ncm_qcapture_gui::home_z()
{
	if (!motors_are_live_)
    return;

  // Emit a signal telling apt_server to go home.
  const ncm_qapt_server::apt_axis axis = ncm_qapt_server::AptZ;
  apt_.set_max_displacement(axis, 1e9);
  emit apt_move_home(axis, /* return_immediately = */ true);
}


//
//: Autofocus procedure, implemented as a finite state machine.
void ncm_qcapture_gui::autoFocus(autoFocus_state_enum state)
{
  autoFocus_state_ = state;

  switch (autoFocus_state_)
  {
    case autoFocus_GetDirectionForward:
    {
      ui.autofocusTextEdit->setPlainText("");

      // Determine the direction of movement for maximum sharpness

      // Move camera until sharpness changes by some (absolute or relative) 
      // amount, or until a fixed distance (e.g. 1mm) has been traversed.
	    const double current_sharpness = processor_.image_sharpness();

      const double lower_limit = current_sharpness * 0.95;
      processor_.set_lower_sharpness_threshold(lower_limit);

      const double upper_limit = current_sharpness / 0.95;
      processor_.set_upper_sharpness_threshold(upper_limit);

      // Set limits on displacement from current position.
      apt_.set_reference_position(ncm_qapt_server::AptZ);

      const float max_displacement_z = 1.0f; // 1mm
      apt_.set_max_displacement(ncm_qapt_server::AptZ, 
                                max_displacement_z);

      vcl_stringstream ss;
      ss << "Searching forwards "
         << current_sharpness << " "
         << "[" << lower_limit << ", " << upper_limit << "]"
         << vcl_endl;
      ui.autofocusTextEdit->appendPlainText(ss.str().c_str());

      // Start searching.
      apt_.Z()->move_at_velocity(0.25f);
      break;
    }

    case autoFocus_GetDirectionReverse:
    {
      // Reverse the direction of the search if there was no improvement 
      // in the forward direction.

      // Redefine lower threshold to avoid triggering due to noise.
 	    const double current_sharpness = processor_.image_sharpness();

      // Decrease lower threshold since the motor has some 'inertia'.
      const double lower_limit = current_sharpness * 0.75;
      processor_.set_lower_sharpness_threshold(lower_limit);

      const double upper_limit = current_sharpness / 0.975;
      processor_.set_upper_sharpness_threshold(upper_limit);

      // Move camera until sharpness changes by some (absolute or relative) 
      // amount, or until a fixed distance has been traversed.
      // Set limits on displacement from current position.
      apt_.set_reference_position(ncm_qapt_server::AptZ);

      // Allow twice the travel in the opposite direction.
      const float max_displacement_z = 2.0f; // 2mm
      apt_.set_max_displacement(ncm_qapt_server::AptZ, 
                                max_displacement_z);

      vcl_stringstream ss;
      ss << "Searching backwards " 
               << current_sharpness << " "
               << "[" << processor_.lower_sharpness_threshold() << ", " 
                      << processor_.upper_sharpness_threshold() << "]"
               << vcl_endl;
      ui.autofocusTextEdit->appendPlainText(ss.str().c_str());

      // Start searching.
      apt_.Z()->move_at_velocity(-0.25f);

      break;
    }

    case autoFocus_FindPeak:
    {
      // Maintain the existing velocity, updating the optimal position every
      // time the upper threshold is exceeded.

      // Store the position at which sharpness has been maximized so far.
      autoFocus_peak_position_ = apt_.Z()->position();

      const double current_sharpness = processor_.image_sharpness();

      // Consider a significant fall as an indicator that the peak has
      // been passed.
      const double lower_limit = current_sharpness * 0.9;
      processor_.set_lower_sharpness_threshold(lower_limit);

      // Consider anything better than the current sharpness to be the peak.
      const double upper_limit = current_sharpness * 1.0;
      processor_.set_upper_sharpness_threshold(upper_limit);

      // Give the camera freedom to move as far as it needs to focus.
      apt_.set_max_displacement(ncm_qapt_server::AptZ, 1e9);

      // Everything else is handled by signals and slots.

      vcl_stringstream ss;
      ss << "FindPeak " 
               << autoFocus_peak_position_ << " "
               << current_sharpness << " "
               << "[" << lower_limit << ", " << upper_limit << "]"
               << vcl_endl;
      ui.autofocusTextEdit->appendPlainText(ss.str().c_str());

      break;
    }

    case autoFocus_GoToPeak:
    {
      // Disable limits as they are no longer needed.
      // (We know where we need to be.)
      processor_.set_lower_sharpness_threshold(-1e9);
      processor_.set_upper_sharpness_threshold(1e9);

      vcl_stringstream ss;
      ss << "Go to peak (" << autoFocus_peak_position_ << ")"
               << vcl_endl;
      ui.autofocusTextEdit->appendPlainText(ss.str().c_str());

      // Don't call apt_.move_to() here - this executes the move_to() function
      // in the GUI thread (i.e. the thread in which it was called).
      // To execute move_to() in the apt_thread_, you must emit a signal that
      // will be received by apt_'s move_to() slot, such that the function
      // is called within the apt thread.
      emit apt_move_to(ncm_qapt_server::AptZ, 
                       autoFocus_peak_position_, false);

      // Defer switching to autoFocus_Inactive until the moveComplete signal
      // is received.
      break;
    }

    case autoFocus_Interrupt:
    {
      // Something wants to interrupt the autoFocus process.
      // (e.g. user manually changes the Z velocity.)

      // Disable limits in sharpness and displacement, then switch to Inactive
      // mode.

      processor_.set_lower_sharpness_threshold(-1e9);
      processor_.set_upper_sharpness_threshold(1e9);

      apt_.set_max_displacement(ncm_qapt_server::AptZ, 1e9);

      autoFocus_state_ = autoFocus_Inactive;
      break;
    }

    case autoFocus_Error:
    {
      // Stop the motor.
      apt_.Z()->stop();

      // Reset the autofocus parameters.
      processor_.set_lower_sharpness_threshold(-1e9);
      processor_.set_upper_sharpness_threshold(1e9);

      apt_.set_max_displacement(ncm_qapt_server::AptZ, 1e9);

      const double reference_sharpness = processor_.image_sharpness();

      vcl_stringstream ss;
      ss << "Sharpness = " << reference_sharpness << vcl_endl;
      ui.autofocusTextEdit->appendPlainText(ss.str().c_str());

      autoFocus_state_ = autoFocus_Inactive;
      break;
    }

    default:
      break;
  }
}
 
//: Set up the scenes that display the camera stream, mosaic and motor positions
void ncm_qcapture_gui::initialize_scene()
{
	// tell main graphics view to use the main camera scene
	ui.graphicsView->setScene(&scene_);

	//also link up the mosaic scene and motor position scene
	ui.mosaicView->setScene(&mosaic_scene_);

	if (motors_are_live_)
	{
		ui.motorsView->setScene(&motors_scene_);
		motors_scene_.setSceneRect(-12.5, -12.5, 25, 25);
		//ui.motorsView->setSceneRect(-12.5, -12.5, 25, 25);

		//Draw motor limits and centre-lines in motor position scene
		motors_scene_.addRect(-12.5, -12.5, 25, 25, QPen(QColor(255,0,0)), QBrush());
		motors_scene_.addLine(-12.5, 0, 12.5, 0, QPen(QColor(0,255,255)));
		motors_scene_.addLine(0, -12.5, 0, 12.5, QPen(QColor(0,255,255)));

		//Create the currrent frame position and add it to the scene (don't worry about
		//its position for now, that will be sorted when redraw scene is called)
		current_frame_position_ = motors_scene_.addRect(-1.0, -0.75, 2, 1.5,
			QPen(QColor(0,255,0)), QBrush());
		mosaic_frame_positions_.clear();

		QTransform t = ui.motorsView->viewportTransform();

		double scale = t.m22();
		double w = ui.motorsView->width() - 10;
		double h = ui.motorsView->height() - 10;

		ui.motorsView->scale( w/(25*scale), h/(25*scale));

	}

	if (preferences_.rotate_frames_)
	{
		QTransform t1 = ui.graphicsView->transform();
		if (t1.m11()==1 && t1.m22()==1)
		{
			ui.graphicsView->rotate(180);
			ui.mosaicView->rotate(180);
			ui.motorsView->rotate(180);
			scene_.rotate_text();
		}
		
	}
	//ui.graphicsView->reverse_zoom();

  // use Pan mode by default, with reverse zoom (wheel up = zoom in)
  //ui.graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);
  //ui.graphicsView->reverse_zoom();
}

//: Initialize the saver class with defaults
void ncm_qcapture_gui::initialize_motors()
{
	apt_.create_controllers();

	motors_are_live_ = (apt_.n_controllers() > 0);

	if (!motors_are_live_)
	{
		QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.setText(
			"Could not connect to motors \n"
			"Are they connected and switched on? \n");
    msgBox.setInformativeText(
			"Click 'Ok' to continue (motor controls will be disabled) \n\n"
			"To reconnect after switching the motors on, \n"
			"please close then restart the software.");
    msgBox.setStandardButtons(QMessageBox::Ok);

    msgBox.exec();
		disable_motor_controls();
		return;
	}

	//If we've got this far motors are live so...
	// Read previous values from the sql database
	get_apt_sql_values();  
	update_apt_comboboxes();
	apt_.X()->set_home_position(12.5f);
  apt_.Y()->set_home_position(12.5f);
  apt_.Z()->set_home_position(12.5f);
	
	initialize_apt_timer();

  //Home the motors on startup - we probably don't *need* to do this every time.
	//But for the sake of 45secs it does seem to make all future driving of the motors
	//more reliable
	if (preferences_.home_motors_on_startup_)
	{	
		//Home the motors on startup
		if (!apt_.X()->homed())
			home_x();
		if (!apt_.Y()->homed())
			home_y();
		if (!apt_.Z()->homed())
			home_z();
	}
	
}

//: Initialize the saver class with defaults
void ncm_qcapture_gui::initialize_saver()
{
  // Get the root dir from the Registry (if it exists).
  saver_.setRootDir( preferences_.root_dir_ );
	
}
 
//: Initialize motor control timer
void ncm_qcapture_gui::initialize_apt_timer()
{
	if (!apt_timer_.isActive())
	{
		const int clicksPerSecond = 10;
		apt_timer_.start(/* interval = */ 1000.0 / clicksPerSecond);
	}
}
 
//: Initialize the processor thread
void ncm_qcapture_gui::initialize_processor_thread()
{
  processor_thread_.setObjectName("processor_thread");
	processor_.moveToThread(&processor_thread_);
	processor_thread_.start();
}
 
//: Initialize the saver thread
void ncm_qcapture_gui::initialize_saver_thread()
{
  saver_thread_.setObjectName("saver_thread");
	saver_.moveToThread(&saver_thread_);
  saver_thread_.start();
}

 
//: Get motor properties from Windows registry.
void ncm_qcapture_gui::get_apt_sql_values()
{
  // Get motor assignments.
  int id = 0;
  int index = -1;
	bool reversed;
	QSqlQuery query("SELECT * from motor_ids WHERE motor_id = 1", database_);
	
	if (query.next())
	{
		id = query.record().value("AptXId").toInt();
		apt_.set_X_id(id);

		id = query.record().value("AptYId").toInt();
		apt_.set_Y_id(id);

		id = query.record().value("AptZId").toInt();
		apt_.set_Z_id(id);

		//Set reversed state for each motor
		reversed = query.record().value("XReversed").toBool();
		apt_.X()->set_reverse(reversed);
		ui.aptXReverse->setChecked(reversed);
		reversed = query.record().value("YReversed").toBool();
		apt_.Y()->set_reverse(reversed);
		ui.aptYReverse->setChecked(reversed);
		reversed = query.record().value("ZReversed").toBool();
		apt_.Z()->set_reverse(reversed);
		ui.aptZReverse->setChecked(reversed);

		/*
    index = ui.aptXCombobox->findText(QString::number(id));
    if (index != -1)
    {
      ui.aptXCombobox->setCurrentIndex(index);

      // Update GUI manually because setCurrentItem does not trigger
      // the activated() signal programmatically.
      on_aptXCombobox_activated(index);
    }

		
		index = ui.aptYCombobox->findText(QString::number(id));
    if (index != -1)
    {
      ui.aptYCombobox->setCurrentIndex(index);

      // Update GUI manually because setCurrentItem does not trigger
      // the activated() signal programmatically.
      on_aptYCombobox_activated(index);
    }

		
		index = ui.aptZCombobox->findText(QString::number(id));
    if (index != -1)
    {
      ui.aptZCombobox->setCurrentIndex(index);
     
      // Update GUI manually because setCurrentItem does not trigger
      // the activated() signal programmatically.
      on_aptZCombobox_activated(index);
    }
		*/
  }
	update_apt_comboboxes();
}
 
//: Initialize the APT (motor) thread
void ncm_qcapture_gui::initialize_apt_thread()
{

  apt_thread_.setObjectName("apt_thread");

  apt_.moveToThread(&apt_thread_);

  // Register the enum so that the correct metafunctions can be called.
  qRegisterMetaType<ncm_qapt_server::apt_axis>("ncm_qapt_server::apt_axis");

  // Start the thread, permitting discrete moves (e.g. homing) but not 
  // continuous ones (e.g. dragging the joystick).
  apt_thread_.start();
}

//Update the display field showing notes for the current subject
void ncm_qcapture_gui::update_subject_display()
{
	ui.subjectNameLabel->setText("<b>Subject: <\b> " + subject_.name_in_study());
	ui.studyNameLabel->setText("<b>Study: <\b> " + study_name_);
	QString subject_text = 
		"Subject ID: " + QString::number(subject_.ncm_id()) + "\n"
		"Dominant hand : " + QString(ncm_hand::toString(subject_.dominant_hand()).c_str()) + "\n"
		"Subject notes: " + subject_.notes();
	ui.subjectDetailsTextEdit->setText(subject_text);

	QString sessions_text;
	QDateTime date = subject_.last_imaged();
	if (date.isNull())
		sessions_text += "(no previous sessions)";
	else
		sessions_text += subject_.previous_sessions_string();

	ui.previousSessionTextEdit->setText(sessions_text);	
}

//Update the display field showing notes for the current subject
void ncm_qcapture_gui::update_session_display()
{
	QString sequences_text;

	if (subject_.live_session()->image_sequences()->empty())
		sequences_text = "(none)";
	else
		sequences_text = subject_.live_session()->sequences_string();
	ui.currentSessionTextEdit->setText(sequences_text);
}

//Update the display field showing notes for the current subject
void ncm_qcapture_gui::update_review_display()
{
	ui.reviewSubjectNameLabel->setText("<b>Subject: <\b> " + subject_.name_in_study());
	ui.reviewStudyNameLabel->setText("<b>Study: <\b> " + study_name_);
	QString subject_text = 
		"Subject ID: " + QString::number(subject_.ncm_id()) + "\n"
		"Dominant hand : " + QString(ncm_hand::toString(subject_.dominant_hand()).c_str()) + "\n"
		"Subject notes: " + subject_.notes();
	ui.reviewSubjectTextEdit->setText(subject_text);
}

//
void ncm_qcapture_gui::update_review_controls()
{

}

void ncm_qcapture_gui::sync_sequence_name(ncm_image_sequence &sequence)
{
	//If sequence name is empty or already synced then do nothing
	if (
		sequence.sequence_id() < 0 ||
		sequence.sequence_name().isEmpty() || 
		sequence.sequence_name() == sequence.images_sub_dir())
		return;

	//Create a database query and try and sync the names
	QSqlQuery query(database_);
	bool success = sequence.sync_name_and_dir(query);

	//If it doesn't work warn the user
	if (!success)	
	{
		QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.setText("Unable to update sequence folder \n"
			"Images remain in " + sequence.images_root_dir() + "/" + sequence.images_sub_dir());
		msgBox.exec();
		return;
	}
	update_session_display();
}

void ncm_qcapture_gui::interpolate_motor_positions(ncm_image_sequence &sequence)
{
	motor_positions_.clear();
	motor_positions_ = vcl_vector<vcl_vector <double>>(0, vcl_vector<double>(3,0));
	motor_times_.clear();
}

 
//: Connect signals to slots (not surprisingly...)
void ncm_qcapture_gui::connect_session_signals_to_slots()
{
  // Signals that trigger slots in the main thread

	//--------------------------------------------------------------------------
	//Processing frames (from initialize_processor_thread())

  // Tag the frame with motor positions before processing.
  QObject::connect( &camera_, SIGNAL(frame_ready()),
                    this, SLOT(tag_frame()) );

  // Camera tells processor there's a frame ready
  QObject::connect( this, SIGNAL(frame_tagged()),
                    &processor_, SLOT(process_frame()) );

	//Processor tells GUI there's a frame ready to draw
  QObject::connect( &processor_, SIGNAL(frame_to_draw()),
                    this, SLOT(redraw_scene()) );

  //Processor tells GUI the saver is done
  QObject::connect( &processor_, SIGNAL(processing_finished(int)),
						        this, SLOT(onProcessingFinished(int)) );

	//--------------------------------------------------------------------------
	//Saving frames (from initialize_save_thread())
	QObject::connect( &processor_, SIGNAL(frame_to_save(ncm_video_frame_header)),
                    this, SLOT(log_frame_properties(ncm_video_frame_header)) );

  // Processor tells saver to save a frame
  QObject::connect( this, SIGNAL(frame_logged(int)),
                    &saver_, SLOT(save_frame(int)) );
  QObject::connect( &saver_, SIGNAL(frame_saved(int)),
                    this, SLOT(onFrameSaved(int)) );
  QObject::connect( &saver_, SIGNAL(saving_finished()),
                    this, SLOT(onSavingFinished()) );

	QObject::connect( this, SIGNAL(finalize_sequence()),
                    this, SLOT(onFinalizeSequence()) );

	//--------------------------------------------------------------------------
	//Aligning/mosaicing frames (from initialize_align_thread())
	QObject::connect( &processor_, SIGNAL(frame_to_align(int)),
                    this, SLOT(onFrameToAlign(int)) ); 

  QObject::connect( this, SIGNAL(frame_to_align()),
                    &aligner_, SLOT(alignNextFrame()) ); 

  QObject::connect( &aligner_, SIGNAL(framesAligned(int, int, double, double, double)),
                    this, SLOT(onFramesAligned(int, int, double, double, double)) );
  QObject::connect( &aligner_, SIGNAL(alignmentFinished()),
                    this, SLOT(onAlignmentFinished()) );

	//--------------------------------------------------------------------------
	//Moving the motors from the GUI (from initialize_apt_thread())
  QObject::connect( this, SIGNAL(apt_move_zero(ncm_qapt_server::apt_axis, bool)),
                    &apt_, SLOT(move_zero(ncm_qapt_server::apt_axis, bool)) );

  QObject::connect( this, SIGNAL(apt_move_home(ncm_qapt_server::apt_axis, bool)),
                    &apt_, SLOT(move_home(ncm_qapt_server::apt_axis, bool)) );

  QObject::connect( this, SIGNAL(apt_move_home(ncm_qapt_server::apt_axis, float, bool)),
                    &apt_, SLOT(move_home(ncm_qapt_server::apt_axis, float, bool)) );

  QObject::connect( this, SIGNAL(apt_move_to(ncm_qapt_server::apt_axis, float, bool)),
                    &apt_, SLOT(move_to(ncm_qapt_server::apt_axis, float, bool)) );

  QObject::connect( this, SIGNAL(apt_move_by(ncm_qapt_server::apt_axis, float, bool)),
                    &apt_, SLOT(move_by(ncm_qapt_server::apt_axis, float, bool)) );

	QObject::connect( &apt_, SIGNAL(move_complete(ncm_qapt_server::apt_axis)),
                    this, SLOT(onAptMoveComplete(ncm_qapt_server::apt_axis)) );

	//Moving the motors from the scene
  QObject::connect( &scene_, SIGNAL(joystickMoved(double, double)),
                    this, SLOT(onSceneJoystickMoved(double, double)) );
  QObject::connect( &scene_, SIGNAL(joystickReleased()),
                    this, SLOT(onSceneJoystickReleased()) );

  QObject::connect( &scene_, SIGNAL(doubleClicked(double, double, Qt::MouseButton)),
										this, SLOT(onSceneDoubleClicked(double, double, Qt::MouseButton)) );
  QObject::connect( &scene_, SIGNAL(wheelMoved( int )),
                    this, SLOT(on_graphicsView_zoomed( int )) );

	//Moving motors during autofocus
  QObject::connect( &processor_, SIGNAL(lower_sharpness_exceeded()),
                    this, SLOT(onLowerSharpnessThresholdCrossed()) );
  QObject::connect( &processor_, SIGNAL(upper_sharpness_exceeded()),
                    this, SLOT(onUpperSharpnessThresholdCrossed()) );

	//--------------------------------------------------------------------------
	//Getting motor information at set time interval
  QObject::connect( &apt_timer_, SIGNAL(timeout()),
                    this, SLOT(onAptTimer()) );
	
}

//: Connect signals to slots (not surprisingly...)
void ncm_qcapture_gui::disconnect_session_signals_from_slots()
{
  // Signals that trigger slots in the main thread

	//--------------------------------------------------------------------------
	//Processing frames (from initialize_processor_thread())

  // Tag the frame with motor positions before processing.
  QObject::disconnect( &camera_, SIGNAL(frame_ready()),
                    this, SLOT(tag_frame()) );

  // Camera tells processor there's a frame ready
  QObject::disconnect( this, SIGNAL(frame_tagged()),
                    &processor_, SLOT(process_frame()) );

	//Processor tells GUI there's a frame ready to draw
  QObject::disconnect( &processor_, SIGNAL(frame_to_draw()),
                    this, SLOT(redraw_scene()) );

  //Processor tells GUI the saver is done
  QObject::disconnect( &processor_, SIGNAL(processing_finished(int)),
						        this, SLOT(onProcessingFinished(int)) );

	//--------------------------------------------------------------------------
	//Saving frames (from initialize_save_thread())
	QObject::disconnect( &processor_, SIGNAL(frame_to_save(ncm_video_frame_header)),
                    this, SLOT(log_frame_properties(ncm_video_frame_header)) );

  // Processor tells saver to save a frame
  QObject::disconnect( this, SIGNAL(frame_logged(int)),
                    &saver_, SLOT(save_frame(int)) );
  QObject::disconnect( &saver_, SIGNAL(frame_saved(int)),
                    this, SLOT(onFrameSaved(int)) );
  QObject::disconnect( &saver_, SIGNAL(saving_finished()),
                    this, SLOT(onSavingFinished()) );

	QObject::disconnect( this, SIGNAL(finalize_sequence()),
                    this, SLOT(onFinalizeSequence()) );

	//--------------------------------------------------------------------------
	//Aligning/mosaicing frames (from initialize_align_thread())
	QObject::disconnect( &processor_, SIGNAL(frame_to_align(int)),
                    this, SLOT(onFrameToAlign(int)) ); 

  QObject::disconnect( this, SIGNAL(frame_to_align()),
                    &aligner_, SLOT(alignNextFrame()) ); 

  QObject::disconnect( &aligner_, SIGNAL(framesAligned(int, int, double, double, double)),
                    this, SLOT(onFramesAligned(int, int, double, double, double)) );
  QObject::disconnect( &aligner_, SIGNAL(alignmentFinished()),
                    this, SLOT(onAlignmentFinished()) );

	//--------------------------------------------------------------------------
	//Moving the motors from the GUI (from initialize_apt_thread())
  QObject::disconnect( this, SIGNAL(apt_move_zero(ncm_qapt_server::apt_axis, bool)),
                    &apt_, SLOT(move_zero(ncm_qapt_server::apt_axis, bool)) );

  QObject::disconnect( this, SIGNAL(apt_move_home(ncm_qapt_server::apt_axis, bool)),
                    &apt_, SLOT(move_home(ncm_qapt_server::apt_axis, bool)) );

  QObject::disconnect( this, SIGNAL(apt_move_home(ncm_qapt_server::apt_axis, float, bool)),
                    &apt_, SLOT(move_home(ncm_qapt_server::apt_axis, float, bool)) );

  QObject::disconnect( this, SIGNAL(apt_move_to(ncm_qapt_server::apt_axis, float, bool)),
                    &apt_, SLOT(move_to(ncm_qapt_server::apt_axis, float, bool)) );

  QObject::disconnect( this, SIGNAL(apt_move_by(ncm_qapt_server::apt_axis, float, bool)),
                    &apt_, SLOT(move_by(ncm_qapt_server::apt_axis, float, bool)) );

	QObject::disconnect( &apt_, SIGNAL(move_complete(ncm_qapt_server::apt_axis)),
                    this, SLOT(onAptMoveComplete(ncm_qapt_server::apt_axis)) );

	//Moving the motors from the scene
  QObject::disconnect( &scene_, SIGNAL(joystickMoved(double, double)),
                    this, SLOT(onSceneJoystickMoved(double, double)) );
  QObject::disconnect( &scene_, SIGNAL(joystickReleased()),
                    this, SLOT(onSceneJoystickReleased()) );

  QObject::disconnect( &scene_, SIGNAL(doubleClicked(double, double, Qt::MouseButton)),
										this, SLOT(onSceneDoubleClicked(double, double, Qt::MouseButton)) );
  QObject::disconnect( &scene_, SIGNAL(wheelMoved( int )),
                    this, SLOT(on_graphicsView_zoomed( int )) );

	//Moving motors during autofocus
  QObject::disconnect( &processor_, SIGNAL(lower_sharpness_exceeded()),
                    this, SLOT(onLowerSharpnessThresholdCrossed()) );
  QObject::disconnect( &processor_, SIGNAL(upper_sharpness_exceeded()),
                    this, SLOT(onUpperSharpnessThresholdCrossed()) );

	//--------------------------------------------------------------------------
	//Getting motor information at set time interval
  QObject::disconnect( &apt_timer_, SIGNAL(timeout()),
                    this, SLOT(onAptTimer()) );
	
}

//: Connect signals to slots (not surprisingly...)
void ncm_qcapture_gui::connect_review_signals_to_slots()
{
  QObject::connect( ui.autoContrast, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setAutoContrast(bool)) );
  QObject::connect( ui.autoContrast, SIGNAL(toggled(bool)),
                    this, SLOT(updateContrastControls()) );

	// Zoom controls
  QObject::connect( ui.zoomFitButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(fit_both()));
  QObject::connect( ui.zoomFitWButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(fit_width()));
  QObject::connect( ui.zoomFitHButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(fit_height()));

	QObject::connect( ui.reviewView, SIGNAL(changed()),
                    this, SLOT(updateZoomControls()) );

  // make connections between UI widgets and local member variables
  // we could possibly make these slots part of the reviewView rather than the
  // scene so that the connections can be set up in QDesigner (but it's not a
  // big deal at the moment).
  /*QObject::connect( ui.show_edges_checkBox, SIGNAL(toggled(bool)),
                    &scene_, SLOT(set_edges_visible(bool)) );
  QObject::connect( ui.show_lines_checkBox, SIGNAL(toggled(bool)),
                    &scene_, SLOT(set_lines_visible(bool)) );
  QObject::connect( ui.show_lines_checkBox, SIGNAL(toggled(bool)),
                    ui.snap_checkBox, SLOT(setEnabled(bool)) );
  QObject::connect( ui.snap_checkBox, SIGNAL(toggled(bool)),
                    &scene_, SLOT(set_snapping(bool)) );
  QObject::connect( ui.show_vessels_checkBox, SIGNAL(toggled(bool)),
                    &scene_, SLOT(set_vessels_visible(bool)) );

  QObject::connect( ui.show_lines_checkBox, SIGNAL(toggled(bool)),
                    &imageProcessor_, SLOT(find_ridges()) );
  QObject::connect( ui.show_edges_checkBox, SIGNAL(toggled(bool)),
                    &imageProcessor_, SLOT(find_edges()) );

  QObject::connect( ui.actionToggleImageControls, SIGNAL(toggled(bool)),
                    ui.imageGroup, SLOT(setVisible(bool)) );
  QObject::connect( ui.actionToggleZoomControls, SIGNAL(toggled(bool)),
                    ui.zoomGroup, SLOT(setVisible(bool)) );
  QObject::connect( ui.actionToggleGrid, SIGNAL(toggled(bool)),
                    &scene_, SLOT(set_grid_visible(bool)) );

  QObject::connect( &imageProcessor_, SIGNAL(imageLoaded(bool)),
                    this, SLOT(onImageProcessor_imageLoaded(bool)) );*/

	/*
  // When vessel combobox changes selection, select that vessel and zoom in
  QObject::connect( &scene_, SIGNAL(selectionChanged(int)),
                    ui.vesselCombo, SLOT(setCurrentIndex(int)) );
  QObject::connect( &scene_, SIGNAL(selectionChanged(int)),
                    this, SLOT(updateVesselControls()) );
  QObject::connect( &scene_, SIGNAL(selectionChanged(int)),
                    this, SLOT(applyVesselZoom(int)) );

  // When a vessel is picked from the list, make it selected in the scene, zoom
  // in on it, and update the vessel controls
  QObject::connect( ui.vesselCombo, SIGNAL(currentIndexChanged(int)),
                    &scene_, SLOT(select_vessel(int)) );
  // Treat any click, regardless of whether it changes the selection, as a
  // reason to redraw and select
  QObject::connect( ui.vesselCombo, SIGNAL(activated(int)),
                    &scene_, SLOT(select_vessel(int)) );
  // This just works nicely - when sweeping over the list, the screen updates
  // in realtime.
  QObject::connect( ui.vesselCombo, SIGNAL(highlighted(int)),
                    &scene_, SLOT(select_vessel(int)) );*/

	/*
  // React to changes in scene
  QObject::connect( &scene_, SIGNAL(vesselAdded()),
                    this, SLOT(updateVesselList()) );
  QObject::connect( &scene_, SIGNAL(vesselDeleted()),
                    this, SLOT(updateVesselList()) );
  QObject::connect( &scene_, SIGNAL(vesselChanged()),
                    this, SLOT(updateVesselList()) );

  QObject::connect( ui.reviewView, SIGNAL(previous()),
                    this, SLOT(on_prevVesselButton_clicked()) );
  QObject::connect( ui.reviewView, SIGNAL(next()),
                    this, SLOT(on_nextVesselButton_clicked()) );

  // Changing the image grade
  QObject::connect( gradeRadioGroup_, SIGNAL(buttonClicked(int)),
                    this, SLOT(setGrade(int)) );

  // Set vessel properties
  QObject::connect( rowRadioGroup_, SIGNAL(buttonClicked(int)),
                    this, SLOT(updateRowButtons(int)) );
  QObject::connect( ui.distalRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setDistal(bool)) );
  QObject::connect( ui.nondistalRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setNondistal(bool)) );

  QObject::connect( sizeRadioGroup_, SIGNAL(buttonClicked(int)),
                    this, SLOT(updateSizeButtons(int)) );
  QObject::connect( ui.sizeNormalRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setSizeNormal(bool)) );
  QObject::connect( ui.sizeEnlargedRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setSizeEnlarged(bool)) );
  QObject::connect( ui.sizeGiantRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setSizeGiant(bool)) );
  QObject::connect( ui.sizeIrregularRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setSizeIrregular(bool)) );
  
  QObject::connect( shapeRadioGroup_, SIGNAL(buttonClicked(int)),
                    this, SLOT(updateShapeButtons(int)) );
  QObject::connect( ui.shapeNormalRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setShapeNormal(bool)) );
  QObject::connect( ui.shapeTortuousRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setShapeTortuousOnly(bool)) );
  QObject::connect( ui.shapeRamifiedRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setShapeRamifiedOnly(bool)) );

  QObject::connect( &imageProcessor_, SIGNAL(rawImageChanged()),
                    ui.reviewView, SLOT(repaint()) );
  QObject::connect( &imageProcessor_, SIGNAL(rawImageChanged()),
                    this, SLOT(updateContrastControls()) );
  QObject::connect( &imageProcessor_, SIGNAL(imageTypeChanged()),
                    this, SLOT(updateToolbar()) );
  QObject::connect( &imageProcessor_, SIGNAL(imageTypeChanged()),
                    this, SLOT(updateStageControls()) );

  // Signals that should update the status bar
  QObject::connect( &imageProcessor_, SIGNAL(statusChanged(const QString&, double)),
                    this, SLOT(updateStatus(const QString&, double)) );
  QObject::connect( &scene_, SIGNAL(apexLengthChanged(double)),
                    this, SLOT(onApexLengthChanged(double)) );*/

  //connectSetFocus();
}

void ncm_qcapture_gui::connectSetFocus()
{
  // After any of the following events, return the focus to the
  // graphics view to process keyboard events
  QObject::connect( ui.brightnessRSlider, SIGNAL(sliderReleased()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.contrastRSlider, SIGNAL(sliderReleased()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.autoContrast, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setFocus()));

  QObject::connect( ui.zoomSlider, SIGNAL(sliderReleased()),
                    ui.reviewView, SLOT(setFocus()));
  //QObject::connect( ui.zoomEdit, SIGNAL(editFinished()),
  //                  ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.zoomInButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.zoomOutButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.zoomFitButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.zoomFitWButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.zoomFitHButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));
	/*
  QObject::connect( ui.vesselCombo, SIGNAL(currentIndexChanged(int)),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.firstVesselButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.prevVesselButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.nextVesselButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.lastVesselButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.deleteAllVessels, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.deleteAllHaemorrhages, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));

  QObject::connect( ui.stackedWidget, SIGNAL(currentChanged(int)),
                    ui.reviewView, SLOT(setFocus()));

  QObject::connect( ui.sizeNormalRadio, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.sizeEnlargedRadio, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.sizeGiantRadio, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setFocus()));
  
  QObject::connect( ui.shapeNormalRadio, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.shapeTortuousRadio, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.shapeRamifiedRadio, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setFocus()));

  QObject::connect( ui.show_edges_checkBox, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setFocus()));

  QObject::connect( ui.show_lines_checkBox, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.compSizeSlider, SIGNAL(sliderReleased()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.snap_checkBox, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.show_vessels_checkBox, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setFocus()));

  QObject::connect( ui.prevStageButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));
  QObject::connect( ui.nextStageButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(setFocus()));*/
										
}

//: Connect signals to slots (not surprisingly...)
void ncm_qcapture_gui::disconnect_review_signals_from_slots()
{
	QObject::disconnect( ui.autoContrast, SIGNAL(toggled(bool)),
                    ui.reviewView, SLOT(setAutoContrast(bool)) );
  QObject::disconnect( ui.autoContrast, SIGNAL(toggled(bool)),
                    this, SLOT(updateContrastControls()) );

	// Zoom controls
  QObject::disconnect( ui.zoomFitButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(fit_both()));
  QObject::disconnect( ui.zoomFitWButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(fit_width()));
  QObject::disconnect( ui.zoomFitHButton, SIGNAL(clicked()),
                    ui.reviewView, SLOT(fit_height()));
}

//: Set the graphicsView zoom factor from a percentage value
void ncm_qcapture_gui::setZoomFrom(int value)
{
	const double minimum_scale = ui.reviewView->fitScale();
  double new_scale = static_cast<double>(value-1) / 100;

  if (new_scale < minimum_scale)
    new_scale = minimum_scale;

  ui.reviewView->setScale(new_scale);
}

void ncm_qcapture_gui::align_full_mosaic(ncm_image_sequence & sequence)
{
	sequence.display_details();

	vcl_string images_datafile = sequence.images_full_dir().toStdString() + "/" + sequence.images_datafile().toStdString();

	//Get filenames
	vcl_vector<vcl_string> filenames;
	sequence.t_read(images_datafile);

	vcl_string image_dir = vul_file::dirname(images_datafile);

	int num_frames = sequence.num_frames();
	//Load motor data - do we really need this?
	vcl_vector<long> timestamps;
	vcl_vector<QPointF> motor_positions;
	vcl_vector<double> motor_zoom;
	vcl_vector<double> sharpness_vec;

	// Conversion from mm (the units of the motor positions) and pixels in the
	// scene. This may not be totally accurate.
	const double pixels_per_mm = preferences_.pixels_per_motor_mm_;

	// These are needed to discard the first few lines if they are not
	// valid numbers.
	bool is_first_valid_position = true;
	QPointF datum(0,0);

	for (int i = 0; i < num_frames; ++i)
	{
		ncm_video_frame_header &header = sequence.frame_header(i);
				
		vcl_string frame_name = image_dir + "/" + header.frame_name_;
		//vcl_cout << frame_name << vcl_endl;
		filenames.push_back(frame_name);

		QPointF pos(header.motor_position_[0],
								header.motor_position_[1]);
			
		bool const position_is_valid = true; //TO DO - what are valid values?

		if (!position_is_valid)
			motor_positions.push_back(QPointF(0.0,0.0));
		else
		{
			// Transform to reflect the reversed motors and differences in scale.
			// (TODO: Apply the negation for reversed motors in the capture 
			//        software, as we can't assume it here and it isn't recorded.)
			pos = -pos * pixels_per_mm;
			pos.setY(-pos.y());

			// Set the datum position to the first valid position.
			if (is_first_valid_position)
			{
				datum = pos;
				is_first_valid_position = false;
			}

			// Store subsequent position relative to datum.
			// Note that the first position will therefore be (0,0).
			motor_positions.push_back(pos - datum);
		}

		// Read the motor positions (zoom).
		motor_zoom.push_back(header.motor_position_[2]);

		// Read the estimated sharpness (not used at the moment).
		sharpness_vec.push_back(header.sharpness_);

		//Read the timestamp
		timestamps.push_back(header.frame_time_.toString("hhmmsszzz").toInt());
	}

	//Initialize the aligner
	offline_aligner_.setFilenames(filenames);

	//Use the motor positions?
	offline_aligner_.resetDisplacements();
	for (unsigned i = 1; i < num_frames; ++i)
	{
		// Compute displacements from motor positions.
		int di = motor_positions[i].x() - motor_positions[i-1].x();
		int dj = motor_positions[i].y() - motor_positions[i-1].y();

		offline_aligner_.setDisplacement(i, QPoint(di, dj));
	}

	int radius = 100;
	offline_aligner_.set_filter(ncm_frame_aligner::filter_g1d);
	offline_aligner_.set_n_levels(0);
	offline_aligner_.set_di_radius(radius);
  offline_aligner_.set_dj_radius(radius);

	offline_aligner_.reset();

	for (unsigned i = 1; i < num_frames; ++i)
	{
		offline_aligner_.alignNextFrame();
		vcl_cout << "Aligned frame " << i << " of " << num_frames << vcl_endl;
	}

	/*
	//offline_aligner_.moveToThread(&offline_thread_);
	//offline_thread_.start();

	// This -> aligner
  QObject::connect( &offline_aligner_, SIGNAL(framesAligned()),
                    &offline_aligner_, SLOT(alignNextFrame()) ); 

  QObject::connect( &offline_aligner_, 
                      SIGNAL(alignmentFinished()),
                    this, 
                      SLOT(make_full_mosaic()) );

	//Now just need to the aligner to align the first frame (note this will run
	//in the main thread, but the subsequent signals and slots should run in the offline frame
	//to align the remaining frames)!
	offline_aligner_.alignNextFrame();*/

	make_full_mosaic(sequence);

}

void ncm_qcapture_gui::make_full_mosaic(ncm_image_sequence &sequence)
{
	vcl_cout << "Making offline mosaic!" << vcl_endl; 
	
	//Instantiate the stitcher, connect up sticther signals and move to the aligner thread
	ncm_qcapture_mosaic_maker stitcher(offline_aligner_.ni(), offline_aligner_.nj());

	//Get aligner displacements and start mosaicing...
	vcl_vector<QPoint> displacements = offline_aligner_.displacements();
	vcl_vector<vcl_string> filenames = offline_aligner_.filenames();
	stitcher.setFilenames(filenames);
	stitcher.makeMosaicFrom(displacements);

	vcl_cout << "Mosaic finished!" << vcl_endl;
	vcl_cout << "Mosaic size, w = " << stitcher.mosaic().ni() << ", h = " << stitcher.mosaic().nj() << vcl_endl;

	//Save the mosaic and map counts
	sequence.set_full_mosaic("full_mosaic.png"); //TO DO make this a user preference...
	vcl_string mosaic_name = sequence.images_full_dir().toStdString() + "/" + sequence.full_mosaic().toStdString();
	vil_save(stitcher.mosaic(), mosaic_name.c_str() );

	sequence.set_full_count_map("full_count_map.png"); //TO DO make this a user preference...
	vcl_string count_name = sequence.images_full_dir().toStdString() + "/" + sequence.full_count_map().toStdString();
	vil_image_view<vxl_uint_16> count_map;
	vil_convert_cast(stitcher.count_image(), count_map);
	vil_save(count_map, count_name.c_str() );

	//Write the sequence to the dabatase
	QSqlQuery query(database_);
	sequence.sql_update(query, "full_mosaic", sequence.full_mosaic());
	sequence.sql_update(query, "full_count_map", sequence.full_count_map());

	/*
	// this -> stitcher_
  QObject::connect( this,  SIGNAL(make_mosaic(vcl_vector<QPoint>)),
                    stitcher_, SLOT(makeMosaicFrom(vcl_vector<QPoint>)) ); 

  // stitcher_ -> this
  QObject::connect( stitcher_, SIGNAL(mosaicUpdated(int)),
                    this, SLOT(onMosaicUpdated(int)) );
  QObject::connect( stitcher_, SIGNAL(mosaicFinished(bool)),
                    this, SLOT(onMosaicFinished(bool)) );

	//Move to thread
	stitcher_->moveToThread(&aligner_thread_);

	//Get aligner displacements and start mosaicing...
	vcl_vector<QPoint> displacements = aligner_.displacements();
	emit make_mosaic(displacements);*/
}
