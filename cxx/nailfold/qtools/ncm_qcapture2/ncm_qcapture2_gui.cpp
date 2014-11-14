#include "ncm_qcapture2_gui.h"

#include <vcl_iomanip.h>
#include <vil/vil_math.h>
#include <vil/vil_load.h>
#include <vil/vil_convert.h>

#include <vsl/vsl_quick_file.h>

#include <qcore/qcore_convert_image.h>

#include <nailfold/ncm_apt_controller.h>

//
//: Constructor
ncm_qcapture2_gui::ncm_qcapture2_gui(QWidget *parent, Qt::WFlags flags)
: QMainWindow(parent, flags), 
	raw_pixmap_item_d_(0),
  brightness_(2),
  contrast_(2),
	current_camera_(0),
	save_from_camera_(2, false),
	diff_frame_(2, false),
	flipped_h_(2, false),
	flipped_v_(2, false),
  // Constants
  motors_are_live_(false),
  pixels_per_mm_(1158.2),
  n_sharpness_(10)
{
	// setup the UI
  ui.setupUi(this);

	//Set up the radio button group to select the current camera
	currentCameraRadioGroup = new QButtonGroup(this);
  currentCameraRadioGroup->addButton(ui.activeCameraButton1, 0);
  currentCameraRadioGroup->addButton(ui.activeCameraButton2, 1);
  currentCameraRadioGroup->setExclusive(true);
  //ui.activeCameraButton1->setChecked(true);
  QObject::connect( currentCameraRadioGroup, SIGNAL(buttonClicked( int )),
                    this, SLOT(change_active_camera( int )) );

	//Initialise the vector to store the two cameras
	cameras_.resize(2);
	cameras_[0] = new ncm_dmk_camera();
	cameras_[1] = new ncm_dmk_camera();

	//Initialise the vector of savers
	savers_.resize(2);
	savers_[0] = new ncm_qcapture2_saver();
	savers_[1] = new ncm_qcapture2_saver();

	//get brightness and contrast from sliders
	brightness_.fill(ui.brightnessSlider->value());
	contrast_.fill(ui.contrastSlider->value());

	//Set up raw pixmaps for the two camera views
	raw_pixmap_items_.resize(2);
	raw_pixmap_items_.fill( 0 );

	//Set up color tables for the two cameras
	raw_colour_tables_.resize(2);
	update_raw_colour_table(brightness_[0], contrast_[0], raw_colour_tables_[0]);
	update_raw_colour_table(brightness_[1], contrast_[1], raw_colour_tables_[1]);
	update_diff_colour_table(ui.brightnessSlider->value(), ui.contrastSlider->value());

	//Initialise the two camera scenes
	scenes_.resize(2);
	scenes_[0] = new QGraphicsScene();
	scenes_[1] = new QGraphicsScene();

	// tell graphics view to use this scene
	ui.cameraView1->setScene(scenes_[ 0 ]);
	ui.cameraView2->setScene(scenes_[ 1 ]);

	ui.lastFrameView->setScene(&scene_d_);

	//Link the data managers queue to the cameras
	cameras_[ 0 ]->set_queue( ncm_qcapture2_data_manager::Instance()->main_queue(0) );
	cameras_[ 1 ]->set_queue( ncm_qcapture2_data_manager::Instance()->main_queue(1) );

  initialize_registry_accessor();
  initialize_saver();

	update_apt_comboboxes();
  initialize_apt_spinboxes();

  /*if (apt_.n_controllers() == 0)
  {
    // Disable the motor tab if no controllers exist.
    int motorsIndex = ui.tabWidget->indexOf(ui.tabMotors);
    ui.tabWidget->setTabEnabled(motorsIndex, false);
  }
  else
  {
    initialize_apt_thread();
  }*/ //MOTOR_SWITCH_OFF

  initialize_processor_thread();
  initialize_saver_thread();

  initializeAlignerThread();

  //vcl_string artefacts_filename = 
  //    "U:/projects/nailfold/capture/2013_02_21/Left/Digit4/x300/15_01_31/_diff.png";
  //vil_image_view<vxl_byte> artefact_image_byte = 
  //    vil_load(artefacts_filename.c_str());
  //vil_convert_cast(artefact_image_byte, artefact_image_);
  //vil_math_scale_and_offset_values(artefact_image_, 1.0, -128);

  initializeStatusbar();

  connect_signals_to_slots();
}

//
//: Destructor
ncm_qcapture2_gui::~ncm_qcapture2_gui()
{
}

 
// Events
void ncm_qcapture2_gui::closeEvent(QCloseEvent *ev)
{
  if (savers_[0]->is_saving() || savers_[1]->is_saving())
  {
    QMessageBox msgBox;
    msgBox.setText("Still saving frames - please wait.");
    msgBox.setIcon(QMessageBox::Information);
    msgBox.exec();
    ev->ignore();
    return;
  }

	// Disconnect the cameras
	if ( cameras_[0]->is_connected( ) )
		cameras_[0]->disconnect_device();

	if ( cameras_[1]->is_connected( ) )
		cameras_[1]->disconnect_device();

	// Close the threads
	saver_thread_.quit();
	processor_thread_.quit();
  //apt_thread_.quit();

  ev->accept();
}

//:Compute new min/max contrasts given new brightness/contrast values
void ncm_qcapture2_gui::update_contrasts(int brightness, int contrast, int & c_min, int & c_max) {
	
	int contrast_min = (255-brightness) - (255-contrast)/2;
  int contrast_max = (255-brightness) + (255-contrast)/2;

	
	if (contrast_min < 0)
    c_min = 0;
  else if (contrast_min > 255)
    c_min = 255;
  else
    c_min = contrast_min;

  if (contrast_max < 0)
    c_max = 0;
  else if (contrast_max > 255)
    c_max = 255;
  else
    c_max = contrast_max;

}

//
//: Change colormap for raw images to effect a change in contrast
void ncm_qcapture2_gui::update_raw_colour_table(int brightness, int contrast, QVector<QRgb> &raw_colour_table)
{
	int contrast_min, contrast_max;
	update_contrasts(brightness, contrast, contrast_min, contrast_max);

  raw_colour_table.resize(256);

  // fill everything from first to constrast_min_ with 0 (black)
  for (int i = 0; i < contrast_min; i++)
    raw_colour_table[i] = qRgb(0,0,0);

  // fill everything from constrast_max_ to end with 255 (white)
  for (int i = contrast_max; i < 256; i++)
    raw_colour_table[i] = qRgb(255,255,255);

  // fill values in between with a linear gradient
  int denom = contrast_max - contrast_min;
  for (int i = contrast_min; i < contrast_max; i++)
  {
    int grey = (255 * (i - contrast_min)) / denom;
    raw_colour_table[i] = qRgb(grey,grey,grey);
  }
}

void ncm_qcapture2_gui::update_diff_colour_table(int brightness, int contrast)
{
	int contrast_min, contrast_max;
	update_contrasts(brightness, contrast, contrast_min, contrast_max);

  raw_colour_table_d_.resize(256);

  // fill everything from first to constrast_min_ with dark blue
  for (int i = 0; i < contrast_min; i++)
    raw_colour_table_d_[i] = qRgb(0,0,128);

  // fill everything from constrast_max_ to end with dark red
  for (int i = contrast_max; i < 256; i++)
    raw_colour_table_d_[i] = qRgb(128,0,0);

  // fill values in between with a linear gradient
  int period = (contrast_max + 1 - contrast_min)>>2;
	int d1 = contrast_min + (period>>1);
	int d2 = d1 + period;
	int d3 = d2 + period;
	int d4 = d3 + period;
	int step = 256 / period;

	int r = 0;
	int g = 0;
	int b = 128;
  for (int i = contrast_min; i < contrast_max; i++)
  {
		if (i < d1) {
			b+= step; b = min(b, 255);
			
		}
		else if (i < d2){
			b = 255;
			g+= step; g = min(g, 255);
			
		}
		else if (i < d3 ) {
			b-= step; b = max(b, 0);			
			g = 255;
			r+= step; r = min(r, 255);
		}
		else if (i < d4 ) {
			b = 0;
			g-= step; g = max(g, 0);
			r = 255;
		}
		else {
			g = 0;
			r-= step;
		}
    raw_colour_table_d_[i] = qRgb(r,g,b);
  }
}

//
//  Private slots
//

//:Action when either load camera button is clicked
void ncm_qcapture2_gui::on_loadCameraButton1_clicked()
{
	if ( load_camera(0) ) {
		ui.useCamera1Button->setEnabled( true );
		ui.useCamera1Button->setChecked( true );
		ui.activeCameraButton1->setEnabled( true );
		ui.activeCameraButton1->setChecked( true );
	}
}

void ncm_qcapture2_gui::on_loadCameraButton2_clicked()
{
	if ( load_camera(1) ) {
		ui.useCamera2Button->setEnabled( true );
		ui.useCamera2Button->setChecked( true );
		ui.activeCameraButton2->setEnabled( true );
		ui.activeCameraButton2->setChecked( true );
	}
}

//:Change which camera is currently active
void ncm_qcapture2_gui::change_active_camera( int cam_num )
{
	//cam_num is camera id from radio button clicked
	current_camera_ = cam_num;

	vcl_cout << "Camera " << cam_num + 1 << " selected" << vcl_endl;

	//Update controls to reflect state of current camera
	update_framerate_controls();
  update_exposure_controls();
  update_gain_controls();

	if ( cameras_[current_camera_]->is_auto_exposure() ) {
		ui.autoExposure->setChecked( true );
		ui.exposureSlider->setEnabled( false );
		ui.exposureSpinBox->setEnabled( false );
	}
	else {
		//switch off camera's auto exposure and get current exposure value
		double exposure, value;
		cameras_[ current_camera_ ]->get_current_exposure( exposure );
		value = ui.exposureSpinBox->valueFromExposure( exposure );

		//enable exposure slider and spin box and set to current auto value
		ui.exposureSlider->setEnabled( true );
		ui.exposureSpinBox->setEnabled( true );
		ui.exposureSpinBox->setValue( value );
		ui.exposureSlider->setValue( int( value ) );
		ui.autoExposure->setChecked( false );
	}

	if ( cameras_[current_camera_]->is_auto_gain() ) {
		ui.autoGain->setChecked( true );
		ui.gainSlider->setEnabled( false );
		ui.gainSpinBox->setEnabled( false );
	}
	else {
		//switch off camera's auto gain and get current gain value
		double gain;
		cameras_[ current_camera_ ]->get_current_gain( gain );

		//enable gain slider and spin box and set to current auto value
		ui.gainSlider->setEnabled( true );
		ui.gainSpinBox->setEnabled( true );
		ui.gainSpinBox->setValue( gain );
		ui.gainSlider->setValue( int( gain ) );
		ui.autoGain->setChecked( false );
	}

	ui.brightnessSlider->setValue( brightness_[ current_camera_ ] );
	ui.contrastSlider->setValue( contrast_[ current_camera_ ] );

	if ( cameras_[current_camera_]->is_live() ) {
		ui.stopLiveButton->setText( QString( "Stop live" ) );
	}
	else {
		ui.stopLiveButton->setText( QString( "Start live" ) );
	}
}

//----------------------------------------------------------------------------
// Display control callbacks
void ncm_qcapture2_gui::on_brightnessSlider_sliderMoved(int value)
{
  brightness_[current_camera_] = value;
	
  // update colormap and redraw
	update_raw_colour_table(brightness_[ current_camera_ ], contrast_[ current_camera_ ], raw_colour_tables_[ current_camera_ ]);
  redraw_scene();	
}

void ncm_qcapture2_gui::on_contrastSlider_sliderMoved(int value)
{
  contrast_[current_camera_] = value;
	
  // update colormap and redraw
  update_raw_colour_table(brightness_[ current_camera_ ], contrast_[ current_camera_ ], raw_colour_tables_[ current_camera_ ]);
  redraw_scene();
}

void ncm_qcapture2_gui::on_stopLiveButton_clicked( )
{
	if ( !ui.stopLiveButton->text().compare( QString( "Start live" ) ) ) {
		cameras_[current_camera_]->start_live();
		ui.stopLiveButton->setText( QString( "Stop live" ) );
	}
	else {
		cameras_[current_camera_]->stop_live();
		ui.stopLiveButton->setText( QString( "Start live" ) );
	}
}

void ncm_qcapture2_gui::on_fliphImage1Button_toggled() {
	ui.cameraView1->scale(-1, 1);
	flipped_h_[0] = !flipped_h_[0];
}

void ncm_qcapture2_gui::on_fliphImage2Button_toggled() {
	ui.cameraView2->scale(-1, 1);
	flipped_h_[1] = !flipped_h_[1];
}

void ncm_qcapture2_gui::on_flipvImage1Button_toggled() {
	ui.cameraView1->scale(1, -1);
	flipped_v_[0] = !flipped_v_[0];
}

void ncm_qcapture2_gui::on_flipvImage2Button_toggled() {
	ui.cameraView2->scale(1, -1);
	flipped_v_[1] = !flipped_v_[1];
}


//----------------------------------------------------------------------------
// Camera control callbacks
void ncm_qcapture2_gui::on_frameRateSelection_currentIndexChanged( const QString &text ) 
{
	double frames_per_sec = text.toDouble();
	vcl_cout << "New FPS = " << frames_per_sec << vcl_endl;
	
	//Stop live mode, change the FPS, then restart live mode
	cameras_[ current_camera_ ]->set_FPS( frames_per_sec );
}

void ncm_qcapture2_gui::on_gainSpinBox_valueChanged( double value )
{
	ui.gainSlider->setValue( int( value ) );
	cameras_[ current_camera_ ]->set_gain( value );

}
void ncm_qcapture2_gui::on_gainSlider_sliderMoved( int value )
{
	ui.gainSpinBox->setValue( double( value ) );
	cameras_[ current_camera_ ]->set_gain( double( value ) );
}
void ncm_qcapture2_gui::on_autoGain_stateChanged( int state )
{
	if ( state )
	{
		//disable slider and spin box
		ui.gainSlider->setEnabled( false );
		ui.gainSpinBox->setEnabled( false );

		//switch on camera's auto gain
		cameras_[ current_camera_ ]->switch_auto_gain( true );

		vcl_cout << "Auto gain switched on" << vcl_endl;
	}
	else
	{
		//switch off camera's auto gain and get current gain value
		double gain;
		cameras_[ current_camera_ ]->switch_auto_gain( false );
		cameras_[ current_camera_ ]->get_current_gain( gain );

		//enable gain slider and spin box and set to current auto value
		ui.gainSlider->setEnabled( true );
		ui.gainSpinBox->setEnabled( true );
		ui.gainSpinBox->setValue( gain );
		ui.gainSlider->setValue( int( gain ) );
		vcl_cout << "Auto gain switched off" << vcl_endl;
	}
}
void ncm_qcapture2_gui::on_exposureSpinBox_valueChanged( double value )
{
	double exposure = ui.exposureSpinBox->exposureFromValue( value );
	ui.exposureSlider->setValue( int( value ) );
	cameras_[ current_camera_ ]->set_exposure( exposure );
}
void ncm_qcapture2_gui::on_exposureSlider_sliderMoved( int value )
{
	double exposure = ui.exposureSpinBox->exposureFromValue( value );
	ui.exposureSpinBox->setValue(double( value ));
	cameras_[ current_camera_ ]->set_exposure( exposure );
}
void ncm_qcapture2_gui::on_autoExposure_stateChanged( int state )
{
	if ( state )
	{
		//disable slider and spin box
		ui.exposureSlider->setEnabled( false );
		ui.exposureSpinBox->setEnabled( false );

		//switch on camera's auto exposure
		cameras_[ current_camera_ ]->switch_auto_exposure( true );

		vcl_cout << "Auto exposure switched on" << vcl_endl;
	}
	else
	{
		//switch off camera's auto exposure and get current exposure value
		double exposure, value;
		cameras_[ current_camera_ ]->switch_auto_exposure( false );
		cameras_[ current_camera_ ]->get_current_exposure( exposure );
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

//----------------------------------------------------------------------------
// Save sequence callbacks
void ncm_qcapture2_gui::on_sequenceLengthSpinBox_valueChanged(double value)
{
}
void ncm_qcapture2_gui::on_saveDirSelect_clicked()
{
	// get filename from dialog box
	QString fileName = QFileDialog::getExistingDirectory(this, 
                       tr("Select folder to save images in"), "", 
                       QFileDialog::ShowDirsOnly);

  // Return if cancelled.
  if (fileName.isEmpty())
    return;

	// Update the save_dir text edit box
	ui.saveDirTextEdit->setText( fileName );

  // Update the save dir in the saver_ object and Registry value
	savers_[0]->setRootDir(fileName.toStdString());
	savers_[1]->setRootDir(fileName.toStdString());
  ra_.write_string("RootDir", fileName.toStdString());
}

void ncm_qcapture2_gui::on_patientNameTextEdit_textChanged()
{
	QString patient_name = ui.patientNameTextEdit->toPlainText();
	savers_[0]->setPatientName( patient_name.toStdString() );
	savers_[1]->setPatientName( patient_name.toStdString() );
  ra_.write_string("PatientName", patient_name.toStdString());
}
void ncm_qcapture2_gui::on_sequenceNameTextEdit_textChanged()
{
	QString sequence_name = ui.sequenceNameTextEdit->toPlainText();
	savers_[0]->setSequenceName( sequence_name.toStdString() );
	savers_[1]->setSequenceName( sequence_name.toStdString() );
}

void ncm_qcapture2_gui::on_handComboBox_currentIndexChanged ( const QString &text ) 
{
	savers_[0]->setHand( text.toStdString() );
	savers_[1]->setHand( text.toStdString() );
}
void ncm_qcapture2_gui::on_digitComboBox_currentIndexChanged ( const QString &text ) 
{
	savers_[0]->setDigit( text.toStdString() );
	savers_[1]->setDigit( text.toStdString() );
}
void ncm_qcapture2_gui::on_magComboBox_currentIndexChanged ( const QString &text ) 
{
	savers_[0]->setMag( text.toStdString() );
	savers_[1]->setMag( text.toStdString() );
}
void ncm_qcapture2_gui::on_camera1SuffixTextEdit_textChanged ( ) 
{
	savers_[0]->setCameraSuffix( ui.camera1SuffixTextEdit->toPlainText().toStdString() );
}
void ncm_qcapture2_gui::on_camera2SuffixTextEdit_textChanged ( ) 
{
	savers_[1]->setCameraSuffix( ui.camera2SuffixTextEdit->toPlainText().toStdString() );
}
void ncm_qcapture2_gui::on_useCamera1Button_clicked()
{
	save_from_camera_[ 0 ] = ui.useCamera1Button->isChecked();
}

void ncm_qcapture2_gui::on_useCamera2Button_clicked()
{
	save_from_camera_[ 1 ] = ui.useCamera2Button->isChecked();
}

void ncm_qcapture2_gui::on_saveButton_clicked()
{
	//If no camera connected do nothing
	if ( !cameras_[0]->is_connected() && !cameras_[1]->is_connected() ) {
		vcl_cout << "No cameras connected" << vcl_endl;
		return;
	}

	//If neither camera selected, do nothing
	if ( !save_from_camera_[0] && !save_from_camera_[1] ) {
		vcl_cout << "Neither camera selected for save" << vcl_endl;
		return;
	}

	//
	// We definitely have some saving to do...
	//

  // Disable the camera and save sequence controls
	ui.tabCamera->setEnabled( false );
	ui.tabCapture->setEnabled( false );

	//Tell both savers to update their frame names
	savers_[0]->update_frame_name();
	savers_[1]->update_frame_name();

	// Make the current sequence directory if it doesn't already exist 
	//(we can do this with either saver, regardless of which cameras are connected)
	vcl_string filename = savers_[0]->makeSaveDir();

	// Define a subdirectory for the current capture
	vcl_string sub_dir = QTime::currentTime().toString("hh_mm_ss").toStdString();

	// Write out a file with capture properties
	filename.append( "sequence_properties_" );
	filename.append( sub_dir );
	filename.append( ".txt" );

	property_fs_.open(filename.c_str());
  write_properties( property_fs_ );
  add_properties_header();

	for (int cam_num = 0; cam_num < 2; cam_num++) {

		if ( save_from_camera_[ cam_num ] ) {
			
			// Set the number of frames to save
			double seq_len = ui.sequenceLengthSpinBox->value();
			double fps = cameras_[ cam_num ]->get_current_FPS();
			processor_.set_num_frames_to_save( int( fps * seq_len ), cam_num );

			//Set the subdir for this saver
			savers_[ cam_num ]->setSubDir(sub_dir);

			// Turn saving flag on
			processor_.start_processing(cam_num);
		}
	}
}
 
//: Add a line to the properties file summarizing the current frame.
void ncm_qcapture2_gui::add_properties_header()
{
  // Write the frame number, time stamp, motor position (x, y, z), sharpness,
  // and so on.
  property_fs_ << vcl_endl
               << vcl_setw(12) << "No."
               << vcl_setw(12) << "Time"
               << vcl_setw(12) << "MotorX"
               << vcl_setw(12) << "MotorY"
               << vcl_setw(12) << "MotorZ"
               << vcl_setw(12) << "Sharpness"
               << vcl_endl;
}
 
//: Add metadata to the most recently captured frame.
void ncm_qcapture2_gui::tag_frame( int cam_num, bool fliph, bool flipv )
{
	if ( diff_frame_[ cam_num ] ) {
		emit frame_tagged( cam_num, true, fliph, flipv );
		diff_frame_[ cam_num ] = false;
	}
	else
		emit frame_tagged( cam_num, false, fliph, flipv );
}
void ncm_qcapture2_gui::tag_frame1() {
	tag_frame(0, flipped_h_[0], flipped_v_[0] );
}
void ncm_qcapture2_gui::tag_frame2() {
	tag_frame(1, flipped_h_[1], flipped_v_[1] );
}

//: Add a line to the properties file summarizing the current frame.
void ncm_qcapture2_gui::log_frame_properties(int frame_number)
{
  if (!property_fs_.is_open())
    return;

  QSharedPointer<ncm_video_frame>& frame = 
    ncm_qcapture2_data_manager::Instance()->save_queue( 0 )->tail();

  float x = 0, y = 0, z = 0;
  frame->get_motor_xyz(x, y, z);

  double sharpness = 0.0;
  sharpness = frame->sharpness();

  vcl_string time_str = "";
  time_str = frame->frame_time().toString("hhmmsszzz").toStdString();

  // Write the frame number, time stamp, motor position (x, y, z), sharpness,
  // and so on.
  property_fs_ << vcl_setw(12) << frame_number
               << vcl_setw(12) << time_str
               << vcl_setw(12) << x
               << vcl_setw(12) << y
               << vcl_setw(12) << z
               << vcl_setw(12) << sharpness
               << vcl_endl;

  double progress = static_cast<double>(frame_number) /
                    processor_.num_frames_to_save( 0 );//bob
  updateProgress(progress);

  emit frame_logged(frame_number);
}
 
//: Update the progress bar when a frame is saved.
void ncm_qcapture2_gui::on_frame_saved(int frame_number, int cam_num)
{
  // Don't update the progress bar if the processor is still processing.
  if (processor_.is_processing(0) || processor_.is_processing(1))
    return;

	
	//Compute progress for this camera
  double progress = static_cast<double>(frame_number) /
                    processor_.num_frames_to_save( cam_num );

	//This may cause the progress bar to jump back and forward
	//but it'll only be stationary when both cameras finished
	updateProgress(progress);
}

//: Respond when all frames have been saved.
void ncm_qcapture2_gui::on_saving_finished(int cam_num)
{
  vcl_cout << "Camera " << cam_num+1 << " finished saving" << vcl_endl;
}
 
//: Redraw items on canvas
void ncm_qcapture2_gui::redraw_scene()
{
	redraw_scene( current_camera_ ); 
}

void ncm_qcapture2_gui::redraw_scene( int cam_num )
{
	QImage* qimage = ncm_qcapture2_data_manager::Instance()->get_qimage( cam_num );

  if (qimage->isNull())
    return;

	/*vil_image_view<vxl_byte> vxl_image;
  qcore_convert_image(vxl_image, *qimage);
  

  if (ui.correctImageCheckbox->isChecked() &&
      artefact_image_.ni() != 0)
  {
    vil_image_view<int> vxl_image_int;
    vil_convert_cast(vil_plane(vxl_image, 0), vxl_image_int);

    vil_math_image_difference(vxl_image_int, artefact_image_, vxl_image_int);

    vil_math_truncate_range(vxl_image_int, 0, 255);
    vil_convert_cast(vxl_image_int, vil_plane(vxl_image, 0));

    qcore_convert_image(*qimage, vil_plane(vxl_image, 0));
  }*/

  draw_pixmap(qimage, cam_num);
  static int di_cumul = 0, dj_cumul = 0;

  if (ui.stabilizeCheckbox->isChecked())
  {
    /*if (!aligner_.is_ready())
      aligner_.set_destination(vxl_image);
    aligner_.set_source(vxl_image);

    int di = 0, dj = 0;
    aligner_.align_src_to_dest(di, dj);
    aligner_.swap_src_with_dest();

    di_cumul += di;
    dj_cumul += dj;
    raw_pixmap_items_[ cam_num ]->setPos(di_cumul, dj_cumul);*/
  }
  else
  {
    di_cumul = 0;
    dj_cumul = 0;
    raw_pixmap_items_[ cam_num ]->setPos(0, 0);
  }

  draw_crosshairs( cam_num );

	scenes_[ cam_num ]->update();

	/*
  // Compute, display and store the sharpness of the image.
  QString msg;
  if (apt_.n_controllers() > 0)
  {
    msg += "X: " + QString::number(apt_.X()->absolute_position()) + " ";
    msg += "Y: " + QString::number(apt_.Y()->absolute_position()) + " ";
    msg += "Zabs: " + QString::number(apt_.Z()->absolute_position()) + " ";
    msg += "Zrel: " + QString::number(apt_.Z()->position()) + " ";
  }
  msg += "Sh: " + QString::number(processor_.image_sharpness());

  updateStatus(msg);
	*/
}
 
//:
void ncm_qcapture2_gui::draw_diff_frame()
{
	QImage *qimage = ncm_qcapture2_data_manager::Instance()->get_qdiffframe();

	qimage->setColorTable(raw_colour_table_d_);

  if (raw_pixmap_item_d_ == 0)
		raw_pixmap_item_d_ = scene_d_.addPixmap( QPixmap::fromImage( *qimage ));
  else
		raw_pixmap_item_d_->setPixmap(QPixmap::fromImage( *qimage ));

	scene_d_.setSceneRect(QRectF(0, 0, qimage->width(), qimage->height() ));
	scene_d_.update();
	ui.lastFrameView->fitInView(scene_d_.sceneRect(), Qt::KeepAspectRatio);
}

void ncm_qcapture2_gui::on_processing_finished()
{
  property_fs_.close();

	// Enable the camera and save sequence controls
	ui.tabCamera->setEnabled( true );
	ui.tabCapture->setEnabled( true );
}
 
//: Update the scene at regular intervals determined by a timer.
void ncm_qcapture2_gui::onAptTimer()
{
  // It's ok to call the move_at_velocity() function directly here because it
  // always returns immediately, unlike move_by().

  /*if (motors_are_live_)
  {
    if (ui.aptXVelocity->value() != 0)
    {
      float velocity = static_cast<float>(ui.aptXVelocity->value()) / 100.0;
      apt_.set_max_displacement(ncm_qapt_server::AptX, 1e9);
      apt_.X()->move_at_velocity(velocity);
    }
    if (ui.aptYVelocity->value() != 0)
    {
      float velocity = static_cast<float>(ui.aptYVelocity->value()) / 100.0;
      apt_.set_max_displacement(ncm_qapt_server::AptY, 1e9);
      apt_.Y()->move_at_velocity(velocity);
    }
    if (ui.aptZVelocity->value() != 0)
    {
      float velocity = static_cast<float>(ui.aptZVelocity->value()) / 100.0;
      apt_.set_max_displacement(ncm_qapt_server::AptZ, 1e9);
      apt_.Z()->move_at_velocity(velocity);
    }
  }

  ui.aptXPosition->setValue(apt_.X()->position());
  ui.aptYPosition->setValue(apt_.Y()->position());
  ui.aptZPosition->setValue(apt_.Z()->position());
	*/ //MOTOR_SWITCH_OFF
}

 
//: Change the controller for each axis
void ncm_qcapture2_gui::on_aptXCombobox_activated(int index)
{/*
  // If combobox is empty (i.e. the list has just been cleared) then do nothing.
  switch (index)
  {
    case -1:
      ui.aptXHome->setEnabled(false);
      return;
    case 0:
      ui.aptXHome->setEnabled(false);
      apt_.set_X_id(-1);
      break;
    default:
      ui.aptXHome->setEnabled(true);
      apt_.set_X_id(ui.aptXCombobox->currentText().toLong());
      ra_.write_numeric("AptXId", apt_.X()->id());
  }

  update_apt_comboboxes();
	*/ //MOTOR_SWITCH_OFF
}
void ncm_qcapture2_gui::on_aptYCombobox_activated(int index)
{
	/*
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
      ra_.write_numeric("AptYId", apt_.Y()->id());
  }

  update_apt_comboboxes();
	*/ //MOTOR_SWITCH_OFF
}
void ncm_qcapture2_gui::on_aptZCombobox_activated(int index)
{
	/*
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
      ra_.write_numeric("AptZId", apt_.Z()->id());
  }

  update_apt_comboboxes();
	*/ //MOTOR_SWITCH_OFF
}

//
//: Reverse directions of relative movements.
void ncm_qcapture2_gui::on_aptXReverse_stateChanged(int state)
{
	/*
  apt_.X()->set_reverse(state == Qt::Checked);
  ra_.write_boolean("XReversed", state == Qt::Checked);
	*/ //MOTOR_SWITCH_OFF
}
void ncm_qcapture2_gui::on_aptYReverse_stateChanged(int state)
{
	/*
  apt_.Y()->set_reverse(state == Qt::Checked);
  ra_.write_boolean("YReversed", state == Qt::Checked);
	*/ //MOTOR_SWITCH_OFF
}
void ncm_qcapture2_gui::on_aptZReverse_stateChanged(int state)
{
	/*
  apt_.Z()->set_reverse(state == Qt::Checked);
  ra_.write_boolean("ZReversed", state == Qt::Checked);
	*/ //MOTOR_SWITCH_OFF
}

//
//: Return velocity sliders to zero when released
void ncm_qcapture2_gui::on_aptXVelocity_sliderReleased()
{
	/*
  ui.aptXVelocity->setValue(0);
  apt_.X()->stop();
	*/ //MOTOR_SWITCH_OFF
}
void ncm_qcapture2_gui::on_aptYVelocity_sliderReleased()
{
	/*
  ui.aptYVelocity->setValue(0);
  apt_.Y()->stop();
	*/ //MOTOR_SWITCH_OFF
}
void ncm_qcapture2_gui::on_aptZVelocity_sliderReleased()
{
	/*
  ui.aptZVelocity->setValue(0);
  apt_.Z()->stop();
	*/ //MOTOR_SWITCH_OFF
}

//
//: Home the motors
void ncm_qcapture2_gui::on_aptXHome_clicked()
{
	/*
  if (!motors_are_live_)
    return;

  // Emit a signal telling apt_server to go home.
  const ncm_qapt_server::apt_axis axis = ncm_qapt_server::AptX;
  apt_.set_max_displacement(axis, 1e9);*/ //MOTOR_SWITCH_OFF
  //emit apt_move_home(axis, /* home_position = */ 12.5f,
  //                         /* return_immediately = */ true);
	//MOTOR_SWITCH_OFF
}
void ncm_qcapture2_gui::on_aptYHome_clicked()
{
	/*
  if (!motors_are_live_)
    return;

  // Emit a signal telling apt_server to go home.
  const ncm_qapt_server::apt_axis axis = ncm_qapt_server::AptY;
  apt_.set_max_displacement(axis, 1e9); */ //MOTOR_SWITCH_OFF
  //emit apt_move_home(axis, /* home_position = */ 12.5f,
  //                         /* return_immediately = */ true);
	//MOTOR_SWITCH_OFF
}
void ncm_qcapture2_gui::on_aptZHome_clicked()
{
	/*
  if (!motors_are_live_)
    return;

  // Emit a signal telling apt_server to go home.
  const ncm_qapt_server::apt_axis axis = ncm_qapt_server::AptZ;
  apt_.set_max_displacement(axis, 1e9);*/ //MOTOR_SWITCH_OFF
  //emit apt_move_home(axis, /* home_position = */ 12.5f,
  //                         /* return_immediately = */ true);
	//MOTOR_SWITCH_OFF
}
 
//: React to mouse clicks on the main graphics view
void ncm_qcapture2_gui::on_graphicsView_clicked(double x, double y)
{
  // Get position of mouse click, relative to centre of sceneview
  double dx = x - 0.5*ui.cameraView1->width();
  double dy = y - 0.5*ui.cameraView1->height();

  set_apt_targetX(dx);
  set_apt_targetY(dy);
}

 
//: React to mouse wheel event on the graphics view
void ncm_qcapture2_gui::on_graphicsView_zoomed(int delta)
{
  // Apply a sensible relative movement in Z.
  // Typically, one 'click' of the mouse wheel will generate a delta of 120.
  // Scale down such that each wheel 'click' moves by 0.1mm
  set_apt_targetZ(static_cast<double>(delta) / 1200.0);
}

void ncm_qcapture2_gui::on_connectToMotors_clicked()
{
  if (apt_timer_.isActive())
  {
    vcl_cout << "stopping" << vcl_endl;
    apt_timer_.stop();
  }
  else
  {
    vcl_cout << "starting" << vcl_endl;
    initialize_apt_timer();
  }
}

void ncm_qcapture2_gui::on_autoFocus_clicked()
{
  if (!motors_are_live_)
    return;

  // Begin the autofocus procedure by searching in the forward direction for
  // a measurable increase in sharpness.
  autoFocus(autoFocus_GetDirectionForward);
}

void ncm_qcapture2_gui::on_calibrate_clicked()
{
  handleCalibrateState();
}
 
//: React to motors completing a move.
void ncm_qcapture2_gui::onAptMoveComplete(ncm_qapt_server::apt_axis)
{
  vcl_cout << "MoveComplete" << vcl_endl;

  if (autoFocus_state_ != autoFocus_Inactive)
    handleAutoFocusState();
  else if (calibrate_state_ != calibrate_Inactive)
    handleCalibrateState();
}
 
//: Handle changes in autofocus_state_.
void ncm_qcapture2_gui::handleAutoFocusState()
{
  switch (autoFocus_state_)
  {
  case autoFocus_Inactive:
    // What are we doing here?
    break;

  case autoFocus_GetDirectionForward:
    // Reverse the direction of search because we've gone far enough that we
    // would have found an improvement had we been going in the right direction.
    /*apt_.set_max_displacement(ncm_qapt_server::AptZ, 1e9);
    autoFocus(autoFocus_GetDirectionReverse);*/ //MOTOR_SWITCH_OFF
    break;

  case autoFocus_GetDirectionReverse:
    // Go to error state since we can't find a direction that generates a 
    // measurable improvement in sharpness. This is likely to happen if 
    // autofocus is started too far from the true position (in the flat region
    // of the function).
    /* apt_.set_max_displacement(ncm_qapt_server::AptZ, 1e9);
    autoFocus(autoFocus_Error); */ //MOTOR_SWITCH_OFF
    break;

  case autoFocus_FindPeak:
    break;

  case autoFocus_GoToPeak:
    break;

  case autoFocus_Error:
    break;

  default:
    break;
  }
}
 
//: Handle changes in calibrate_state_.
void ncm_qcapture2_gui::handleCalibrateState()
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
void ncm_qcapture2_gui::onLowerSharpnessThresholdCrossed()
{/*
  vcl_cout << "LowerThresholdCrossed" << vcl_endl;

  switch (autoFocus_state_)
  {
  case autoFocus_Inactive:
    // What are we doing here?
    break;

  case autoFocus_GetDirectionForward:
    {
    // Reverse the direction of search because we're going the wrong way.

 	  const double current_sharpness = processor_.image_sharpness();
    const double lower_limit = current_sharpness * 0.9;
    const double upper_limit = current_sharpness * 1.1;

    // Need to move the lower threshold a little first, so that it doesn't
    // immediately trigger an event.
    processor_.set_lower_sharpness_threshold(lower_limit);
    autoFocus(autoFocus_GetDirectionReverse);

    vcl_cout << "GetDirectionForward: "
             << "[" << lower_limit << ", " << upper_limit << "]"
             << vcl_endl;

    break;
    }

  case autoFocus_GetDirectionReverse:
    // Go to error state as this shouldn't, in theory, happen: we're in this
    // state because going in the opposite direction decreased sharpness, so 
    // this direction *must* increase sharpness.
    autoFocus(autoFocus_Error);
    break;

  case autoFocus_FindPeak:
    // Stop looking for the peak once the lower threshold has been crossed, as
    // we're now 'over the hump' and getting worse again.
    autoFocus(autoFocus_GoToPeak);
    break;

  case autoFocus_GoToPeak:
    break;

  case autoFocus_Error:
    break;

  default:
    break;
  }*/
}

void ncm_qcapture2_gui::onUpperSharpnessThresholdCrossed()
{/*
  vcl_cout << "UpperThresholdCrossed" << vcl_endl;

  switch (autoFocus_state_)
  {
  case autoFocus_Inactive:
    // What are we doing here?
    break;

  case autoFocus_GetDirectionForward:
    // Fall through
  case autoFocus_GetDirectionReverse:
    // Find the peak, since we know we're now going in the right direction.
    autoFocus(autoFocus_FindPeak);
    break;

  case autoFocus_FindPeak:
    {
    // Store the position associated with the highest sharpness encountered 
    // so far.
    autoFocus_peak_position_ = apt_.Z()->position();

    // Increase both the upper and lower thresholds.
	  const double current_sharpness = processor_.image_sharpness();
    const double lower_limit = current_sharpness * 0.9;
    const double upper_limit = current_sharpness * 1.1;

    processor_.set_lower_sharpness_threshold(lower_limit);
    processor_.set_upper_sharpness_threshold(upper_limit);

    vcl_cout << "FindPeak: " << autoFocus_peak_position_ 
             << "[" << lower_limit << ", " << upper_limit << "]"
             << vcl_endl;

    break;
    }

  case autoFocus_GoToPeak:
    break;

  case autoFocus_Error:
    break;

  default:
    break;
  }*/
}

//Action if compute diff image is clicked
void ncm_qcapture2_gui::on_computeDiffButton_clicked() {
	if (cameras_[0]->is_connected() && cameras_[1]->is_connected()) {
		diff_frame_[0] = true;
		diff_frame_[1] = true;
	}
}

void ncm_qcapture2_gui::on_diffBrightnessSlider_sliderMoved(int value)
{
  int contrast = ui.diffContrastSlider->value();
  int brightness = value;

  // update colormap and redraw
  update_diff_colour_table(brightness, contrast);
  draw_diff_frame();	
}

void ncm_qcapture2_gui::on_diffContrastSlider_sliderMoved(int value)
{
  int contrast = value;
  int brightness = ui.diffBrightnessSlider->value();

  // update colormap and redraw
  update_diff_colour_table(brightness, contrast);
  draw_diff_frame();	
}
 
//
//  Private methods
//

//Load camera
bool ncm_qcapture2_gui::load_camera(int cam_num)
{
	//Try and connect the device
	if ( !cameras_[ cam_num ]->connect_device() ) {
		vcl_cout << "Camera failed to load" << vcl_endl;
		return false;
	}
	//Make this camera the current camera
	change_active_camera( cam_num );

	//Set this camera to be saved from
	save_from_camera_[ cam_num ] = true;

	// clear old scene
	raw_pixmap_items_[ cam_num ] = 0;
	raw_pixmap_item_d_ = 0;
	scenes_[ cam_num ]->clear();
	scene_d_.clear();

	// Use highest frame rate by default
  ui.frameRateSelection->setCurrentIndex(0);

	//-----------------------------------------------------------------------------

	//Start the device's live mode
	cameras_[ cam_num ]->start_live();
	ui.stopLiveButton->setText( QString( "Stop live" ) );
	ui.stopLiveButton->setEnabled( true );

	// draw scene and fit image in view
  redraw_scene( cam_num );
	if (cam_num) //(==1)
		ui.cameraView2->fitInView(scenes_[cam_num]->sceneRect(),Qt::KeepAspectRatio);
	else
		ui.cameraView1->fitInView(scenes_[cam_num]->sceneRect(),Qt::KeepAspectRatio);

	//Enable the camera and save sequence controls
	ui.tabDisplay->setEnabled( true );
	ui.tabCamera->setEnabled( true );
	ui.tabCapture->setEnabled( true );

	return true;
}

//: Create the status bar's labels and progress bar
void ncm_qcapture2_gui::initializeStatusbar()
{
  statusBar_main_.setFrameStyle(QFrame::Panel & QFrame::Sunken);
  statusBar_main_.setLineWidth(1);
  statusBar_main_.setText("");
  statusBar_main_.setContentsMargins(4, 0, 4, 0);
  statusBar()->addWidget(&statusBar_main_, 1);

  statusBar_progress_.setFixedWidth(160);
  statusBar()->addWidget(&statusBar_progress_);
}

//
//: Update status bar.
//  A progress value in [0..1] sets the progress bar value.
//  progress < 0 makes no change.
void ncm_qcapture2_gui::updateStatus(
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
void ncm_qcapture2_gui::updateProgress(
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
void ncm_qcapture2_gui::initializeAlignerThread()
{
  aligner_thread_.setObjectName("aligner_thread");
  aligner_.set_n_levels(2);

  aligner_.moveToThread(&aligner_thread_);

  // Register data types before they are used in signals
  qRegisterMetaType< vcl_vector<vcl_string> >("vcl_vector<vcl_string>");

  //// This -> aligner
  //QObject::connect( this, SIGNAL(alignSequence(vcl_vector<vcl_string>)),
  //                  &aligner_, SLOT(alignSequence(vcl_vector<vcl_string>)) ); 

  //// Aligner -> this
  //QObject::connect( &aligner_, SIGNAL(framesAligned(int, int)),
  //                  this, SLOT(onFramesAligned(int, int)) ); 
  //QObject::connect( &aligner_, SIGNAL(alignmentFinished()),
  //                  this, SLOT(onAlignmentFinished()) ); 

  //aligner_thread_.start();
}

 
//: Add the qimage to the canvas.
void ncm_qcapture2_gui::draw_pixmap(QImage* qimage, int cam_num)
{
  if ( (qimage == NULL) || 
       (qimage->isNull()) )
  {
    return;
  }
	qimage->setColorTable(raw_colour_tables_[ cam_num ]);

  if (raw_pixmap_items_[ cam_num ] == NULL)
		raw_pixmap_items_[ cam_num ] = scenes_[ cam_num ]->addPixmap( QPixmap::fromImage( *qimage ));
  else
		raw_pixmap_items_[ cam_num ]->setPixmap(QPixmap::fromImage( *qimage ));

	scenes_[ cam_num ]->setSceneRect(QRectF(0,0, qimage->width(), qimage->height() ));

}

 
//: Add crosshairs to scene.
void ncm_qcapture2_gui::draw_crosshairs(int cam_num)
{
  QPen crosshair_pen(QColor(0, 255, 0));

  scenes_[ cam_num ]->addLine(scenes_[ cam_num ]->width()*0.5, 0, scenes_[ cam_num ]->width()*0.5, scenes_[ cam_num ]->height(), 
                 crosshair_pen);
  scenes_[ cam_num ]->addLine(0, scenes_[ cam_num ]->height()*0.5, scenes_[ cam_num ]->width(), scenes_[ cam_num ]->height()*0.5,
                 crosshair_pen);
}

//
//: Write sequence properties to a text file.
void ncm_qcapture2_gui::write_properties( vcl_ostream& tfs )
{
	tfs << vsl_indent() << "NCM QCapture Sequence Properties" << vcl_endl;
  vsl_indent_inc(tfs);

	tfs << vsl_indent() << "Date: " << QDate::currentDate().toString( "ddd dd MMM yyyy" ).toStdString() << vcl_endl;
	tfs << vsl_indent() << "Time: " << QTime::currentTime().toString( "hh:mm:ss" ).toStdString() << vcl_endl;
	tfs << vcl_endl;
	tfs << vsl_indent() << "Patient name: " << ui.patientNameTextEdit->toPlainText().toStdString() << vcl_endl;
	tfs << vsl_indent() << "Sequence name: " << ui.sequenceNameTextEdit->toPlainText().toStdString() << vcl_endl;
	tfs << vsl_indent() << "Sequence length: " << ui.sequenceLengthSpinBox->value() << " s" << vcl_endl;
	tfs << vcl_endl;

	for (int cam_num = 0; cam_num < 2; cam_num++) {
		if ( save_from_camera_[ cam_num ] ) {
			write_camera_properties( tfs , cam_num);
		}
		else {
			tfs << vsl_indent() << "Camera " << cam_num+1 << ", not active" << vcl_endl;
		}
	}
  tfs << vcl_flush;
}

void ncm_qcapture2_gui::write_camera_properties( vcl_ostream& tfs , int cam_num)
{
	tfs << vsl_indent() << "Camera " << cam_num+1 << ", Type: " << cameras_[cam_num]->get_camera_type() << vcl_endl;
	tfs << vsl_indent() << "Camera ID: " << cameras_[cam_num]->get_camera_id() << vcl_endl;
	vsl_indent_inc(tfs);
	tfs << vsl_indent() << "Frame rate: " << cameras_[cam_num]->get_current_FPS() << " fps" << vcl_endl;
	tfs << vsl_indent() << "Exposure time: ";
	if ( cameras_[cam_num]->is_auto_exposure() )
		tfs << "auto" << vcl_endl;
	else
		tfs << cameras_[cam_num]->get_current_exposure() << vcl_endl;

	tfs << vsl_indent() << "Gain: ";
	if ( cameras_[cam_num]->is_auto_gain() )
		tfs << "auto" << vcl_endl;
	else
		tfs << cameras_[cam_num]->get_current_gain() << vcl_endl;

	tfs << vsl_indent() << "Orientation: H flipped " << flipped_h_[ cam_num ] << ", V flipped " << flipped_v_[ cam_num ] << vcl_endl;
}

void ncm_qcapture2_gui::update_framerate_controls()
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

	if ( cameras_[current_camera_]->get_available_FPS( available_FPS ) &&
       cameras_[current_camera_]->get_current_FPS( fps ) )
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

void ncm_qcapture2_gui::update_exposure_controls()
{
	//Get exposure range of camera and set spin box and slider
	double min_exposure, max_exposure;
	if ( cameras_[current_camera_]->get_exposure_range( min_exposure, max_exposure ) )
	{
		ui.exposureSpinBox->setRange( 0, 100 );
		ui.exposureSlider->setRange( 0, 100 );
		ui.exposureSpinBox->setExposureMin( min_exposure );
		ui.exposureSpinBox->setExposureMax( max_exposure );
	}
}

void ncm_qcapture2_gui::update_gain_controls()
{
	//Get gain range of camera and set spin box and slider	
	double min_gain, max_gain;
	if ( cameras_[current_camera_]->get_gain_range( min_gain, max_gain ) )
	{
		ui.gainSpinBox->setRange( min_gain, max_gain );
		ui.gainSlider->setRange( int( min_gain ), int( max_gain ) );
	}
}

//
//: Update controller comboboxes that enable you to assign controllers to axes.
void ncm_qcapture2_gui::update_apt_comboboxes()
{
	/*
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
	*/ //MOTOR_SWITCH_OFF
}

//
//: Set the target position for X axis.
void ncm_qcapture2_gui::set_apt_targetX(double distance_pixels)
{
	/*
  const double mm_per_pixel = 1.0 / pixels_per_mm_;
  double distance_mm = mm_per_pixel * distance_pixels;

  if (motors_are_live_)
  {
    apt_.X()->set_velocity(1.0f);
    apt_.X()->move_by(distance_mm);
  }
	*/ //MOTOR_SWITCH_OFF
}

void ncm_qcapture2_gui::set_apt_targetY(double distance_pixels)
{
	/*
  const double mm_per_pixel = 1.0 / pixels_per_mm_;
  double distance_mm = mm_per_pixel * distance_pixels;
  
  if (motors_are_live_)
  {
    apt_.Y()->set_velocity(1.0f);
    apt_.Y()->move_by(distance_mm);
  }
	*/ //MOTOR_SWITCH_OFF
}

void ncm_qcapture2_gui::set_apt_targetZ(double z)
{
	/*
  if (motors_are_live_)
  {
    apt_.Z()->set_velocity(1.0f);
    apt_.Z()->move_by(z);
  }
	*/ //MOTOR_SWITCH_OFF
}

//
//: Autofocus procedure, implemented as a finite state machine.
void ncm_qcapture2_gui::autoFocus(autoFocus_state_enum state)
{
	/*
  autoFocus_state_ = state;

  switch (state)
  {
  case autoFocus_GetDirectionForward:
    {
    vcl_cout << "GetDirectionForward" << vcl_endl;
    // Determine the direction of movement for maximum sharpness

	  const double reference_sharpness = processor_.image_sharpness();
    const double lower_limit = reference_sharpness * 0.9;
    const double upper_limit = reference_sharpness * 1.1;

    // Move camera until sharpness changes by some (absolute or relative) 
    // amount, or until a fixed distance (e.g. 1mm) has been traversed.
    processor_.set_lower_sharpness_threshold(lower_limit);
    processor_.set_upper_sharpness_threshold(upper_limit);

    vcl_cout << "GetDirectionForward: ["
             << lower_limit << ", " << upper_limit << "]"
             << vcl_endl;

    // Set limits on displacement from current position.
    apt_.set_reference_position(ncm_qapt_server::AptZ);
    apt_.set_max_displacement(ncm_qapt_server::AptZ, 1.0f);

    // Start searching.
    apt_.Z()->move_at_velocity(0.25f);
    break;
    }

  case autoFocus_GetDirectionReverse:
    vcl_cout << "GetDirectionReverse" << vcl_endl;
    // Reverse the direction of the search if there was no improvement in the 
    // forward direction.

    // Reuse previous thresholds instead of defining new ones.

    // Move camera until sharpness changes by some (absolute or relative) 
    // amount, or until a fixed distance has been traversed.
    // Set limits on displacement from current position.
    apt_.set_reference_position(ncm_qapt_server::AptZ);
    apt_.set_max_displacement(ncm_qapt_server::AptZ, 2.0f);

    // Start searching.
    apt_.Z()->move_at_velocity(-0.25f);
    break;

  case autoFocus_FindPeak:
    vcl_cout << "FindPeak" << vcl_endl;
    // Maintain the existing velocity, updating the optimal position every
    // time the upper threshold is exceeded.

    // When the peak has been passed, the lower threshold will be crossed
    // causing a transition in the finite state machine.

    // No need to do anything here - it's all handled by signals and slots.
    apt_.set_max_displacement(ncm_qapt_server::AptZ, 1e9);
    break;

  case autoFocus_GoToPeak:
    vcl_cout << "GoToPeak" << vcl_endl;
    // Reverse direction and go back more slowly until a maximum is found.

    processor_.set_lower_sharpness_threshold(-1e9);
    processor_.set_upper_sharpness_threshold(1e9);

    vcl_cout << "Peak pos = " << autoFocus_peak_position_ << vcl_endl;
    apt_.Z()->move_to(autoFocus_peak_position_, true);
    //emit apt_move_to(ncm_qapt_server::AptZ, autoFocus_peak_position_, false);

    autoFocus_state_ = autoFocus_Inactive;
    break;

  default:
    break;
  }*/
}

//
//: Initialize the registry accessor and default parameters if they don't
//  already exist.
void ncm_qcapture2_gui::initialize_registry_accessor()
{
  ra_.get_current_user_key();

  long result = 
      ra_.get_subkey("Software\\University of Manchester\\NCM QCapture");

  if (result != ERROR_SUCCESS)
    ra_.create_subkey("Software\\University of Manchester\\NCM QCapture");
}

//
//: Initialize the saver class with defaults
void ncm_qcapture2_gui::initialize_saver()
{
  // Get the root dir from the Registry (if it exists).
  vcl_string root_dir;
  long result = ra_.read_string("RootDir", root_dir);
  if (result != ERROR_SUCCESS)
  {
    root_dir = "C:/isbe/nailfold/playground/dmk/capture_software";
    ra_.write_string("RootDir", root_dir);
  }
	ui.saveDirTextEdit->setText( root_dir.c_str() );
	ui.saveDirTextEdit->setReadOnly( true );
	savers_[0]->setRootDir( root_dir );
	savers_[1]->setRootDir( root_dir );

  // This will change with every session so don't store in the Registry.
	QString sequence_name = QDate::currentDate().toString( "yyyy_MM_dd" );
	ui.sequenceNameTextEdit->setText( sequence_name );
	savers_[0]->setSequenceName( sequence_name.toStdString() );
	savers_[1]->setSequenceName( sequence_name.toStdString() );

  // Get the frame prefix from the Registry (if it exists).
  vcl_string patient_name;
  result = ra_.read_string("PatientName", patient_name);
  if (result != ERROR_SUCCESS)
  {
    patient_name = "PatientX";
    ra_.write_string("PatientName", patient_name);
  }
	ui.patientNameTextEdit->setText( patient_name.c_str() );
	savers_[0]->setPatientName( patient_name );
	savers_[1]->setPatientName( patient_name );

	//Make sure the other saver parameters match the gui
	QString hand = ui.handComboBox->currentText();
	savers_[0]->setHand( hand.toStdString() );
	savers_[1]->setHand( hand.toStdString() );
	
	QString digit = ui.digitComboBox->currentText();
	savers_[0]->setDigit( digit.toStdString() );
	savers_[1]->setDigit( digit.toStdString() );

	QString mag = ui.magComboBox->currentText();
	savers_[0]->setMag( mag.toStdString() );
	savers_[1]->setMag( mag.toStdString() );
}

//
//: Initialize motor position spinboxes
void ncm_qcapture2_gui::initialize_apt_spinboxes()
{
	/*
  ncm_apt_controller_base* c = NULL;
  QDoubleSpinBox* sp = NULL;

  c = apt_.X();
  sp = ui.aptXPosition;
  //sp->setMinimum(c->min_position());
  //sp->setMaximum(c->max_position());
  //sp->setValue(c->position());
  //sp->setSingleStep(x_precision_);
  //sp->setDecimals(-vcl_log10(x_precision_));

  c = apt_.Y();
  sp = ui.aptYPosition;
  //sp->setMinimum(c->min_position());
  //sp->setMaximum(c->max_position());

  c = apt_.Z();
  sp = ui.aptZPosition;
  //sp->setMinimum(c->min_position());
  //sp->setMaximum(c->max_position());
	*/ //MOTOR_SWITCH_OFF
}

//
//: Initialize motor control timer
void ncm_qcapture2_gui::initialize_apt_timer()
{/*
  apt_timer_.start(50);
  apt_thread_.start();

  // Make a note of the current position and set a threshold upon which to
  // emit a signal.
  apt_.set_reference_position(ncm_qapt_server::AptZ);
  apt_.set_max_displacement(ncm_qapt_server::AptZ, 1.0f);*/
}

//
//: Initialize the processor thread
void ncm_qcapture2_gui::initialize_processor_thread()
{
  processor_thread_.setObjectName("processor_thread");

	processor_.moveToThread(&processor_thread_);
	
  // Connect signals to slots that will run in the thread


  // Camera tells processor there's a frame ready
  QObject::connect( cameras_[0], SIGNAL(frame_ready()),
                    this, SLOT(tag_frame1()) );
	QObject::connect( cameras_[1], SIGNAL(frame_ready()),
                    this, SLOT(tag_frame2()) );

	QObject::connect( this, SIGNAL(frame_tagged(int, bool, bool, bool)),
                    &processor_, SLOT(process_frame(int, bool, bool, bool)) );

	//Processor tells GUI there's a frame ready to draw
  QObject::connect( &processor_, SIGNAL(frame_to_draw( int )),
                    this, SLOT(redraw_scene( int )) );

  //Processor tells GUI the saver is done
  QObject::connect( &processor_, SIGNAL(processing_finished()),
						        this, SLOT(on_processing_finished()) );

	//Processor tells GUI a diff frame is ready to draw
  QObject::connect( &processor_, SIGNAL(diff_to_draw()),
						        this, SLOT(draw_diff_frame()) );

	processor_thread_.start();
}

//
//: Initialize the saver thread
void ncm_qcapture2_gui::initialize_saver_thread()
{
  saver_thread_.setObjectName("saver_thread");
	savers_[0]->moveToThread(&saver_thread_);
	savers_[1]->moveToThread(&saver_thread_);

	savers_[0]->setCameraSuffix( ui.camera1SuffixTextEdit->toPlainText().toStdString() );
	savers_[1]->setCameraSuffix( ui.camera2SuffixTextEdit->toPlainText().toStdString() );

  // Connect signals to slots that will run in the thread

  // Processor tells saver to save a frame
  QObject::connect( &processor_, SIGNAL(frame_to_save1(int, int)),
                    savers_[0], SLOT(save_frame(int, int)) );
	QObject::connect( &processor_, SIGNAL(frame_to_save2(int, int)),
                    savers_[1], SLOT(save_frame(int, int)) );

	//Saver tells main gui it has saved a frame
  QObject::connect( savers_[0], SIGNAL(frame_saved(int, int)),
                    this, SLOT(on_frame_saved(int, int)) );
	QObject::connect( savers_[1], SIGNAL(frame_saved(int, int)),
                    this, SLOT(on_frame_saved(int, int)) );

	//Saver tells main gui it has finished saving
  QObject::connect( savers_[0], SIGNAL(saving_finished(int)),
                    this, SLOT(on_saving_finished(int)) );	 
  QObject::connect( savers_[1], SIGNAL(saving_finished(int)),
                    this, SLOT(on_saving_finished(int)) );

  saver_thread_.start();
}

void ncm_qcapture2_gui::get_apt_registry_values()
{
  long result = ERROR_SUCCESS;

  // Get motor assignments.
  int id = 0;
  int item = -1;
  
  result = ra_.read_numeric("AptXId", id);
  if (result == ERROR_SUCCESS)
  {
    item = ui.aptXCombobox->findText(QString::number(id));
    if (item != -1)
      ui.aptXCombobox->setCurrentIndex(item);
  }

  ra_.read_numeric("AptYId", id);
  if (result == ERROR_SUCCESS)
  {
    item = ui.aptYCombobox->findText(QString::number(id));
    if (item != -1)
      ui.aptYCombobox->setCurrentIndex(item);
  }

  ra_.read_numeric("AptZId", id);
  if (result == ERROR_SUCCESS)
  {
    item = ui.aptZCombobox->findText(QString::number(id));
    if (item != -1)
      ui.aptZCombobox->setCurrentIndex(item);
  }

  // Find out whether axes were reversed.
  bool reversed = false;

  result = ra_.read_boolean("XReversed", reversed);
  if (result == ERROR_SUCCESS)
    ui.aptXReverse->setChecked(reversed);

  result = ra_.read_boolean("YReversed", reversed);
  if (result == ERROR_SUCCESS)
    ui.aptYReverse->setChecked(reversed);

  result = ra_.read_boolean("ZReversed", reversed);
  if (result == ERROR_SUCCESS)
    ui.aptZReverse->setChecked(reversed);
}
 
//
//: Initialize the APT (motor) thread
void ncm_qcapture2_gui::initialize_apt_thread()
{
  // Read previous values from the registry
  get_apt_registry_values();

  //apt_thread_.setObjectName("apt_thread");

  //apt_.moveToThread(&apt_thread_);
  //apt_timer_.moveToThread(&apt_thread_);

  // Starting the timer does nothing until the thread is started.
  //apt_timer_.start(/* interval = */ 50);

	/*
  // Register the enum so that the correct metafunctions can be called.
  qRegisterMetaType<ncm_qapt_server::apt_axis>("ncm_qapt_server::apt_axis");

  // Connect signals to slots that will run in the thread
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

  // Don't start the thread yet.
	*/ //MOTOR_SWITCH_OFF
}

//
//: Connect signals to slots (not surprisingly...)
void ncm_qcapture2_gui::connect_signals_to_slots()
{
  // Signals that trigger slots in the main thread

	/* Don't worry about the autofocus stuff for now (it won't be used in the two camera setup)
  QObject::connect( &apt_timer_, SIGNAL(timeout()),
                    this, SLOT(onAptTimer()) );

  QObject::connect( &apt_, SIGNAL(move_complete(ncm_qapt_server::apt_axis)),
                    this, SLOT(onAptMoveComplete(ncm_qapt_server::apt_axis)) );

  QObject::connect( &processor_, SIGNAL(lower_sharpness_exceeded()),
                    this, SLOT(onLowerSharpnessThresholdCrossed()) );
  QObject::connect( &processor_, SIGNAL(upper_sharpness_exceeded()),
                    this, SLOT(onUpperSharpnessThresholdCrossed()) );
										*/
}

