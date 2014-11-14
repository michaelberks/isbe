#include "ncm_qcamera_gui.h"


ncm_qcamera_gui::ncm_qcamera_gui()
    : contrast_min_(0),
      contrast_max_(255)
{
}

void ncm_qcamera_gui::set_ui_controls(
		QPushButton *_loadCameraButton,
		QGroupBox *_cameraControls,
		QSlider *_exposureSlider,
		QExposureSpinBox *_exposureSpinBox,
		QCheckBox *_autoExposure,
		QSlider *_gainSlider,
		QDoubleSpinBox *_gainSpinBox,
		QCheckBox *_autoGain, 
		QComboBox *_frameRateSelection,
		QGroupBox *_displayControls,
		QSlider *_brightnessSlider,
		QSlider *_contrastSlider,
		QPushButton *_stopLiveButton )
{
	loadCameraButton_ = _loadCameraButton;
	cameraControls_ = _cameraControls;
	exposureSlider_ = _exposureSlider;
	exposureSpinBox_ = _exposureSpinBox;
	autoExposure_ = _autoExposure;
	gainSlider_ = _gainSlider;
	gainSpinBox_ = _gainSpinBox;
	autoGain_ = _autoGain; 
	frameRateSelection_ = _frameRateSelection;
	displayControls_ = _displayControls;
	brightnessSlider_ = _brightnessSlider;
	contrastSlider_ = _contrastSlider;
	stopLiveButton_ = _stopLiveButton;
}

void ncm_qcamera_gui::set_camera( ncm_camera *_camera )
{
	camera_ = _camera;
}

void ncm_qcamera_gui::set_scene_controls(
									nailfold_graphics_view *_graphicsView,
									QGraphicsPixmapItem *_raw_pixmap_item,
									QVector<QRgb> *_raw_colour_table,
									QGraphicsScene *_scene,
									QImage *_qimage )
{
	graphicsView_ = _graphicsView;
	raw_pixmap_item_ = _raw_pixmap_item;
	raw_colour_table_ = _raw_colour_table;
	scene_ = _scene;
	qimage_ = _qimage;

	// tell graphics view to use this scene
	graphicsView_->setScene( scene_ );

  // define default colormaps
	update_raw_colour_table();
}

//--------------------------------------------------------------------------------
//UI callback controls
//--------------------------------------------------------------------------------
void ncm_qcamera_gui::loadCameraButton_clicked()
{
	//Check if device is already connected and if so, tidy up
	

	//Try and connect the device
	if ( !camera_->connect_device() ) {
		vcl_cout << "Camera failed to load" << vcl_endl;
		return;
	}

	// clear old scene
	raw_pixmap_item_ = 0;
	scene_->clear();

	//---------------------------------------------------------------------
	// Get camera property ranges to fill control widget values

	//Get available frame rates and update frame selection combo box
	vcl_vector<double> available_FPS;
	double fps;
	if ( camera_->get_available_FPS( available_FPS ) & camera_->get_current_FPS( fps ) )
	{
		QString text_double;
		for (int i = 0; i < available_FPS.size(); i++)
		{
			double fps_i = available_FPS[ i ];
			text_double.setNum( fps_i );
			frameRateSelection_->addItem( text_double );
		
			if ( fps_i == fps )
				frameRateSelection_->setCurrentIndex( i );
		}
	}

	//Get gain range of camera and set spin box and slider	
	double min_gain, max_gain;
	if ( camera_->get_gain_range( min_gain, max_gain ) )
	{
		gainSpinBox_->setRange( min_gain, max_gain );
		gainSlider_->setRange( int( min_gain ), int( max_gain ) );
	}

	//Get exposure range of camera and set spin box and slider
	double min_exposure, max_exposure;
	if ( camera_->get_exposure_range( min_exposure, max_exposure ) )
	{
		exposureSpinBox_->setRange( 0, 100 );
		exposureSlider_->setRange( 0, 100 );
		exposureSpinBox_->setExposureMin( min_exposure );
		exposureSpinBox_->setExposureMax( max_exposure );
	}

	//-----------------------------------------------------------------------------

	//Start the device's live mode
	camera_->start_live();

	// draw scene and fit image in view
  redraw_scene();
  graphicsView_->fitInView(scene_->sceneRect(),Qt::KeepAspectRatio);

	//Enable the camera and save sequence controls
	displayControls_->setEnabled( true );
	brightnessSlider_->setEnabled( true );
	contrastSlider_->setEnabled( true );

	cameraControls_->setEnabled( true );
	gainSlider_->setEnabled( false ); //Start in auto mode
	exposureSlider_->setEnabled( false ); //Start in auto mode
	gainSpinBox_->setEnabled( false ); //Start in auto mode
	exposureSpinBox_->setEnabled( false ); //Start in auto mode
}

void ncm_qcamera_gui::stopLiveButton_clicked()
{
	camera_->stop_live();
}

//----------------------------------------------------------------------------
// Display control callbacks
void ncm_qcamera_gui::brightnessSliderMoved(int value)
{
  int contrast = contrastSlider_->value();
  int brightness = value;
  update_display(contrast_min_, contrast_max_, contrast, brightness);
	
}

void ncm_qcamera_gui::contrastSliderMoved(int value)
{
  int contrast = value;
  int brightness = brightnessSlider_->value();
	update_display(contrast_min_, contrast_max_, contrast, brightness);
  
}

//----------------------------------------------------------------------------
// Camera control callbacks
void ncm_qcamera_gui::frameRateSelectionActivated ( const QString &text ) 
{
	double frames_per_sec = text.toDouble();
	vcl_cout << "New FPS = " << frames_per_sec << vcl_endl;
	
	//Stop live mode, change the FPS, then restart live mode
	camera_->set_FPS( frames_per_sec );
}

void ncm_qcamera_gui::gainSpinBoxChanged( double value )
{
	gainSlider_->setValue( int( value ) );
	camera_->set_gain( value );

}
void ncm_qcamera_gui::gainSliderMoved( int value )
{
	gainSpinBox_->setValue( double( value ) );
	camera_->set_gain( double( value ) );
}
void ncm_qcamera_gui::autoGainChanged( int state )
{
	if ( state )
	{
		//disable slider and spin box
		gainSlider_->setEnabled( false );
		gainSpinBox_->setEnabled( false );

		//switch on camera's auto gain
		camera_->switch_auto_gain( true );

		vcl_cout << "Auto gain switched on" << vcl_endl;
	}
	else
	{
		//switch off camera's auto gain and get current gain value
		double gain;
		camera_->switch_auto_gain( false );
		camera_->get_current_gain( gain );

		//enable gain slider and spin box and set to current auto value
		gainSlider_->setEnabled( true );
		gainSpinBox_->setEnabled( true );
		gainSpinBox_->setValue( gain );
		gainSlider_->setValue( int( gain ) );
		vcl_cout << "Auto gain switched off" << vcl_endl;
	}
}
void ncm_qcamera_gui::exposureSpinBoxChanged( double value )
{
	double exposure = exposureSpinBox_->exposureFromValue( value );
	exposureSlider_->setValue( int( value ) );
	camera_->set_exposure( exposure );
}
void ncm_qcamera_gui::exposureSliderMoved( int value )
{
	double exposure = exposureSpinBox_->exposureFromValue( value );
	exposureSpinBox_->setValue(double( value ));
	camera_->set_exposure( exposure );
}
void ncm_qcamera_gui::autoExposureChanged( int state )
{
	if ( state )
	{
		//disable slider and spin box
		exposureSlider_->setEnabled( false );
		exposureSpinBox_->setEnabled( false );

		//switch on camera's auto exposure
		camera_->switch_auto_exposure( true );

		vcl_cout << "Auto exposure switched on" << vcl_endl;
	}
	else
	{
		//switch off camera's auto exposure and get current exposure value
		double exposure, value;
		camera_->switch_auto_exposure( false );
		camera_->get_current_exposure( exposure );
		value = exposureSpinBox_->valueFromExposure( exposure );

		//enable exposure slider and spin box and set to current auto value
		exposureSlider_->setEnabled( true );
		exposureSpinBox_->setEnabled( true );
		exposureSpinBox_->setValue( value );
		exposureSlider_->setValue( int( value ) );
		vcl_cout << "Auto exposure switched off" << vcl_endl;
		vcl_cout << "Current exposure: " << exposure << vcl_endl;
	}
}

void ncm_qcamera_gui::update_raw_colour_table()
{
  raw_colour_table_->resize(256);

  // fill everything from first to constrast_min_ with 0 (black)
  for (int i = 0; i < contrast_min_; i++)
    raw_colour_table_->replace(i, qRgb(0,0,0) );

  // fill everything from constrast_max_ to end with 255 (white)
  for (int i = contrast_max_; i < 256; i++)
    raw_colour_table_->replace(i, qRgb(255,255,255) );

  // fill values in between with a linear gradient
  int denom = contrast_max_ - contrast_min_;
  for (int i = contrast_min_; i < contrast_max_; i++)
  {
    int grey = (255 * (i - contrast_min_)) / denom;
    raw_colour_table_->replace(i, qRgb(grey,grey,grey) );
  }
}

//Private helper functions
//: Redraw items on canvas
void ncm_qcamera_gui::redraw_scene()
{
	qimage_->setColorTable(*raw_colour_table_);

  if (raw_pixmap_item_ == 0)
		raw_pixmap_item_ = scene_->addPixmap( QPixmap::fromImage( *qimage_ ));
  else
		raw_pixmap_item_->setPixmap(QPixmap::fromImage( *qimage_ ));

	scene_->setSceneRect(QRectF(0,0, qimage_->width(), qimage_->height() ));
	scene_->update();
}

void ncm_qcamera_gui::update_display(int & contrast_min, int & contrast_max, int contrast, int brightness)
{
	int c_min = (255-brightness) - (255-contrast)/2;
  int c_max = (255-brightness) + (255-contrast)/2;

	
	if (c_min < 0)
    contrast_min = 0;
  else if (c_min > 255)
    contrast_min = 255;
  else
    contrast_min = c_min;

  if (c_max < 0)
    contrast_max = 0;
  else if (c_max > 255)
    contrast_max = 255;
  else
    contrast_max = c_max;
	
  // update colormap and redraw
  update_raw_colour_table();
  redraw_scene();
}