#include "ncm_dmk_camera.h"

//: Define static class members here
// ncm_dmk_camera::member_ = 0;

//: Define static const class members here
// ncm_dmk_camera::const_member_;

//: Define class member functions


ncm_dmk_camera::ncm_dmk_camera()
{
	// Set the sink 
	dmk_sink_ = FrameHandlerSink::create( &frame_filter_, DShowLib::eY800, 0 );
	
	//dmk_sink_->setSnapMode( true ); //start in snap mode so a single frame can be snapped
	dmk_sink_->setSnapMode( false ); //start in grab mode to stream from device
	dmk_grabber_.setSinkType( dmk_sink_ );
}

ncm_dmk_camera::~ncm_dmk_camera()
{
}

void ncm_dmk_camera::capture_image(vil_image_view<vxl_byte>& dest)
{
}

bool ncm_dmk_camera::is_connected( )
{
	return dmk_grabber_.isDevValid();
}
bool ncm_dmk_camera::connect_device( )
{
	return connect_device(true);
}
bool ncm_dmk_camera::connect_device(bool use_existing )
{
	//Check if device already connected, if so, disconnect
	if ( is_connected( ) )
		disconnect_device( );

	//First try loading from existing device page
	if (!use_existing || !dmk_grabber_.loadDeviceStateFromFile("device.xml", true))
	{
		//if that didn't work, open a device load page
		//Display the device load page
		if( !dmk_grabber_.showDevicePage() )
		{
			vcl_cout << "Device page wouldn't load: " 
							 << dmk_grabber_.getLastError().c_str() 
							 << vcl_endl;
			return 0;
		}
	
		// If we have selected a valid device, save it to the file "device.xml", so
		// the application can load it automatically when it is started the next time.
		if( dmk_grabber_.isDevValid())
			dmk_grabber_.saveDeviceStateToFile("device.xml");
		else
		{
			vcl_cout << "Couldn't find any devices" << vcl_endl;
			return 0;
		}	
	}

	//Connect the listener frame_ready signal to the camera's frame_ready signal
  QObject::connect( &frame_filter_, SIGNAL( frame_ready() ),
		                this, SIGNAL( frame_ready() ));

	//Get the available VCD properties for the device
	pVCD_props_ = dmk_grabber_.getAvailableVCDProperties();

	return true;
}

bool ncm_dmk_camera::disconnect_device( )
{
	// If live video is running, stop it.
	if(dmk_grabber_.isLive())
		dmk_grabber_.stopLive();

	// Disconnect the listener frame_ready signal from the camera's frame_ready 
  // signal
  QObject::disconnect(&frame_filter_, SIGNAL( frame_ready() ),
							        this, SIGNAL( frame_ready() ));

	return 0;
}

bool ncm_dmk_camera::start_live( )
{
	return dmk_grabber_.startLive( /* showVideoWindow = */ false );
}

bool ncm_dmk_camera::stop_live( )
{
	return dmk_grabber_.stopLive();
}

bool ncm_dmk_camera::is_live( )
{
	return dmk_grabber_.isLive();
}

vcl_string ncm_dmk_camera::get_camera_type()
{
	return "Imaging Source DMK USB camera";
}

vcl_string ncm_dmk_camera::get_camera_id()
{
	//See also getBaseName(), getName()
	return dmk_grabber_.getDev().getUniqueName();
}

bool ncm_dmk_camera::get_available_FPS(vcl_vector<double>& available_FPS)
{
	Grabber::tFPSListPtr FPS_list = dmk_grabber_.getAvailableFPS();

	int n = FPS_list->size();

	if (n <= 0) {
    available_FPS.resize(0);
		return false;
	}
  else
  {
    if (available_FPS.size() != n) 
      available_FPS.resize(n);

		for (int i = 0; i < n; i++)
			available_FPS[i] = FPS_list.operator *() [ i ];

		return true;
  }
}

bool ncm_dmk_camera::get_exposure_range(double& min_exposure, double& max_exposure)
{
	// Retrieve the absolute value interface for exposure.
	tIVCDAbsoluteValuePropertyPtr pAbsVal = 0;    

	if( pVCD_props_->findInterfacePtr( VCDID_Exposure, VCDElement_Value, pAbsVal ) != NULL )
	{
		min_exposure = pAbsVal->getRangeMin();
		max_exposure = pAbsVal->getRangeMax();
		vcl_cout << "Exposure range: [" 
             << min_exposure << ", " << max_exposure << "]" 
             << vcl_endl;
		return true;
	}
	else
	{
		vcl_cout << "DMK: Couldn't find exposure interface" << vcl_endl;
		return false;
	}
}

bool ncm_dmk_camera::get_gain_range(double& min_gain, double& max_gain)
{
	tsPropertyRange gain_range = 
      dmk_grabber_.getPropertyRange( VideoProcAmp_Gain );

	min_gain = double( gain_range.min );
	max_gain = double( gain_range.max );

	return ((0 < min_gain) && (min_gain <= max_gain));
}

bool ncm_dmk_camera::get_current_FPS(double& curr_FPS)
{
	curr_FPS = dmk_grabber_.getFPS();
	return curr_FPS > 0;
}

bool ncm_dmk_camera::get_current_exposure(double& curr_exposure)
{
	// Retrieve the absolute value interface for exposure.
	tIVCDAbsoluteValuePropertyPtr pAbsVal = 0;    

	if( pVCD_props_->findInterfacePtr( VCDID_Exposure, VCDElement_Value, pAbsVal ) != 0 )
	{
		pAbsVal->get_Value( &curr_exposure );
		return true;
	}
	else
	{
		vcl_cout << "DMK: Couldn't find exposure interface" << vcl_endl;
		return false;
	}
}

bool ncm_dmk_camera::get_current_gain(double& curr_gain)
{
	curr_gain = dmk_grabber_.getProperty( VideoProcAmp_Gain );
	return true;
}

double ncm_dmk_camera::get_current_FPS()
{
	double fps = 0;
	get_current_FPS( fps );
	return fps;
}
double ncm_dmk_camera::get_current_exposure()
{
	double exposure = 0;
	get_current_exposure( exposure );
	return exposure;
}
double ncm_dmk_camera::get_current_gain()
{
	double gain = 0;
	get_current_gain( gain );
	return gain;
}

bool ncm_dmk_camera::is_auto_gain()
{
	return dmk_grabber_.isPropertyAutomationEnabled( VideoProcAmp_Gain );
}

bool ncm_dmk_camera::is_auto_exposure()
{
	return dmk_grabber_.isPropertyAutomationEnabled( CameraControl_Exposure );
}

bool ncm_dmk_camera::set_FPS(double new_FPS)
{
	dmk_grabber_.stopLive();

  if (dmk_grabber_.setFPS( new_FPS ))
    vcl_cout << "Frame rate changed to " << new_FPS << vcl_endl;
  else
    vcl_cout << "Frame rate not changed" << vcl_endl;

	dmk_grabber_.startLive( false );
	return true;
}
bool ncm_dmk_camera::set_exposure(double new_exposure)
{
  // Retrieve the absolute value interface for exposure.
	tIVCDAbsoluteValuePropertyPtr pAbsVal = 0;    
	if( pVCD_props_->findInterfacePtr( VCDID_Exposure, VCDElement_Value, pAbsVal ) != 0 )
	{
		pAbsVal->setValue( new_exposure );
		get_current_exposure( new_exposure );
		vcl_cout << "Exposure set to " << new_exposure << vcl_endl;
		return true;
	}
	else
	{
		vcl_cout << "DMK: Couldn't find exposure interface" << vcl_endl;
		return false;
	}
}
bool ncm_dmk_camera::set_gain(double new_gain)
{
	if ( dmk_grabber_.setProperty( VideoProcAmp_Gain, long( new_gain ) ) )
	{
		//vcl_cout << "Gain set to " << dmk_grabber_.getProperty(  VideoProcAmp_Gain ) << vcl_endl;
		return true;
	}
	else
	{
		vcl_cout << dmk_grabber_.getLastError().c_str() << vcl_endl;
		return false;
	}
}

bool ncm_dmk_camera::switch_auto_gain(bool auto_on)
{
	return dmk_grabber_.setProperty( VideoProcAmp_Gain, auto_on );
}

bool ncm_dmk_camera::switch_auto_exposure(bool auto_on)
{
	return dmk_grabber_.setProperty( CameraControl_Exposure, auto_on );
}

void ncm_dmk_camera::set_queue(ncm_video_frame_queue *_queue)
{
	queue_ = _queue;

	//Link the camera's queue to the filter
	frame_filter_.set_queue( queue_ );
}
