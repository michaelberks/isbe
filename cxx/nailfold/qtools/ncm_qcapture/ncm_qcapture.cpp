#include "ncm_qcapture_gui.h"
#include "ncm_qcapture_preferences.h"

#include <QtGui\QApplication>
//#include <initguid.h>
//#include <dbt.h>

//DEFINE_GUID(GUID_DEVINTERFACE_USB_DEVICE, 0xA5DCBF10L, 0x6530, 0x11D2, 0x90, 0x1F, 0x00, 0xC0, 0x4F, 0xB9, 0x51, 0xED);

int main(int argc, char *argv[])
{

	QApplication a(argc, argv);
	qRegisterMetaType<ncm_video_frame_header>();
	if( !DShowLib::InitLibrary() )
	{
		vcl_cout << "The IC Imaging Control Class Library could not be initialized.\n(invalid license key?)" << vcl_endl;
		exit( 1 );
	}

	
  /*DEV_BROADCAST_DEVICEINTERFACE NotificationFilter;
     
  ZeroMemory( &NotificationFilter, sizeof(NotificationFilter) );
  NotificationFilter.dbcc_size = sizeof(DEV_BROADCAST_DEVICEINTERFACE);
  NotificationFilter.dbcc_devicetype = DBT_DEVTYP_DEVICEINTERFACE;
  NotificationFilter.dbcc_classguid  = GUID_DEVINTERFACE_USB_DEVICE;*/
     
  

	ncm_qcapture_preferences preferences;
	
  /*if (a.arguments().count() > 1)
  {
    // Application called with extra parameter(s)
    QString arg1 = a.arguments().at(1);

    // If called with -local argument then don't use the remote server
    if (arg1.compare("-local") == 0)
    {
			preferences.use_sql_database_ = false;
      vcl_cout << "No database access" << vcl_endl;
    }
  }*/

	ncm_qcapture_gui main_window(preferences);
	a.installEventFilter(&main_window);
	//HDEVNOTIFY hDevNotify = RegisterDeviceNotification(main_window.winId(), &NotificationFilter, DEVICE_NOTIFY_ALL_INTERFACE_CLASSES);
	main_window.show();

	return a.exec();
}
