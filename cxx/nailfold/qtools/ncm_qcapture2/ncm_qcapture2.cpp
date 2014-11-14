#include "ncm_qcapture2_gui.h"

#include <QtGui\QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	if( !DShowLib::InitLibrary() )
	{
		vcl_cout << "The IC Imaging Control Class Library could not be initialized.\n(invalid license key?)" << vcl_endl;
		exit( 1 );
	}
	ncm_qcapture2_gui main_window;
  main_window.show();

	return a.exec();
}
