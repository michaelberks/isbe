#include "ncm_qtool_tmp.h"

#include <QtGui\QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	ncm_qtool_tmp main_window;
	main_window.show();
	return a.exec();
}
