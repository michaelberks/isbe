#include "ncm_qstitch_main.h"

#include <QtGui\QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	ncm_qstitch_main main_window;
  main_window.show();

	return a.exec();
}
