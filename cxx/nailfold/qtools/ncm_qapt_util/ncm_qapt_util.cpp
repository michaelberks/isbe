#include "ncm_qapt_util_gui.h"

#include <vcl_iostream.h>

#include <QtGui\QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	ncm_qapt_util_gui main_window;

#if _DEBUG
  //main_window.setFixedSize(1024, 720); 
  main_window.show();
#else
  //main_window.showMaximized();
  main_window.show();
#endif 

  // Get site name and such
  if (!main_window.initialize())
    return 1;

  int return_value = a.exec();

  return return_value;
}

#if _DEBUG
#else
// Include these lines to compile the application as a Windows app rather
// than a Console app, thus hiding the console window at startup. This is
// useful for the Release build. You'll also need to set the /SUBSYSTEM linker 
// option to WINDOWS (Project->Properties->Linker->System->SubSystem)

#include <Windows.h>
int WINAPI WinMain(HINSTANCE hInstance, 
                   HINSTANCE hPrevInstance, 
                   LPSTR lpCmdLine, 
                   int nShowCmd)
{
  return main(__argc, __argv);
}
#endif