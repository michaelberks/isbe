#include "ncm_qmarkup_gui.h"

#include <QtGui\QApplication>

#include <nailfold/qtools/ncm_qmarkup/ncm_qmarkup_preferences.h>

int main(int argc, char *argv[])
{
  ncm_qmarkup_preferences preferences;

	QApplication a(argc, argv);
  if (a.arguments().count() > 1)
  {
    // Application called with extra parameter(s)
    QString arg1 = a.arguments().at(1);

    // If called with -local argument then don't use the remote server
    if (arg1.compare("-local") == 0)
    {
      preferences.set_use_remote_server(false);
      vcl_cout << "Local access only" << vcl_endl;
    }
  }

	nailfold_qmarkup_gui main_window(preferences);

#if _DEBUG
  main_window.setFixedSize(1024, 720); 
  main_window.show();
#else
  main_window.showMaximized();
#endif 

  // Get site name and such
  if (!main_window.initialize())
    return 1;

  if (argc>1)
    main_window.setDefaultPath(argv[1]);

  int return_value = a.exec();
  // preferences.write_to_registry();

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