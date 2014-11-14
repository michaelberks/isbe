#ifndef NAILFOLD_QMARKUP_GUI_H
#define NAILFOLD_QMARKUP_GUI_H

#include <QMainWindow>
#include <QWheelEvent>
#include <QProgressBar>
#include <QThread>
#include <QDir>

#include "ui_ncm_qseries_gui.h"
#include "ui_ncm_qseries_options.h"

#include <vil/vil_image_view.h>

#include <nailfold/qtools/ncm_qseries/ncm_qseries_grade.h>

#include <nailfold/qtools/ncm_qscene.h>
#include <nailfold/qtools/ncm_qimagehandler.h>

#include <nailfold/qtools/ncm_qseries/ncm_qseries_image_manager.h>
#include <nailfold/qtools/ncm_qseries/ncm_qseries_preferences.h>
//#include <nailfold/qtools/dmk_listener.h>

class nailfold_qseries_gui : public QMainWindow
{
  Q_OBJECT

public:
	enum GradeStage { Stage1 = 0,
                    Stage2 = 1};

  nailfold_qseries_gui(ncm_qseries_preferences& preferences,
                       QWidget *parent = 0, Qt::WFlags flags = 0);
  ~nailfold_qseries_gui();

  //: Non UI members manipulated by the UI
  void getFrame();

  bool initialize();
  bool wasConstructed() const;

public slots:

protected:
  void closeEvent(QCloseEvent* ev);

private slots:

  bool loadImage();

  // Menu/Toolbar actions
  void updateToolbar();
  void on_actionSaveMarkup_triggered();
  void on_actionPreferences_triggered();

  // Brightness and contrast controls
  void updateContrastControls();
  void on_brightnessSlider_sliderMoved(int value);
  void on_contrastSlider_sliderMoved(int value);
	void on_brightnessSlider_2_sliderMoved(int value);
  void on_contrastSlider_2_sliderMoved(int value);

	//Choosing grade and selecting reasons
	void setGrade(int normal);
	void selectReasons(int reason_index);

	//Moving the progression slider
	void on_progressionSlider_sliderMoved(int value);

  // Zoom controls
  void updateZoomControls();
  void updateZoomLimits();
  void on_zoomSlider_sliderMoved(int value);
  void on_zoomEdit_editingFinished();
  void on_zoomOutButton_clicked();
  void on_zoomInButton_clicked();
	void on_gridButton_toggled();
  void applyUserZoom();

  // Help/Thumbnails panel
  void updateHelpPanel();
  void on_hidePanelButton_clicked();

  void updateStageControls();
	void on_goBackButton_clicked();
	void on_continueButton_clicked();
  
  void requestPreviousImage();
  void requestNextImage();

  // imageProcessor_ can't be added as a child of the main window because it
  // needs to be moved to its own thread later. Therefore, it can't be parsed
  // by the autoconnect mechanism.
  //void on_imageProcessor_imageLoaded(bool success);
  void onImageProcessor_imageLoaded(bool success);
	void onImageProcessor2_imageLoaded(bool success);

  void on_graphicsView_viewModeChanged(int);

  void on_unmarkedCheckbox_toggled(bool checked);

  void updateStatus(const QString& status, double progress = 0.0);
  void updateStatusLeft();
  

private:
  //
  // Private member functions
  //

  QString windowTitle() const;
  void hideScene();
  void showScene();
  bool getSiteName();
  bool closeIfFinished();

  void createGradeRadioGroup();
	void createReasonsRadioGroup();

  void applyUserLevelChange();
  void requestStageChange(int newStage);
  bool warnOnStageChange();

	void goToStage(int stage);

  //: Set icons to standard ones provided by Trolltech
  void setIcons();

  //: Set up the image processing thread
  void setupImageThread();

  //: Connect signals and slots between components
  void connectSignalsToSlots();
  void connectSetFocus();

  //: Set up the status bar how we want it
  void setupStatusbar();

  //: Set the graphicsView zoom factor from a percentage value
  void setZoomFrom(int value);

  //: Delete the temporary backup file
  //  Returns true if successful
  bool deleteBackup() const;

  //: Load image from a file
  bool loadImageDialog();

  //: Write markup data to a text file
  int saveMarkup();
  bool saveMarkupDialog();
  bool saveIfModified();

  //: Apply zoom as specified in user preferences
  void applyImageZoom(ncm_qseries_preferences::ImageZoom zoom);

  void updateGridSpacing();

  void updateGradeRadios();
  void uncheckGradeRadios();

	void updateReasonsRadios();
  void uncheckReasonsRadios();

  //
  // Private member variables
  //

  //: User interface object
  Ui::NailfoldMainWindow ui;

  //: User interface object
  Ui::optionsDialog optionsUi;

  //: Flag set to true if the window was constructed completely
  bool constructed_;

  QStyle* style_;

  //: File server that returns image filenames from local folder and
  ncm_qseries_image_manager imageManager_;

  //: Thread for image processing operations
  QThread imageThread_;

  //: Class to handle image I/O and processing
  ncm_qimagehandler imageProcessor_;
	ncm_qimagehandler imageProcessor2_;

  //: Alias to application-level preferences
  ncm_qseries_preferences& preferences_;

	//: Current stage of grading for this image
	GradeStage currentStage_;

  //: annotation class that stores image properties
  ncm_qseries_grade grade_;

  //: canvas on which to draw stuff
  //  Note that this appears after markup_ in the declaration since scene_ needs
  //  the data in markup_ to destruct itself and must therefore be destroyed
  //  first (i.e. declared later than markup_)
  ncm_qscene scene_;
  ncm_qscene scene2_;

  //: Labels to go in status bar. The status bar will contain one main label,
  //  plus two panels at the right hand side to indicate the function of the 
  //  left and right mouse buttons at any given time
  QLabel statusBar_main_;
  QProgressBar statusBar_progress_;
  QLabel statusBar_LPanel_;
  QLabel statusBar_RPanel_;


  //: Widgets added manually
  QButtonGroup* gradeRadioGroup_;
  QButtonGroup* reasonsRadioGroup_;

  vcl_vector<QLabel*> progressLabels_;
};

#endif // NAILFOLD_QMARKUP_GUI_H
