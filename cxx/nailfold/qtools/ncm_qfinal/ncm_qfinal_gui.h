#ifndef NAILFOLD_QFINAL_GUI_H
#define NAILFOLD_QFINAL_GUI_H

#include <QMainWindow>
#include <QWheelEvent>
#include <QProgressBar>
#include <QThread>
#include <QDir>
#include <QFileSystemModel>
#include <QSortFilterProxyModel>

#include "ui_ncm_qfinal_gui.h"
#include "ui_ncm_qfinal_options.h"

#include <vil/vil_image_view.h>

#include <nailfold/ncm_annotation.h>
#include <nailfold/ncm_image_server.h>

#include <nailfold/qtools/ncm_qscene.h>
#include <nailfold/qtools/ncm_qimagehandler.h>

#include <nailfold/qtools/ncm_qfinal/ncm_qfinal_preferences.h>
#include <nailfold/ncm_registry_accessor.h>
//#include <nailfold/qtools/dmk_listener.h>

class MySortFilterProxyModel : public QSortFilterProxyModel
{
protected:
    virtual bool MySortFilterProxyModel::filterAcceptsRow(
            int source_row, const QModelIndex &source_parent) const{
        
				QFileSystemModel *sm = qobject_cast<QFileSystemModel*>(sourceModel());
				QModelIndex idx = sourceModel()->index(source_row, 0 , source_parent);
				if (!sm->isDir(idx) )
					return QSortFilterProxyModel::filterAcceptsRow(source_row, source_parent);
        
				//if (source_parent == sm->index(sm->rootPath())) {    
        //    
        //} 
        return true;            
    }
};

class nailfold_qfinal_gui : public QMainWindow
{
  Q_OBJECT

public:
  nailfold_qfinal_gui(ncm_qfinal_preferences& preferences,
                       QWidget *parent = 0, Qt::WFlags flags = 0);
  ~nailfold_qfinal_gui();

  //: Non UI members manipulated by the UI
  void getFrame();

  void setDefaultPath(vcl_string image_path);

  bool initialize();
  bool wasConstructed() const;

public slots:

protected:
  void closeEvent(QCloseEvent* ev);

private slots:

	//Selecting cases
	void on_caseTreeView_clicked( const QModelIndex & index );
	void on_caseTreeView_doubleClicked( const QModelIndex & index );
	void on_actionSelectCase_triggered();

  // Menu/Toolbar actions
  void updateToolbar();	
  void on_actionLoadImage_triggered();
	void on_actionLoadMarkup_triggered();
  void on_actionSaveMarkup_triggered();
  void on_actionPreferences_triggered();
  void on_actionHelp_triggered();

  void on_actionPan_triggered();
  void on_actionLabelImage_triggered();
  void on_actionAddVessels_triggered();
  void on_actionDefineVesselSize_triggered();
  void on_actionDefineVesselShape_triggered();
  void on_actionLabelApices_triggered();
  void on_actionDrawVesselPath_triggered();
  void on_actionAddHaemorrhages_triggered();

  void on_actionDebug_triggered();
  void on_actionCrash_triggered();
  void on_actionWhatsThis_triggered();

  // Brightness and contrast controls
  void updateContrastControls();
  void on_brightnessSlider_sliderMoved(int value);
  void on_contrastSlider_sliderMoved(int value);
  void on_compSizeSlider_sliderMoved(int value);
  void on_compSizeSlider_sliderReleased();

  // Zoom controls
  void updateZoomControls();
  void updateZoomLimits();
  void on_zoomSlider_sliderMoved(int value);
  void on_zoomEdit_editingFinished();
  void on_zoomOutButton_clicked();
  void on_zoomInButton_clicked();
  void applyUserZoom();
  void applyVesselZoom(int vessel_index);

  // Help/Thumbnails panel
  void updateHelpPanel();
  void on_hidePanelButton_clicked();

	//Choosing what's visible in display mode
	void on_show_auto_haemorrhages_checkBox_toggled();
	void on_show_auto_paths_checkBox_toggled();
	void on_show_auto_placeholders_checkBox_toggled();
	void on_show_auto_apices_checkBox_toggled();
	void on_show_vessel_haemorrhages_checkBox_toggled();
	void on_show_vessel_paths_checkBox_toggled();
	void on_show_vessel_placeholders_checkBox_toggled();
	void on_show_vessel_apices_checkBox_toggled();
 
  // Choosing a vessel
  void updateVesselList();
  void clearVesselListSelection();
  void updateVesselControls();
  void on_firstVesselButton_clicked();
  void on_prevVesselButton_clicked();
  void on_nextVesselButton_clicked();
  void on_lastVesselButton_clicked();
  void on_deleteAllVessels_clicked();

  //: Set the grade of the disease
  void setGrade(int grade_index);

  void updateRowButtons(int selected);
  void on_prevRow_clicked();
  void on_nextRow_clicked();

  void updateSizeButtons(int selected);
  void on_prevSize_clicked();
  void on_nextSize_clicked();

  void updateShapeButtons(int selected);
  void on_prevShape_clicked();
  void on_nextShape_clicked();

  void updateStageControls();
  void on_prevStageButton_clicked();
  void on_nextStageButton_clicked();

  void updateImageControls();
  void on_prevImage_clicked();
  void on_nextImage_clicked();

  void on_applySizeToAll_clicked();
  void on_applyShapeToAll_clicked();

  void on_deleteAllHaemorrhages_clicked();

  // imageProcessor_ can't be added as a child of the main window because it
  // needs to be moved to its own thread later. Therefore, it can't be parsed
  // by the autoconnect mechanism.
  //void on_imageProcessor_imageLoaded(bool success);
  void onImageProcessor_imageLoaded(bool success);
  void onApexLengthChanged(double apex_width_pixels);

  void onThumbnail_clicked(vcl_string);

  void on_scene_editModeChanged(int edit_mode);
  void on_scene_mouseMoved(double, double);
  void on_graphicsView_viewModeChanged(int);

  void on_unmarkedCheckbox_toggled(bool checked);

  void updateStatus(const QString& status, double progress = 0.0);
  void updateStatusLeft();
  

private:
  //
  // Private member functions
  //
	bool loadImageFrom(vcl_string filename);
	bool loadMarkupFrom(vcl_string filename);

	void updateRootPath(vcl_string dir_name);
	void updateCasePath(vcl_string dir_name);
	void updateImagePath(vcl_string dir_name);
	void updateMarkupPath(vcl_string dir_name);

  void updateProgressFrame();
  void applyAnnotationChange();
  void updateItemVisibility();

  QString windowTitle() const;
  void hideScene();
  void showScene();
  bool getSiteName();
  void emptyOutbox();
  void pullImagesFromServer();
  void pullMarkupFromServer();
  bool closeIfFinished();

  void createModeActionGroup();
  void createGradeRadioGroup();
  void createRowRadioGroup();
  void createSizeRadioGroup();
  void createShapeRadioGroup();

  void applyUserLevelChange();
  void requestEditModeChange(int newMode);
  bool warnOnEditModeChange();

  //: Set icons to standard ones provided by Trolltech
  void setIcons();

  //: Set up the image processing thread
  void setupImageThread();

  //: Connect signals and slots between components
  void connectSignalsToSlots();
  void connectSetFocus();

  //: Set up the status bar how we want it
  void setupStatusbar();

  //: Set up the frame indicating progress in the labelling process
  void setupProgressFrame();

  //: Set the graphicsView zoom factor from a percentage value
  void setZoomFrom(int value);

  //: Add thumbnail images from same directory to display
  void getThumbnails();

  //: Delete the temporary backup file
  //  Returns true if successful
  bool deleteBackup() const;

  //: Load image from a file
	bool selectCaseDialog();
  bool loadImageDialog();
	bool loadMarkupDialog();

  //: Write markup data to a text file
  int saveMarkup();
  bool saveMarkupDialog();
  bool saveIfModified();

  //: Apply zoom as specified in user preferences
  void applyImageZoom(ncm_qfinal_preferences::ImageZoom zoom);

  void updateGridSpacing();

  void updateGradeRadios();
  void uncheckGradeRadios();

  //:Read directory paths from the registry
	void initialize_registry_accessor();
  void get_registry_paths();

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

  //: Handler for grabbing images from remote server
  ncm_server_handler* serverHandler_;

  //: File server that returns image filenames from local folder and
  ncm_file_server fileServer_;

  //: Thread for image processing operations
  QThread imageThread_;

  //: Class to handle image I/O and processing
  ncm_qimagehandler imageProcessor_;

  //: Alias to application-level preferences
  ncm_qfinal_preferences& preferences_;

  //: annotation class that stores image properties
  ncm_annotation markup_;

  //: canvas on which to draw stuff
  //  Note that this appears after markup_ in the declaration since scene_ needs
  //  the data in markup_ to destruct itself and must therefore be destroyed
  //  first (i.e. declared later than markup_)
  ncm_qscene scene_;

  //: Current image folder
	vcl_string rootCasePath_;
	vcl_string currentCasePath_;
  vcl_string currentImagePath_;
	vcl_string currentMarkupPath_;

  //: Directory listing of currentImagePath_
  QDir currentImageDir_;
	QDir currentMarkupDir_;

  //: Registry accessor
  ncm_registry_accessor ra_;

  //: Labels to go in status bar. The status bar will contain one main label,
  //  plus two panels at the right hand side to indicate the function of the 
  //  left and right mouse buttons at any given time
  QLabel statusBar_main_;
  QProgressBar statusBar_progress_;
  QLabel statusBar_LPanel_;
  QLabel statusBar_RPanel_;


  //: Widgets added manually
  QActionGroup* modeGroup_;
  QButtonGroup* gradeRadioGroup_;
  QButtonGroup* rowRadioGroup_;
  QButtonGroup* sizeRadioGroup_;
  QButtonGroup* shapeRadioGroup_;

  vcl_vector<QLabel*> progressLabels_;

	//:File server model to browse for cases
	QFileSystemModel* baseFileModel_;
	MySortFilterProxyModel* imageFilterModel_;
	MySortFilterProxyModel* markupFilterModel_;
};

#endif // NAILFOLD_QFINAL_GUI_H
