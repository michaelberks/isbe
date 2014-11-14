#include "ncm_qfinal_gui.h"

#include <QFileDialog>
#include <QBitmap>
#include <QGraphicsPixmapItem>
#include <QMessageBox>
#include <QActionGroup>
#include <QButtonGroup>
#include <QWindowsXPStyle>
#include <QWhatsThis>
#include <QFileInfo>
#include <QDateTime>
#include <QDesktopServices>
#include <QRegExp>
#include <QLibraryInfo>
#include <Qt/QSqlDatabase.h>
#include <Qt/QSqlQuery.h>

//#include <vcl_cmath.h>
#include <vcl_iostream.h>

//#include <vnl/vnl_random.h>

#include <vul/vul_file.h>

#include <nailfold/ncm_vessel.h>
#include <nailfold/ncm_scp_handler.h>
#include <nailfold/ncm_sftp_handler.h>

#include <nailfold/qtools/ncm_qthumbnail.h>
#include <nailfold/qtools/ncm_qvesselitem.h>
#include <nailfold/qtools/ncm_qvesselitem_appearance.h>
#include <nailfold/qtools/ncm_qapexitem.h>
#include <nailfold/qtools/ncm_qapexitem_appearance.h>
#include <nailfold/qtools/ncm_qhaemorrhageitem.h>
#include <nailfold/qtools/ncm_qhaemorrhageitem_appearance.h>
#include <nailfold/qtools/ncm_qsite_selector.h>

#include "ncm_qfinal_options.h"

//
//:
void nailfold_qfinal_gui::on_actionDebug_triggered()
{
  int dummy = 0;
}
//
//: Cause the application to crash - allows us to test what happens in this
//  eventuality (e.g. test that autosave works)
void nailfold_qfinal_gui::on_actionCrash_triggered()
{
  assert(false);
}

//
//: Constructor
nailfold_qfinal_gui::nailfold_qfinal_gui(ncm_qfinal_preferences& preferences,
                                           QWidget *parent, Qt::WFlags flags)
: QMainWindow(parent, flags),
	rootCasePath_(""),
	currentCasePath_(""),
  currentImagePath_(""),
	currentMarkupPath_(""),
  preferences_(preferences),
  constructed_(false),
  scene_(this)
{
  // Name objects to allow autoconnection of signals to slots
  scene_.setObjectName("scene");
  imageProcessor_.setObjectName("imageProcessor");

  // setup the UI
  ui.setupUi(this);

  setWindowTitle(windowTitle());

  style_ = new QWindowsXPStyle;

  // Get standard icons from the platform-dependent theme
  setIcons();

  // associate the scene with the markup and image processor
  scene_.set_annotation(&markup_);
  scene_.set_image_processor(&imageProcessor_);
  scene_.set_grid_visible(ui.actionToggleGrid->isChecked());

  const double border = preferences_.normal_vessel_zoom() *
                        preferences_.enlarged_zoom_relative() * 
                        preferences_.giant_zoom_relative();
  scene_.set_border_size(border);

  // use Pan mode by default, with reverse zoom (wheel up = zoom in)
  ui.graphicsView->useMarkupMode();
  ui.graphicsView->reverse_zoom();
  ui.graphicsView->setFocus();

  // tell graphics view to use this scene
  ui.graphicsView->setScene(&scene_);
  ui.zoomedView->setScene(&scene_);
  ui.zoomedView->scale(2,2);

  // Make first item selected on startup
  ui.stackedWidget->setCurrentIndex(0);

	//Instantiate file model
	baseFileModel_ = new QFileSystemModel();
	
	//Link case browser tree view with file model
	imageFilterModel_ = new MySortFilterProxyModel();
	imageFilterModel_->setSourceModel( baseFileModel_ );
	ui.caseTreeView->setModel(imageFilterModel_);
	QRegExp rx(".*.png$|.*.bmp$|.*.jpg$", Qt::CaseInsensitive);
	imageFilterModel_->setFilterRegExp(rx);

	//Link markup browser list view with file model
	markupFilterModel_ = new MySortFilterProxyModel();
	markupFilterModel_->setSourceModel( baseFileModel_ );
	ui.markupListView->setModel(markupFilterModel_);
  
  // Create and initialize various objects
  createModeActionGroup();
  createGradeRadioGroup();
  createRowRadioGroup();
  createSizeRadioGroup();
  createShapeRadioGroup();

  setupProgressFrame();
  setupStatusbar();
  setupImageThread();

  updateGridSpacing();

  // Set 'normal' apexitem line width to 2 microns
  const double pixels_per_micron = preferences_.grid_pixels_per_mm() / 1000.0;
  QGraphicsApexItem::appearance().set_line_width(2.0 * pixels_per_micron);
  QGraphicsApexItem::appearance().set_width_relative_to_image_size(false);

  // Set visible UI elements to reflect the level of the user
  applyUserLevelChange();

  connectSignalsToSlots();

	//Just tryopening a connecting to the nailfold database here... we'll play around later
	//with exactly what and where this should be done
	vcl_cout << "Path to plugins " << QLibraryInfo::location (QLibraryInfo::PluginsPath).toStdString() << vcl_endl;
	QSqlDatabase db = QSqlDatabase::addDatabase("QMYSQL", "nailfold_connection");
  db.setHostName("");
  db.setDatabaseName("nailfold_initial_test");
  db.setUserName("root");
  db.setPassword("Gooners12!");
  bool ok = db.open();

	if (ok)
	{
		vcl_cout << "Yay, we connected to the database!" << vcl_endl;
		QSqlQuery query(db);
		query.exec("SELECT username FROM user");
     while (query.next()) {
         QString user = query.value(0).toString();
				 vcl_cout << "User " << user.toStdString() << vcl_endl;
     }

		 query.exec("SELECT study_name FROM study");
     while (query.next()) {
         QString study = query.value(0).toString();
				 vcl_cout << "Study " << study.toStdString() << vcl_endl;
     }

		 query.exec("SELECT subject_study_name FROM subject");
     while (query.next()) {
         QString subject = query.value(0).toString();
				 vcl_cout << "Subject " << subject.toStdString() << vcl_endl;
     }
	}

#if _DEBUG
#else
  // hide the Debug button in the menubar for the Release version
  ui.actionDebug->setVisible(false);
  ui.actionCrash->setVisible(false);
#endif 
  
  // Flag to show that constructor ended successfully
  constructed_ = true;
}

//
//: Destructor
nailfold_qfinal_gui::~nailfold_qfinal_gui()
{
  // Empty the outbox again before you go (just in case).
  emptyOutbox();

  imageThread_.quit();

  // Stop taking snapshots now that we're finished.
  // Otherwise, when the scene cleans itself up by deleting the ncm_vesselitems,
  // the subsequent deletion of ncm_vessels causes a new snapshot to be taken
  // such that the backup file is created again.
  markup_.set_snapshots(false);

  // If backup file exists then delete it
  vcl_string backup_filename = markup_.filename() + '~';
  QFile(backup_filename.c_str()).remove();
}

bool nailfold_qfinal_gui::initialize()
{
	initialize_registry_accessor();

  if (preferences_.use_remote_server())
  {
    // Use faculty's SFTP server
    serverHandler_ = new ncm_sftp_handler(/* server = */ "³æôð®ò³³®­è³®­¡î®¡£®µ«",
                                          /* username = */ "ðôò¡î",
                                          /* password = */ "î¡©ìæ¯ìä²°°°");

    // Use ISBE's SCP server (nimrod)
    //serverHandler_ = new ncm_scp_handler(/* server = */ "qs°®¸¸®²´¶®qu",
    //                                     /* username = */ "ðôò¡î",
    //                                     /* password = */ "ÖuyðÌ¡");

    const bool connected = serverHandler_->connect();
    if (connected)
    {
      // Get the name of the site (stored in <root>/site.txt), send any unsent
      // markup, and retrieve any new images for marking up.
      getSiteName();
      emptyOutbox();

      statusBar_main_.setText("Getting data from server...");

      QMessageBox msgBox;
      msgBox.setText("Click OK to retrieve images and markup from remote server.");
      msgBox.setInformativeText("(This may take a few minutes.)");
      msgBox.setStandardButtons(QMessageBox::Ok);
      msgBox.setDefaultButton(QMessageBox::Ok);
      msgBox.setModal(false);
      msgBox.exec();

      pullMarkupFromServer();
      pullImagesFromServer();

      msgBox.close();
      statusBar_main_.setText("");

      // If there are no images left to mark and the user agrees to quit, return
      // immediately so that the application can close down.
      if (closeIfFinished())
        return false;
    }
    else
    {
      // Issue a warning if connect failed
      QMessageBox msgBox;
      msgBox.setIcon(QMessageBox::Warning);
      msgBox.setText("Unable to connect to remote server.");
      msgBox.exec();
    }
  }
	else
	{
		QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.setText("Not using remote server.");
    msgBox.exec();

		//Get paths to images on local machine
		get_registry_paths();

		//All the images should be on the local machine, just need to initialise the inbox
		//fileServer_.parse_images();
	}

	//By deafult we want to show images regardless of markup
	ui.unmarkedCheckbox->setChecked(false);
	fileServer_.use_unmarked_only(false);

	/*
  // Go to first image in Inbox and load it
  bool first_found = fileServer_.first();
  if (first_found) {
    loadImageFrom(fileServer_.current_image());
	}
  else
    hideScene();
	*/
  return true;
}

QString nailfold_qfinal_gui::windowTitle() const
{
  QString version_string = QString::number(preferences_.version());
  return QString("NCM Markup v" + version_string);
}

bool nailfold_qfinal_gui::wasConstructed() const
{
  return constructed_;
}

void nailfold_qfinal_gui::hideScene()
{
  ui.graphicsView->setScene(NULL);

  ui.imageGroup->setEnabled(false);
  ui.zoomGroup->setEnabled(false);

  ui.vesselsGroup->setEnabled(false);

  ui.stackedWidget->setEnabled(false);
  ui.prevStageButton->setEnabled(false);
  ui.nextStageButton->setEnabled(false);
}

void nailfold_qfinal_gui::showScene()
{
  ui.graphicsView->setScene(&scene_);
  ui.graphicsView->fit_both();

  ui.imageGroup->setEnabled(true);
  ui.zoomGroup->setEnabled(true);

  updateVesselControls();
  updateStageControls();
}

bool nailfold_qfinal_gui::getSiteName()
{
  int error_code = -1;

  error_code = fileServer_.read_site_name();

  if (error_code == 0)
    return true;

  if (error_code == 1)
  {
    // File not found - create a new one using the site selector
    ncm_qsite_selector siteSelector(this);
    siteSelector.exec();

    // Try once more - if it worked, great
    error_code = fileServer_.read_site_name();
    if (error_code == 0)
      return true;
  }

  // Otherwise report an error to the user
  QMessageBox msgBox;
  msgBox.setIcon(QMessageBox::Warning);
  switch (error_code)
  {
  case 1:
    msgBox.setText("site.txt not found");
    msgBox.exec();
    return false;

  case 4:
    msgBox.setText("Unknown user found in site.txt");
    msgBox.exec();
    return true; // we can still run with this

  default:
    assert(false);
  }

  return true; // Shouldn't get this far
}

void nailfold_qfinal_gui::emptyOutbox()
{
  if (!preferences_.use_remote_server())
    return;

  // Empty the Outbox to the remote server if possible
  statusBar_main_.setText("Sending outbox to server...");
  fileServer_.push_outbox_to(*serverHandler_);
  statusBar_main_.setText("");
}

void nailfold_qfinal_gui::pullImagesFromServer()
{
  if (!preferences_.use_remote_server())
    return;

  int list_error = fileServer_.pull_imagelist_from(*serverHandler_);
  int pull_error = fileServer_.pull_images_from(*serverHandler_);
  fileServer_.parse_inbox();
  
  if (list_error == -1)
  {
    QMessageBox msgBox;
    msgBox.setText("Unable to get list of images from remote server.");
    msgBox.exec();
  }
  else if (pull_error > 0)
  {
    QMessageBox msgBox;
    msgBox.setText(QString::number(pull_error) + 
                   " images failed to copy from remote server.");
    msgBox.exec();
  }
}

void nailfold_qfinal_gui::pullMarkupFromServer()
{
  if (!preferences_.use_remote_server())
    return;

  statusBar_main_.setText("Getting existing markup from server...");
  int pull_error = fileServer_.pull_markup_from(*serverHandler_);
  statusBar_main_.setText("");

  if (pull_error != 0)
  {
    QMessageBox msgBox;
    msgBox.setText("Unable to get existing annotations from remote server.");
    msgBox.exec();
  }
}

//
//: Optionally close the application if there are no more images to mark up
bool nailfold_qfinal_gui::closeIfFinished()
{
  // If that was the last image then ask user if they want to close
  if (fileServer_.is_empty())
  {
    QMessageBox msgBox;
    msgBox.setIcon(QMessageBox::Warning);
    msgBox.setText("No more images to mark up");
    msgBox.setInformativeText("Close QMarkup?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::Yes);

    const int response = msgBox.exec();

    if (response == QMessageBox::Yes)
    {
      close();
      return true;
    }
  }

  return false;
}

void nailfold_qfinal_gui::setIcons()
{
  // fromTheme doesn't work automatically with Windows
  //ui.actionLoadImage->setIcon(QIcon::fromTheme("document-open"));
  //ui.actionSaveMarkup->setIcon(QIcon::fromTheme("document-save"));

  if (style_ != NULL)
  {
    ui.actionLoadImage->setIcon(style_->standardIcon(QStyle::SP_DialogOpenButton));
    ui.actionSaveMarkup->setIcon(style_->standardIcon(QStyle::SP_DialogSaveButton));
    ui.actionPreferences->setIcon(style_->standardIcon(QStyle::SP_FileDialogContentsView));
    ui.actionWhatsThis->setIcon(style_->standardIcon(QStyle::SP_MessageBoxQuestion));
  }
}

void nailfold_qfinal_gui::setupImageThread()
{
  // Assign all image processing to a separate thread to avoid the GUI hanging
  // up on large images
  imageProcessor_.moveToThread(&imageThread_);

  // Have the thread running at all times to handle brightness and contrast
  // changes
  imageThread_.start();
}

//
//: Create the status bar's labels and progress bar
void nailfold_qfinal_gui::setupStatusbar()
{
  statusBar_main_.setFrameStyle(QFrame::Panel & QFrame::Sunken);
  statusBar_main_.setLineWidth(1);
  statusBar_main_.setText("");
  statusBar_main_.setContentsMargins(4, 0, 4, 0);
  statusBar()->addWidget(&statusBar_main_, 1);

  statusBar_progress_.setFixedWidth(160);
  statusBar()->addWidget(&statusBar_progress_);

  statusBar_LPanel_.setFrameStyle(QFrame::Panel & QFrame::Sunken);
  statusBar_LPanel_.setLineWidth(1);
  statusBar_LPanel_.setFixedWidth(120);
  statusBar_LPanel_.setText("");
  statusBar()->addWidget(&statusBar_LPanel_);

  statusBar_RPanel_.setFrameStyle(QFrame::Panel & QFrame::Sunken);
  statusBar_RPanel_.setLineWidth(1);
  statusBar_RPanel_.setFixedWidth(120);
  statusBar_RPanel_.setText("");
  statusBar()->addWidget(&statusBar_RPanel_);
}

//
//: Update status bar
void nailfold_qfinal_gui::updateStatus(const QString& status, 
                                         double progress)
{
  statusBar_main_.setText(status);

  if (progress >= 0.0)
  {
    const int progress_min = statusBar_progress_.minimum();
    const int progress_range = statusBar_progress_.maximum() - progress_min;
    statusBar_progress_.setValue(progress_min + progress*progress_range);
  }
  else
    statusBar_progress_.reset();

  statusBar()->update();
}

void nailfold_qfinal_gui::updateStatusLeft()
{
  QString statusString = QString::number(fileServer_.n_unmarked()) + 
                         " images left to mark";
  statusBar_LPanel_.setText(statusString);
}

//
//: Update the status bar when an apex changes length
void nailfold_qfinal_gui::onApexLengthChanged(double apex_width_pixels) 
{
  if (apex_width_pixels >= 0.0)
  {
    const double apex_width_mm = apex_width_pixels / 
                                 preferences_.grid_pixels_per_mm();

    statusBar_main_.setText("Apex width = " + 
                            QString::number(1000.0 * apex_width_mm) +
                            " micron");
  }
  else
  {
    statusBar_main_.setText("");
  }
}


//
//: Set up frame that indicates progress in the labelling process
void nailfold_qfinal_gui::setupProgressFrame()
{
  /*QHBoxLayout* pbl = new QHBoxLayout();
  pbl->setSpacing(0);
  QHBoxLayout* pll = new QHBoxLayout();
  pll->setSpacing(19);
  pll->setContentsMargins(19, 0, 19, 0);

  QSizePolicy joinPolicy(QSizePolicy::Fixed, QSizePolicy::Minimum);
  QSizePolicy centrePolicy(QSizePolicy::Expanding, QSizePolicy::Minimum);
  
  // Add body of 'start' block
  QLabel* leftEndLabel = new QLabel();
  leftEndLabel->setSizePolicy(joinPolicy);
  leftEndLabel->setPixmap(QPixmap(":/images/progress/left_end.png"));
  pbl->addWidget(leftEndLabel);

  QLabel* centreLabel = new QLabel();
  centreLabel->setPixmap(QPixmap(":/images/progress/centre.png"));
  centreLabel->setSizePolicy(centrePolicy); // expand the label
  centreLabel->setScaledContents(true); // and its contents, too
  pbl->addWidget(centreLabel);

  QLabel* textLabel = new QLabel("Start");
  textLabel->setSizePolicy(centrePolicy);
  textLabel->setAlignment(Qt::AlignCenter);
  textLabel->setEnabled(false);
  pll->addWidget(textLabel);

  // At least one join required
  QLabel* joinLabel = new QLabel();
  joinLabel->setSizePolicy(joinPolicy);
  joinLabel->setPixmap(QPixmap(":/images/progress/join.png"));
  pbl->addWidget(joinLabel);

  const int nControls = ui.stackedWidget->count();
  for (int i = 0; i < nControls; ++i)
  {
    centreLabel = new QLabel();
    centreLabel->setPixmap(QPixmap(":/images/progress/centre.png"));
    centreLabel->setSizePolicy(centrePolicy);
    centreLabel->setScaledContents(true);
    pbl->addWidget(centreLabel);

    // Match the label to the title of the stacked widget panel.
    QWidget* w = ui.stackedWidget->widget(i);
    QString label_string("Stage");
    for (int iChild = 0; iChild < w->children().size(); ++iChild)
    {
      QGroupBox* group = qobject_cast<QGroupBox*>(w->children()[iChild]);
      if (group != NULL)
        label_string = group->title();
    }
    textLabel = new QLabel(label_string);
    textLabel->setSizePolicy(centrePolicy);
    textLabel->setAlignment(Qt::AlignCenter);

    // Disable the label if the stackedWidget page is disabled.
    const bool stageIsEnabled = ui.stackedWidget->widget(i)->isEnabled();
    textLabel->setEnabled(stageIsEnabled);

    pll->addWidget(textLabel);
    progressLabels_.push_back(textLabel);

    joinLabel = new QLabel();
    joinLabel->setSizePolicy(joinPolicy);
    joinLabel->setPixmap(QPixmap(":/images/progress/join.png"));
    pbl->addWidget(joinLabel);
  }

  // Finally, add the body of the right hand label
  centreLabel = new QLabel("C");
  centreLabel->setPixmap(QPixmap(":/images/progress/centre.png"));
  centreLabel->setScaledContents(true);
  centreLabel->setSizePolicy(centrePolicy);
  pbl->addWidget(centreLabel);

  textLabel = new QLabel("Finish");
  textLabel->setSizePolicy(centrePolicy);
  textLabel->setAlignment(Qt::AlignCenter);
  textLabel->setEnabled(false);
  pll->addWidget(textLabel);

  QLabel* rightEndLabel = new QLabel("R");
  rightEndLabel->setSizePolicy(joinPolicy);
  rightEndLabel->setPixmap(QPixmap(":/images/progress/right_end.png"));
  pbl->addWidget(rightEndLabel);

  // Make gridLayout the chosen layout for the frame
  QGridLayout* gridLayout = new QGridLayout();
  gridLayout->setSpacing(0);
  gridLayout->setContentsMargins(0, 0, 0, 0);
  gridLayout->addLayout(pbl, 0, 0);
  gridLayout->addLayout(pll, 0, 0);
  ui.progressFrame->setLayout(gridLayout);*/
}

//
//: Update paths for the case, images and markup
void nailfold_qfinal_gui::updateRootPath(vcl_string dir_name)
{
	if (dir_name != rootCasePath_)
	{
		//Update the root path and also set the current case to this
		rootCasePath_ = dir_name;
		currentCasePath_ = dir_name;

		//Write the new value to the registry
		ra_.write_string("ncmRootCaseDir", dir_name);
	}
}

void nailfold_qfinal_gui::updateCasePath(vcl_string dir_name)
{
	if (dir_name != currentCasePath_)
	{
		currentCasePath_ = dir_name;

		baseFileModel_->setRootPath(QString(dir_name.c_str()));	
		ui.caseTreeView->setRootIndex(imageFilterModel_->mapFromSource(baseFileModel_->index(QString(dir_name.c_str()))));
		
	}
}

void nailfold_qfinal_gui::updateImagePath(vcl_string dir_name)
{
	if (dir_name != currentImagePath_)
	{
		currentImagePath_ = dir_name;
		fileServer_.set_image_dir(dir_name);
		fileServer_.parse_images();

		// store image path and basename of image
		currentImageDir_.setPath(currentImagePath_.c_str());
		QStringList filters;
		filters << "*.png" << "*.bmp" << "*.jpg";
		currentImageDir_.setNameFilters(filters);

		getThumbnails();
	}
}

void nailfold_qfinal_gui::updateMarkupPath(vcl_string dir_name)
{
	if (dir_name != currentMarkupPath_)
	{
		currentMarkupPath_ = dir_name;
		fileServer_.set_markup_dir(dir_name);

		// store image path and basename of image
		currentMarkupDir_.setPath(currentMarkupPath_.c_str());
		QStringList filters;
		filters << "*.txt";
		currentMarkupDir_.setNameFilters(filters);

		//Write the new value to the registry
		ra_.write_string("ncmMarkupDir", dir_name);
	}
}

void nailfold_qfinal_gui::updateProgressFrame()
{
  /*int edit_mode = scene_.edit_mode();

  // Update progress indicator at top of screen
  for (int i = 0; i < ui.stackedWidget->count(); ++i)
  {
    QWidget* w = ui.stackedWidget->widget(i);
    QString label_string("Oops!");
    for (int iChild = 0; iChild < w->children().size(); ++iChild)
    {
      QGroupBox* group = qobject_cast<QGroupBox*>(w->children()[iChild]);
      if (group != NULL)
        label_string = group->title();
    }
    if (i == edit_mode)
      progressLabels_[i]->setText("<b>"+label_string+"</b>");
    else
      progressLabels_[i]->setText(label_string);

    // Disable the label if the stackedWidget page is disabled.
    const bool stageIsEnabled = ui.stackedWidget->widget(i)->isEnabled();
    progressLabels_[i]->setEnabled(stageIsEnabled);
  }*/
}

//
//: Apply a specified image zooming to the graphics view
void nailfold_qfinal_gui::applyImageZoom(ncm_qfinal_preferences::ImageZoom zoom)
{
  switch (zoom)
  {
    case ncm_qfinal_preferences::ZoomNoResize:
      // do nothing
      break;
    case ncm_qfinal_preferences::ZoomFit:
      ui.graphicsView->fit_both();
      break;
    case ncm_qfinal_preferences::ZoomFitWidth:
      ui.graphicsView->fit_width();
      break;
    case ncm_qfinal_preferences::ZoomFitHeight:
      ui.graphicsView->fit_height();
      break;
    default:
      // throw error?
      break;
  }
}

//
//: Apply a specified image zooming to the graphics view
void nailfold_qfinal_gui::applyVesselZoom(int vessel_index)
{
  ncm_vessel* vessel = markup_.vessel(vessel_index);
  if (vessel == NULL)
    return;

  // Get normal zoom first
  double zoom = preferences_.normal_vessel_zoom(); // default
  ncm_vessel_properties vessel_props = vessel->properties();
  
  // Zoom further out for larger vessels
  if (vessel_props.is_size_enlarged())
    zoom *= preferences_.enlarged_zoom_relative();
  else if (vessel_props.is_size_giant() ||
           vessel_props.is_size_irregular())
  {
    zoom *= preferences_.enlarged_zoom_relative();
    zoom *= preferences_.giant_zoom_relative();
  }

  // Apply a shift if we're drawing the vessel path
  double y_shift = 0.0;
  if (scene_.edit_mode() == ncm_qscene::ModeDrawVesselPath)
    y_shift = preferences_.vessel_path_yshift();

  ui.graphicsView->focus_on_vessel(vessel_index, zoom, y_shift);
}

void nailfold_qfinal_gui::applyUserZoom()
{
  switch (scene_.edit_mode())
  {
    case ncm_qscene::ModeClassifyImage:
      applyImageZoom(preferences_.label_image_zoom());
      break;
    case ncm_qscene::ModeAddVessels:
      applyImageZoom(preferences_.add_vessels_zoom());
      break;
    case ncm_qscene::ModeSetVesselSize: // fall through
    case ncm_qscene::ModeSetVesselShape:
      applyImageZoom(preferences_.set_vessel_props_zoom());
      break;
    case ncm_qscene::ModeLabelApices: // fall through
    case ncm_qscene::ModeDrawVesselPath:
      applyVesselZoom(ui.vesselCombo->currentIndex());
      break;
    case ncm_qscene::ModeAddHaemorrhages:
      applyImageZoom(preferences_.add_haemorrhages_zoom());
      break;
		case ncm_qscene::ModeDisplayAll:
			applyImageZoom(preferences_.label_image_zoom());
      break;
    default:
      assert(false);
  }
}

//
//: Update toolbar at top of screen
void nailfold_qfinal_gui::updateToolbar()
{
  // Enable/disable buttons depending on whether the image is a dermatogram or
  // capillarogram
  const bool image_is_capillarogram = 
      (imageProcessor_.image_type() == ncm_qimagehandler::TypeCapillarogram);

  const bool user_is_advanced = preferences_.is_advanced_user();

  ui.actionAddVessels->setEnabled(image_is_capillarogram);
  ui.actionDefineVesselSize->setEnabled(image_is_capillarogram);
  ui.actionDefineVesselShape->setEnabled(image_is_capillarogram);
  ui.actionLabelApices->setEnabled(image_is_capillarogram);
  ui.actionDrawVesselPath->setEnabled(image_is_capillarogram && 
                                      user_is_advanced);
  ui.actionAddHaemorrhages->setEnabled(image_is_capillarogram);

  if (ui.graphicsView->isPanning())
  {
    ui.actionPan->setChecked(true);
  }
	else if (scene_.edit_mode() == ncm_qscene::ModeDisplayAll)
	{
	}
	else
  {
    // Action should be offset by 1 to account for PanMode being the first
    int actionIndex = static_cast<int>(scene_.edit_mode());
    QAction* action = modeGroup_->actions()[actionIndex + 1];
    action->setChecked(true);
  }
}

//
//: Update what graphics items are displayed in the scene
void nailfold_qfinal_gui::updateItemVisibility()
{
  // alises that define what is and is not visible in the scene
  ncm_qvesselitem_appearance& v_app = QGraphicsVesselItem::appearance();
  ncm_qapexitem_appearance& a_app = QGraphicsApexItem::appearance();
  ncm_qhaemorrhageitem_appearance& h_app = QGraphicsHaemorrhageItem::appearance();

  bool show_vessel_placeholders = false;
  bool show_vessel_paths = false;
  bool show_apices = false;
  bool show_haemorrhage_placeholders = false;

  switch (scene_.edit_mode())
  {
    case ncm_qscene::ModeClassifyImage:
      // change nothing
      break;
    case ncm_qscene::ModeLabelApices:
      show_apices = true;
      // fall through
    case ncm_qscene::ModeAddVessels:
      // fall through
    case ncm_qscene::ModeSetVesselSize:
      // fall through
    case ncm_qscene::ModeSetVesselShape:
      show_vessel_placeholders = true;
      break;
    case ncm_qscene::ModeDrawVesselPath:
      show_vessel_paths = true;
      break;
    case ncm_qscene::ModeAddHaemorrhages:
      show_haemorrhage_placeholders = true;
      break;

		case ncm_qscene::ModeDisplayAll:
			if (markup_.version() == markup_.latest_auto_version())
			{
				show_haemorrhage_placeholders = ui.show_auto_haemorrhages_checkBox->isChecked();
				show_vessel_paths = ui.show_auto_paths_checkBox->isChecked();
				show_vessel_placeholders = ui.show_auto_placeholders_checkBox->isChecked();
				show_apices = ui.show_auto_apices_checkBox->isChecked();
			}
			else
			{
				show_haemorrhage_placeholders = ui.show_vessel_haemorrhages_checkBox->isChecked();
				show_vessel_paths = ui.show_vessel_paths_checkBox->isChecked();
				show_vessel_placeholders = ui.show_vessel_placeholders_checkBox->isChecked();
				show_apices = ui.show_vessel_apices_checkBox->isChecked();
			}
      
      break;
    default:
      assert(false);
      break;
  }
	vcl_cout << "Vessel placeholders " << show_vessel_placeholders << vcl_endl;
  vcl_cout << "Vessel paths " << show_vessel_paths << vcl_endl;
  vcl_cout << "Vessel apices " << show_apices << vcl_endl;
  vcl_cout << "Vessel haemorrhages " << show_haemorrhage_placeholders << vcl_endl;

  v_app.set_show_placeholder(show_vessel_placeholders);
  v_app.set_show_path(show_vessel_paths);
  a_app.set_show_apex(show_apices);
  h_app.set_show_placeholder(show_haemorrhage_placeholders);
}

//
//: Update brightness and contrast that may change as a result of autocontrast
void nailfold_qfinal_gui::updateContrastControls()
{
  const int min_contrast = imageProcessor_.min_contrast();
  const int max_contrast = imageProcessor_.max_contrast();
  const int brightness = 255 - (min_contrast + max_contrast)/2;
  const int contrast = 255 - max_contrast + min_contrast;

  ui.brightnessSlider->setValue(brightness);
  ui.contrastSlider->setValue(contrast);

  ui.brightnessSlider->setEnabled(!ui.autoContrast->isChecked());
  ui.contrastSlider->setEnabled(!ui.autoContrast->isChecked());
}
 
//: Update vessel selection controls, depending on how many items are in the
//  list and which one is selected
void nailfold_qfinal_gui::updateVesselControls()
{
  // Enable/disable the entire set of controls accordingly
  switch (scene_.edit_mode())
  {
    case ncm_qscene::ModeLabelApices: 
      // fall through
    case ncm_qscene::ModeDrawVesselPath:
      ui.vesselsGroup->setEnabled(true);
      break;

    default:
      ui.vesselsGroup->setEnabled(false);
  }

  const int nVessels = markup_.n_vessels();
  const bool vesselsExist = (nVessels > 0);
  ui.deleteAllVessels->setEnabled(vesselsExist);

  const int nDistal = markup_.n_distal();
  const bool distalVesselsExist = (nDistal > 0);
  ui.applySizeToAll->setEnabled(distalVesselsExist);
  ui.applyShapeToAll->setEnabled(distalVesselsExist);

  // if no vessels then disable all controls
  if (!distalVesselsExist)
  {
    ui.vesselCombo->setEnabled(false);
    ui.firstVesselButton->setEnabled(false);
    ui.prevVesselButton->setEnabled(false);
    ui.nextVesselButton->setEnabled(false);
    ui.lastVesselButton->setEnabled(false);
    return;
  }

  const int firstVessel = 0;
  const int lastVessel = nDistal-1;
  const int currentVessel = ui.vesselCombo->currentIndex();

  // All buttons enabled by default
  bool firstEnabled = true;
  bool prevEnabled = true;
  bool nextEnabled = true;
  bool lastEnabled = true;

  // If item selected is first then disable 'first' and 'previous' buttons
  if (currentVessel == firstVessel)
  {
    firstEnabled = false;
    prevEnabled = false;
  }

  // If item selected is last then disable 'last' and 'next' buttons
  if (currentVessel == lastVessel)
  {
    lastEnabled = false;
    nextEnabled = false;
  }

  // Apply settings
  ui.vesselCombo->setEnabled(true);
  ui.firstVesselButton->setEnabled(firstEnabled);
  ui.prevVesselButton->setEnabled(prevEnabled);
  ui.nextVesselButton->setEnabled(nextEnabled);
  ui.lastVesselButton->setEnabled(lastEnabled);

  // Q: Do we disable the combobox if there is only one vessel? It technically
  // wouldn't be needed (there's only one vessel you can select) but still has
  // the effect of zooming in on the vessel.
}

//
//: Update limits on zoom controls
void nailfold_qfinal_gui::updateZoomLimits()
{
  const double minimum_scale = ui.graphicsView->fitHeightScale();
  ui.zoomSlider->setMinimum(static_cast<int>(100*minimum_scale + 1));

  const double maximum_scale = 10 * minimum_scale;
  ui.zoomSlider->setMaximum(static_cast<int>(100*maximum_scale - 1));

  ui.graphicsView->setMaxScale(maximum_scale);
}

//
//: Update zoom controls
void nailfold_qfinal_gui::updateZoomControls()
{
  const double scale = ui.graphicsView->transform().m11();
  ui.zoomSlider->setValue(static_cast<int>(100*scale + 1));

  QString zoom_string = QString::number(static_cast<int>(100*scale + 1));
  ui.zoomEdit->setText(zoom_string+"%");
}

//
//: Make sure help panel is up to date
void nailfold_qfinal_gui::updateHelpPanel()
{
  // Display correct help page for this mode
  if (ui.graphicsView->isPanning())
  {
    const int lastPage = ui.helpStackedWidget->count()-1;
    ui.helpStackedWidget->setCurrentWidget(ui.panHelpPage);
    return;
  }

  switch (scene_.edit_mode())
  {
  case ncm_qscene::ModeClassifyImage:
    //if (ui.gradeUndefinedRadio->isChecked())
    //  ui.helpStackedWidget->setCurrentWidget(ui.labelImageHelpPage);
    //else 
    if (ui.gradeNormalRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.labelImageHelpPage);
    else if (ui.gradeEarlyRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.labelImageHelpPage);
    else if (ui.gradeActiveRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.labelImageHelpPage);
    else if (ui.gradeLateRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.labelImageHelpPage);
    else if (ui.gradeExtremeRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.labelImageHelpPage);
    else if (ui.gradePoorQualityRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.labelImageHelpPage);
    else
      ui.helpStackedWidget->setCurrentWidget(ui.labelImageHelpPage);
    break;

  case ncm_qscene::ModeAddVessels:
    if (ui.distalRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.addDistalVesselsHelpPage);
    else if (ui.nondistalRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.addNondistalVesselsHelpPage);
    else
      assert(false);
    break;

  case ncm_qscene::ModeSetVesselSize:
    if (ui.sizeNormalRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.normalSizeVesselHelpPage);
    else if (ui.sizeEnlargedRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.enlargedVesselHelpPage);
    else if (ui.sizeGiantRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.giantVesselHelpPage);
    else if (ui.sizeIrregularRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.irregularVesselHelpPage);
    else
      assert(false);
    break;

  case ncm_qscene::ModeSetVesselShape:
    if (ui.shapeNormalRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.normalShapeVesselHelpPage);
    else if (ui.shapeTortuousRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.meanderingVesselHelpPage);
    else if (ui.shapeRamifiedRadio->isChecked())
      ui.helpStackedWidget->setCurrentWidget(ui.angiogenicVesselHelpPage);
    else
      assert(false);
    break;

  case ncm_qscene::ModeLabelApices:
    ui.helpStackedWidget->setCurrentWidget(ui.labelApicesHelpPage);
    break;

  case ncm_qscene::ModeDrawVesselPath:
    ui.helpStackedWidget->setCurrentWidget(ui.drawVesselPathHelpPage);
    break;

  case ncm_qscene::ModeAddHaemorrhages:
    ui.helpStackedWidget->setCurrentWidget(ui.addHaemorrhagesHelpPage);
    break;

	case ncm_qscene::ModeDisplayAll:
		break;

  default:
    assert(false);
  }
}

//
//: Update remaining widgets based on information in the markup
void nailfold_qfinal_gui::applyAnnotationChange()
{
  // Update the list of vessels
  updateVesselList();

  updateGradeRadios();

  const bool vesselsExist = (markup_.n_vessels() > 0);
  ui.deleteAllVessels->setEnabled(vesselsExist);

  const bool distalVesselsExist = (markup_.n_distal() > 0);
  ui.applySizeToAll->setEnabled(distalVesselsExist);
  ui.applyShapeToAll->setEnabled(distalVesselsExist);

  const bool haemorrhages_exist = (markup_.n_haemorrhages() > 0);
  ui.deleteAllHaemorrhages->setEnabled(haemorrhages_exist);
}

//
//: Define the default path in which to look for images (can be set from
//  command line parameters)
void nailfold_qfinal_gui::setDefaultPath(vcl_string image_path)
{
  currentImagePath_ = image_path;
}

//
//  Private slots
//

//: Update buttons that cycle through vessel row
void nailfold_qfinal_gui::updateRowButtons(int selected)
{
  int nRows = rowRadioGroup_->buttons().size();

  const bool selected_is_first = (selected == 0);
  ui.prevRow->setEnabled(!selected_is_first);

  // If selection is the last in the list then treat the "Next *" button as if
  // it were "Next Stage".
  const bool selected_is_last = (selected == nRows-1);
  if (selected_is_last)
  {
    ui.nextRow->setText("Next Stage");
    QObject::disconnect(ui.nextRow, 0, 0, 0);
    QObject::connect(ui.nextRow, SIGNAL(clicked()), 
                     this, SLOT(on_nextStageButton_clicked()) );
  }
  else
  {
    ui.nextRow->setText("Next Row");
    QObject::disconnect(ui.nextRow, 0, 0, 0);
    QObject::connect(ui.nextRow, SIGNAL(clicked()), 
                     this, SLOT(on_nextRow_clicked()) );
  }

  updateHelpPanel();
}
void nailfold_qfinal_gui::on_prevRow_clicked()
{
  // Alias for convenience
  QButtonGroup*& radioGroup = rowRadioGroup_;

  // selected should never be zero, as the button should be disabled in this
  // instance
  int selected = radioGroup->checkedId();
  QAbstractButton* prevRadio = radioGroup->button(selected-1);
  prevRadio->click();
}
void nailfold_qfinal_gui::on_nextRow_clicked()
{
  // Alias for convenience
  QButtonGroup*& radioGroup = rowRadioGroup_;

  // selected should never be zero, as the button should be disabled in this
  // instance
  int selected = radioGroup->checkedId();
  QAbstractButton* nextRadio = radioGroup->button(selected+1);
  nextRadio->click();
}
//
//: Update buttons that cycle through vessel sizes
void nailfold_qfinal_gui::updateSizeButtons(int selected)
{
  const bool selected_is_first = (selected == 0);
  ui.prevSize->setEnabled(!selected_is_first);

  // If selection is the last in the list then treat the "Next *" button as if
  // it were "Next Stage".
  const int lastSize = sizeRadioGroup_->buttons().size() - 1;
  const bool selected_is_last = (selected == lastSize);
  if (selected_is_last)
  {
    ui.nextSize->setText("Next Stage");
    QObject::disconnect(ui.nextSize, 0, 0, 0);
    QObject::connect(ui.nextSize, SIGNAL(clicked()), 
                     this, SLOT(on_nextStageButton_clicked()) );
  }
  else
  {
    ui.nextSize->setText("Next Size");
    QObject::disconnect(ui.nextSize, 0, 0, 0);
    QObject::connect(ui.nextSize, SIGNAL(clicked()), 
                     this, SLOT(on_nextSize_clicked()) );
  }

  updateHelpPanel();
}
void nailfold_qfinal_gui::on_prevSize_clicked()
{
  // selected should never be zero, as the button should be disabled in this
  // instance
  int selected = sizeRadioGroup_->checkedId();
  QAbstractButton* prevRadio = sizeRadioGroup_->button(selected-1);
  prevRadio->click();
}
void nailfold_qfinal_gui::on_nextSize_clicked()
{
  // Selected should never be last, as the button should then associated with
  // "Next Stage" action.
  int selected = sizeRadioGroup_->checkedId();
  QAbstractButton* nextRadio = sizeRadioGroup_->button(selected+1);
  nextRadio->click();
}
//
//: Update buttons that cycle through vessel shapes
void nailfold_qfinal_gui::updateShapeButtons(int selected)
{
  const bool selected_is_first = (selected == 0);
  ui.prevShape->setEnabled(!selected_is_first);

  // If selection is the last in the list then treat the "Next *" button as if
  // it were "Next Stage".
  const int lastShape = shapeRadioGroup_->buttons().size() - 1;
  const bool selected_is_last = (selected == lastShape);
  if (selected_is_last)
  {
    ui.nextShape->setText("Next Stage");
    QObject::disconnect(ui.nextShape, 0, 0, 0);
    QObject::connect(ui.nextShape, SIGNAL(clicked()), 
                     this, SLOT(on_nextStageButton_clicked()) );
  }
  else
  {
    ui.nextShape->setText("Next Shape");
    QObject::disconnect(ui.nextShape, 0, 0, 0);
    QObject::connect(ui.nextShape, SIGNAL(clicked()), 
                     this, SLOT(on_nextShape_clicked()) );
  }

  updateHelpPanel();
}

void nailfold_qfinal_gui::on_prevShape_clicked()
{
  // selected should never be zero, as the button should be disabled in this
  // instance
  int selected = shapeRadioGroup_->checkedId();
  QAbstractButton* prevRadio = shapeRadioGroup_->button(selected-1);
  prevRadio->click();
}
void nailfold_qfinal_gui::on_nextShape_clicked()
{
  // selected should never be zero, as the button should be disabled in this
  // instance
  int selected = shapeRadioGroup_->checkedId();
  QAbstractButton* nextRadio = shapeRadioGroup_->button(selected+1);
  nextRadio->click();
}
//
//: Update stage controls
void nailfold_qfinal_gui::updateStageControls()
{
	if (scene_.edit_mode() == ncm_qscene::ModeDisplayAll)
	{
		//TO DO: what will we do?
		return;
	}
  // Display correct set of controls for this mode
  const int edit_mode = scene_.edit_mode();
  ui.stackedWidget->setCurrentIndex(edit_mode);

  // If we're in pan mode then disable all stage controls altogether.
  if (ui.graphicsView->isPanning())
  {
    ui.stackedWidget->setEnabled(false);
    ui.prevStageButton->setEnabled(false);
    ui.nextStageButton->setEnabled(false);
    return;
  }

  // Otherwise, stacked widget should be enabled
  ui.stackedWidget->setEnabled(true);

  // If there are no images or the image is a dermatogram then disable 
  // prev/next stage buttons altogether.
  if (fileServer_.is_empty() ||
      (imageProcessor_.image_type() == ncm_qimagehandler::TypeDermatogram))
  {
    ui.prevStageButton->setEnabled(false);
    ui.nextStageButton->setEnabled(false);
  }
  else
  {
    // Always disable 'Previous Stage' button if in first mode
    ui.prevStageButton->setEnabled(edit_mode != ncm_qscene::ModeFirst);

    // If in last mode (Add Haemorrhages), 'Next Stage' should mimic the 
    // behaviour of the 'Next Image' button so that the user can go straight 
    // to the next image.
    if (edit_mode == ncm_qscene::ModeLast)
    {
      ui.nextStageButton->setText(ui.nextImage->text());
      ui.nextStageButton->setEnabled(ui.nextImage->isEnabled());
    }
    else
    {
      ui.nextStageButton->setText("Next Stage");
      ui.nextStageButton->setEnabled(true);
    }
  }
}

void nailfold_qfinal_gui::on_prevStageButton_clicked()
{
  const int first_stage = 0;
  if (ui.stackedWidget->currentIndex() != first_stage)
  {
    // Skip any stages whose stackedWidget page is disabled.
    // Used for skipping DrawVesselPath as basic user.
    do 
    {
      requestEditModeChange(ui.stackedWidget->currentIndex()-1);
    } while (!ui.stackedWidget->currentWidget()->isEnabled());
  }
}
void nailfold_qfinal_gui::on_nextStageButton_clicked()
{
  const int last_stage = ui.stackedWidget->count() - 1;

  // In the last stage of marking up an image, change the behaviour of this
  // button (its text is modified in updateStageControls()) to mimic that of
  // the 'Next Image' button for quicker progress.
  if (ui.stackedWidget->currentIndex() == last_stage)
    on_nextImage_clicked();
  else
  {
    // Skip any stages whose stackedWidget page is disabled.
    // Used for skipping DrawVesselPath as basic user.
    do 
    {
      requestEditModeChange(ui.stackedWidget->currentIndex()+1);
    } while (!ui.stackedWidget->currentWidget()->isEnabled());
  }
}
//
//: Update the buttons that progress through the current image directory
void nailfold_qfinal_gui::updateImageControls()
{
  // 'Previous Image' button should be enabled only when there is a previous
  // image to go to (i.e. when the current image is not the first).
  ui.prevImage->setEnabled(!fileServer_.is_first());

  // Whether the 'Next Image' button is enabled (and what it says) depends on
  // how many images are left to mark up.
  const unsigned n_images = fileServer_.count();
  if (n_images > 1)
  {
    // If there are other images then set the button properties accordingly.
    showScene();
    ui.nextImage->setText("Next Image");
    ui.nextImage->setEnabled(!fileServer_.is_last());
  }
  else if (n_images == 1)
  {
    // If this is the last image to be marked (or it has already been marked)
    // then change the button text and make it always enabled.
    showScene();
    ui.nextImage->setText("Save and Finish");
    ui.nextImage->setEnabled(true);
  }
  else // if (n_images == 0)
  {
    // If this is the last image to be marked (or it has already been marked)
    // then change the button text and make it always enabled.
    hideScene();
    ui.nextImage->setText("Next Image");
    ui.nextImage->setEnabled(false);
  }
}

void nailfold_qfinal_gui::on_prevImage_clicked()
{
  // Scene edit mode is changed (without warnings) if the image is
  // successfully loaded. Apply warnings here, before the image is loaded.
  const bool dismiss_warning = warnOnEditModeChange();
  if (!dismiss_warning)
    return;

  // This is here so that the image is flagged as 'marked' before going to the 
  // next file. Not ideal but works in practice.
  if (!saveIfModified())
    return;

  bool found_previous = fileServer_.previous();
  if (found_previous)
    loadImageFrom(fileServer_.current_image());
  else
  {
    // This was the first image. Just update the image controls.
    updateImageControls();
  }
}
void nailfold_qfinal_gui::on_nextImage_clicked()
{
  // Scene edit mode is changed (without warnings) if the image is
  // successfully loaded. Apply warnings here, before the image is loaded.
  const bool dismiss_warning = warnOnEditModeChange();
  if (!dismiss_warning)
    return;

  if (!fileServer_.is_last())
  {
    // There are more images to come, so check if this needs saving then go
    // to next image.
    if (!saveIfModified())
      return;

    bool found_next = fileServer_.next();
    loadImageFrom(fileServer_.current_image());
  }
  else
  {
    // This was the last image. Close the application.
    close();
  }
}

void nailfold_qfinal_gui::on_unmarkedCheckbox_toggled(bool checked)
{
  fileServer_.use_unmarked_only(checked);

  if (fileServer_.current_image() == "")
  {
    // If no image was loaded then load the first one.
    // This can happen when all images have been marked then we switch to show
    // all images.
    fileServer_.first();
    loadImageFrom(fileServer_.current_image());
  }
  else if (fileServer_.is_using_unmarked_only() &&
           fileServer_.is_current_marked())
  {
    // If we switch to unmarked images and the current image is marked then
    // look for an unmarked image to show. If there are no such images then
    // clear the image altogether.
    if (fileServer_.is_empty())
    {
      updateImageControls();
      updateStageControls();
    }
    else
    {
      // Look forward.
      bool found = fileServer_.next();

      // If nothing there then look back.
      if (!found)
        fileServer_.previous();

      // Load the image found (there must be one as we've already checked that
      // the image queue is not empty
      loadImageFrom(fileServer_.current_image());
    }
  }
  else
  {
    // Just update the stage and image controls
    updateImageControls();
    updateStageControls();
  }

	//Need to update thumbnails even though we've not changed folders
  getThumbnails();
}

//
//: Toolbar actions
void nailfold_qfinal_gui::on_actionPan_triggered()
{
  ui.graphicsView->usePanMode();
}
void nailfold_qfinal_gui::on_actionLabelImage_triggered()
{
  ui.graphicsView->useMarkupMode();
  requestEditModeChange(ncm_qscene::ModeClassifyImage);
}
void nailfold_qfinal_gui::on_actionAddVessels_triggered()
{
  ui.graphicsView->useMarkupMode();
  requestEditModeChange(ncm_qscene::ModeAddVessels);
}
void nailfold_qfinal_gui::on_actionDefineVesselSize_triggered()
{
  ui.graphicsView->useMarkupMode();
  requestEditModeChange(ncm_qscene::ModeSetVesselSize);
}
void nailfold_qfinal_gui::on_actionDefineVesselShape_triggered()
{
  ui.graphicsView->useMarkupMode();
  requestEditModeChange(ncm_qscene::ModeSetVesselShape);
}
void nailfold_qfinal_gui::on_actionLabelApices_triggered()
{
  ui.graphicsView->useMarkupMode();
  requestEditModeChange(ncm_qscene::ModeLabelApices);
}
void nailfold_qfinal_gui::on_actionDrawVesselPath_triggered()
{
  ui.graphicsView->useMarkupMode();
  requestEditModeChange(ncm_qscene::ModeDrawVesselPath);
}
void nailfold_qfinal_gui::on_actionAddHaemorrhages_triggered()
{
  ui.graphicsView->useMarkupMode();
  requestEditModeChange(ncm_qscene::ModeAddHaemorrhages);
}

//
//:
void nailfold_qfinal_gui::on_applySizeToAll_clicked()
{
  unsigned n_distal = markup_.n_distal();
  for (unsigned i = 0; i < n_distal; ++i)
  {
    markup_.vessel(i)->properties().copy_size_from(
        scene_.new_vessel_properties()
    );
  }

  scene_.update();
  updateVesselList();
}
void nailfold_qfinal_gui::on_applyShapeToAll_clicked()
{
  // FIXME

  for (unsigned i = 0; i < markup_.n_distal(); ++i)
  {
    markup_.vessel(i)->properties().copy_shape_from(
        scene_.new_vessel_properties()
    );
  }

  scene_.update();
  updateVesselList();
}

//
//: Update list of vessels created
void nailfold_qfinal_gui::updateVesselList()
{
  // Annoyingly, adding an item to a combobox also sets the currentIndex
  // (unlike in a list, that leaves the currentRow unchanged). Therefore,
  // whenever you add a new vessel, currentIndexChanged() tells the 
  // graphicsView to zoom into the new vessel. To avoid this, I temporarily
  // disconnect the signals and slots until I have a chance to revert to the 
  // previously selected item.
  QObject::disconnect( ui.vesselCombo, SIGNAL(currentIndexChanged(int)),
                       &scene_, SLOT(select_vessel(int)) );

  // Store currently selected index in combobox.
  const int currentIndex = ui.vesselCombo->currentIndex();

  ui.vesselCombo->clear();
  for (unsigned i = 0; i < markup_.n_distal(); ++i)
  {
    QString s = "Vessel " + QString::number(i+1) + " (" +
                markup_.vessel(i)->properties().size_string().c_str() + "; " +
                markup_.vessel(i)->properties().shape_string().c_str() + ")";

    ui.vesselCombo->addItem(s);
  }

  // Restore previous selection. If this no longer exists, it will return to 
  // a null selection (currentIndex = -1)
  ui.vesselCombo->setCurrentIndex(currentIndex);

  // Reconnect signals
  QObject::connect( ui.vesselCombo, SIGNAL(currentIndexChanged(int)),
                    &scene_, SLOT(select_vessel(int)) );

  // If we're in a zoomed view then reselect the vessel and update the view
  // accordingly
  if (scene_.edit_mode() == ncm_qscene::ModeLabelApices ||
      scene_.edit_mode() == ncm_qscene::ModeDrawVesselPath)
  {
    scene_.select_vessel(currentIndex);
    applyVesselZoom(currentIndex);
  }

  updateVesselControls();
}

//
//: Select nothing
void nailfold_qfinal_gui::clearVesselListSelection()
{
  ui.vesselCombo->setCurrentIndex(-1);
}

//
//: Change the vessel selection in the combobox
void nailfold_qfinal_gui::on_firstVesselButton_clicked()
{
  const int nVessels = markup_.n_distal();

  if (nVessels > 0)
    scene_.select_vessel(0);
}

void nailfold_qfinal_gui::on_prevVesselButton_clicked()
{
  const int nVessels = markup_.n_distal();
  const int currentIndex = ui.vesselCombo->currentIndex();

  if (currentIndex > 0)
    scene_.select_vessel(currentIndex-1);
}

void nailfold_qfinal_gui::on_nextVesselButton_clicked()
{
  const int nVessels = markup_.n_distal();
  const int currentIndex = ui.vesselCombo->currentIndex();

  if (currentIndex < nVessels-1)
    scene_.select_vessel(currentIndex+1);
}

void nailfold_qfinal_gui::on_lastVesselButton_clicked()
{
  const int nVessels = markup_.n_distal();

  if (nVessels > 0)
    scene_.select_vessel(nVessels-1);
}

//
//  Zoom controls
//
void nailfold_qfinal_gui::setZoomFrom(int value)
{
  const double minimum_scale = ui.graphicsView->fitScale();
  double new_scale = static_cast<double>(value-1) / 100;

  if (new_scale < minimum_scale)
    new_scale = minimum_scale;

  ui.graphicsView->setScale(new_scale);
}
void nailfold_qfinal_gui::on_zoomSlider_sliderMoved(int value)
{
  setZoomFrom(value);
}

void nailfold_qfinal_gui::on_zoomEdit_editingFinished()
{
  // Remove last character (%) from string and convert
  QString value_string = ui.zoomEdit->text();
  value_string.chop(1);
  const int value = value_string.toInt();

  setZoomFrom(value);
}

void nailfold_qfinal_gui::on_zoomOutButton_clicked()
{
  // Avoid calling functions that generate more signals, such as this:
  //ui.zoomSlider->setValue(ui.zoomSlider->value()-ui.zoomSlider->pageStep());

  // Instead, set the graphicsView zoom directly, and let the resulting signal
  // trigger an update of the GUI controls
  int new_zoom = ui.zoomSlider->value() - ui.zoomSlider->pageStep();
  setZoomFrom(new_zoom);
}
void nailfold_qfinal_gui::on_zoomInButton_clicked()
{
  // Avoid calling functions that generate more signals, such as this:
  //ui.zoomSlider->setValue(ui.zoomSlider->value()+ui.zoomSlider->pageStep());

  // Instead, set the graphicsView zoom directly, and let the resulting signal
  // trigger an update of the GUI controls
  int new_zoom = ui.zoomSlider->value() + ui.zoomSlider->pageStep();
  setZoomFrom(new_zoom);
}

//
//: Button to hide/show side panel with help and thumbnails
void nailfold_qfinal_gui::on_hidePanelButton_clicked()
{
  if (ui.tabControl->isVisible())
  {
    ui.tabControl->hide();
    ui.hidePanelButton->setText("<");
  }
  else
  {
    ui.tabControl->show();
    ui.hidePanelButton->setText(">");
  }
}

//
//: Change what is visible when displaying (not editing) markup
void nailfold_qfinal_gui::on_show_auto_haemorrhages_checkBox_toggled()
{
	updateItemVisibility();
	scene_.update();
}

void nailfold_qfinal_gui::on_show_auto_paths_checkBox_toggled()
{
	updateItemVisibility();
	scene_.update();
}

void nailfold_qfinal_gui::on_show_auto_placeholders_checkBox_toggled()
{
	updateItemVisibility();
	scene_.update();
}

void nailfold_qfinal_gui::on_show_auto_apices_checkBox_toggled()
{
	updateItemVisibility();
	scene_.update();
}

void nailfold_qfinal_gui::on_show_vessel_haemorrhages_checkBox_toggled()
{
	updateItemVisibility();
}

void nailfold_qfinal_gui::on_show_vessel_paths_checkBox_toggled()
{
	updateItemVisibility();
}

void nailfold_qfinal_gui::on_show_vessel_placeholders_checkBox_toggled()
{
	updateItemVisibility();
}

void nailfold_qfinal_gui::on_show_vessel_apices_checkBox_toggled()
{
	updateItemVisibility();
}


//
//: Sound warning when making a transition from a given scene edit mode.
//  This is used to avoid accidentally missing annotations.
//  Returns true unless user elects to heed the warning.
bool nailfold_qfinal_gui::warnOnEditModeChange()
{
  // Issue warnings where appropriate before applying any change.
  QMessageBox msgBox;
  msgBox.setIcon(QMessageBox::Warning);
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  msgBox.setDefaultButton(QMessageBox::No);
  msgBox.setInformativeText("Continue?");

  unsigned missing_flags = markup_.missing_data();

  switch (scene_.edit_mode())
  {
    case ncm_qscene::ModeClassifyImage:
      if (preferences_.warn_if_ungraded)
      {
        // Warn if image has not been graded.
        if ((missing_flags & ncm_annotation::MissingData_Grade) !=
              ncm_annotation::MissingData_None)
        {
          msgBox.setStandardButtons(QMessageBox::Ok);
          msgBox.setDefaultButton(QMessageBox::Ok);
          msgBox.setInformativeText("");
          msgBox.setText("You must give the image a grade.");
          msgBox.exec();
          return false; // Do not allow user to continue if image is ungraded.
        }
      }
      break;

    case ncm_qscene::ModeAddVessels:
      if (preferences_.warn_if_vessels_missing)
      {
        // Warn if no vessels were added.
        if ((missing_flags & ncm_annotation::MissingData_Vessels) !=
              ncm_annotation::MissingData_None)
        {
          msgBox.setText("No vessels tagged");
          return (msgBox.exec() == QMessageBox::Yes);
        }
      }
      break;

    case ncm_qscene::ModeSetVesselSize:
      if (preferences_.warn_if_sizes_missing)
      {
        // Warn if some vessels not labelled with a size.
        if ((missing_flags & ncm_annotation::MissingData_Sizes) !=
              ncm_annotation::MissingData_None)
        {
          msgBox.setText("Vessels have undefined size");
          return (msgBox.exec() == QMessageBox::Yes);
        }
      }
      break;

    case ncm_qscene::ModeSetVesselShape:
      if (preferences_.warn_if_shapes_missing)
      {
        // Warn if some vessels not labelled with a shape.
        if ((missing_flags & ncm_annotation::MissingData_Shapes) !=
              ncm_annotation::MissingData_None)
        {
          msgBox.setText("Vessels have undefined shape");
          return (msgBox.exec() == QMessageBox::Yes);
        }
      }
      break;

    case ncm_qscene::ModeLabelApices:
      if (preferences_.warn_if_apices_missing)
      {
        // Warn if some vessels not labelled with an apex.
        if ((missing_flags & ncm_annotation::MissingData_Apices) !=
              ncm_annotation::MissingData_None)
        {
          msgBox.setText("Vessels have undefined apices");
          return (msgBox.exec() == QMessageBox::Yes);
        }
      }
      break;

    case ncm_qscene::ModeDrawVesselPath:
      // Break if not an advanced user since we shouldn't even be here in the
      // first place.
      if (!preferences_.is_advanced_user())
        break;

      if (preferences_.warn_if_paths_missing)
      {
        // Warn if some vessels not labelled with a path.
        if ((missing_flags & ncm_annotation::MissingData_Paths) !=
              ncm_annotation::MissingData_None)
        {
          msgBox.setText("Vessels have undefined paths");
          return (msgBox.exec() == QMessageBox::Yes);
        }
      }
      break;

    case ncm_qscene::ModeAddHaemorrhages:
      if (preferences_.warn_if_haemorrhages_missing)
      {
        // Warn if no haemorrhages tagged
        if ((missing_flags & ncm_annotation::MissingData_Haemorrhages) !=
              ncm_annotation::MissingData_None)
        {
          msgBox.setText("No haemorrhages tagged");
          return (msgBox.exec() == QMessageBox::Yes);
        }
      }
      break;

		case ncm_qscene::ModeDisplayAll:
				//Do nothing
				break;

    default:
      assert(false);
  }

  return true;
}


//
//: Process a request to change edit mode.
void nailfold_qfinal_gui::requestEditModeChange(int newMode)
{
  // If the new mode is the same as the current one then do nothing.
  if (static_cast<ncm_qscene::editMode>(newMode) == scene_.edit_mode())
    return;

  const bool dismiss_warning = warnOnEditModeChange();
  if (dismiss_warning)
    scene_.set_edit_mode(newMode);
}

//
//: Make widgets consistent with scene's edit mode (typically called after a
//  confirmed change in edit mode).
void nailfold_qfinal_gui::on_scene_editModeChanged(int edit_mode)
{
  // Set up mode-specific controls and select an appropriate vessel (usually
  // none).
  switch (scene_.edit_mode())
  {
    case ncm_qscene::ModeClassifyImage:
      //ui.gradeUndefinedRadio->setChecked(true);
      updateGradeRadios();
      scene_.select_vessel(-1);
      break;

    case ncm_qscene::ModeAddVessels:
      ui.distalRadio->click();
      scene_.new_vessel_properties().setSizeUndefined(true);
      scene_.new_vessel_properties().setShapeUndefined(true);
      scene_.select_vessel(-1);
      break;

    case ncm_qscene::ModeSetVesselSize:
      ui.sizeNormalRadio->click();
      scene_.new_vessel_properties().setSizeNormal(true);
      scene_.select_vessel(-1);
      break;

    case ncm_qscene::ModeSetVesselShape:
      ui.shapeNormalRadio->click();
      scene_.new_vessel_properties().setShapeNormal(true);
      scene_.select_vessel(-1);
      break;

    case ncm_qscene::ModeLabelApices:
      // fall through
    case ncm_qscene::ModeDrawVesselPath:
      // Select the first vessel when labelling its specific properties
      scene_.select_vessel(0);
      break;

    case ncm_qscene::ModeAddHaemorrhages:
      scene_.select_vessel(-1);
      break;

		case ncm_qscene::ModeDisplayAll:
      scene_.select_vessel(-1);
			ui.graphicsView->usePanMode();
      break;

    default:
      // Covered all options - we should never reach this point.
      assert(false);
  }

  updateItemVisibility();
	updateToolbar();
	updateHelpPanel();

  // Update various GUI elements
	if (scene_.edit_mode() != ncm_qscene::ModeDisplayAll)
	{
		
		updateProgressFrame();		
		updateVesselControls();
		updateStageControls();
	}

  applyUserZoom();
}

void nailfold_qfinal_gui::on_scene_mouseMoved(double x, double y)
{
  ui.zoomedView->centerOn(x, y);
}

//
// Events
//
void nailfold_qfinal_gui::closeEvent(QCloseEvent *ev)
{
  // only close if (a) markup is not modified or (b) user declines to save or 
  // (c) user successfully saves to a file.
  if (saveIfModified())
    ev->accept();
  else
    ev->ignore();
}

void nailfold_qfinal_gui::getThumbnails()
{
  // Create alias for convenience
  QBoxLayout* thumbs = ui.thumbnails_layout;

  // Clear existing thumbnails
  QLayoutItem *item = thumbs->takeAt(0);
  while (item != 0)
  {
    // We need to delete both the widget and item separately here.
    // This may be a result of creating our own widgets - I'm not sure.
    delete item->widget();
    delete item;
    item = thumbs->takeAt(0);
  }

  // Create new thumbnails from currentImageDir_
  const int pixmap_width = 160;

  vcl_vector<vcl_string> images = fileServer_.valid_images();
  unsigned n_images = images.size();
	unsigned max_images = 100;

	if (n_images > max_images)
		n_images = max_images;

  updateStatus("Creating thumbnails...", 0.0);
  for (unsigned i = 0; i < n_images; ++i)
  {
    QImage qimg(images[i].c_str());

    QPixmap qpixmap;
    qpixmap.convertFromImage(qimg.scaledToWidth(pixmap_width));

    ncm_qthumbnail *qlabel = new ncm_qthumbnail(qpixmap, images[i].c_str());
    thumbs->addWidget(qlabel);
    QObject::connect( qlabel, SIGNAL(clicked(vcl_string)),
                      this, SLOT(onThumbnail_clicked(vcl_string)) );
    updateStatus("Creating thumbnails..", static_cast<double>(i+1)/n_images);
  }
  thumbs->addStretch();
  updateStatus("");
}

//
//: Load an image when its thumbnail is clicked
void nailfold_qfinal_gui::onThumbnail_clicked(vcl_string filename)
{
  // Scene edit mode is changed (without warnings) if the image is
  // successfully loaded. Apply warnings here, before the image is loaded.
  const bool dismiss_warning = warnOnEditModeChange();
  if (!dismiss_warning)
    return;

  // There are more images to come, so check if this needs saving then go
  // to next image.
  if (!saveIfModified())
    return;

  loadImageFrom(filename);
}

//
//  File I/O
//

//
//: Delete the temporary backup file
bool nailfold_qfinal_gui::deleteBackup() const
{
  vcl_string backup_filename = markup_.snapshot_filename();
  return QFile(backup_filename.c_str()).remove();
}

//
//: Load an image from file
bool nailfold_qfinal_gui::loadImageFrom(vcl_string image_filename)
{
  // Ignore if no image name is supplied
  if (image_filename == "")
    return false;

  // If the markup has been modified, ask user if they want to save
  // If they 'Cancel', return without further action
  if (!saveIfModified())
    return false;

	//Update the image path
	vcl_string dirname = vul_file::dirname(image_filename);
	updateImagePath(dirname);

  // Update the window title
  setWindowTitle(windowTitle() + " - " + image_filename.c_str());

  // Clear scene (including markup) before we load anything
  // Turn off snapshots temporarily otherwise we'll create a new backup file
  markup_.set_snapshots(false);
  scene_.clear();
  markup_.set_snapshots(true);

  // Tell the image processor which file to load then load it.
  imageProcessor_.set_filename(image_filename);
  imageProcessor_.load_image();

  fileServer_.set_current(image_filename);

  // Define markup filename and check if it exists
  const vcl_string markup_filename = 
      vul_file::strip_extension(image_filename) + "_markup.txt";
  const bool markup_exists = vul_file::exists(markup_filename);

  // Define backup filename check if it exists
  const vcl_string backup_filename = markup_filename + '~';
  const bool backup_exists = vul_file::exists(backup_filename);

  bool markup_loaded = false;

  // Check for an autosaved markup that still exists (indicative of a crash)
  if (backup_exists)
  {
    // if markup exists then compare timestamps to check if backup is newer
    bool backup_is_newer = false;
    if (markup_exists)
    {
      QFileInfo fileinfo;

      fileinfo.setFile(QString(backup_filename.c_str()));
      QDateTime backup_timestamp = fileinfo.lastModified();
      fileinfo.setFile(QString(markup_filename.c_str()));
      QDateTime markup_timestamp = fileinfo.lastModified();

      backup_is_newer = (markup_timestamp < backup_timestamp);
    }

    // If backup is there but not the 'non-backup', or if backup is newer than
    // the 'non-backup' then ask user if they want to restore from the backup
    // instead
    if (!markup_exists || backup_is_newer)
    {
      QMessageBox msgBox;
      msgBox.setText("Restore From Backup?");
      msgBox.setInformativeText("Restore previous session from backup?");
      msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
      msgBox.setDefaultButton(QMessageBox::Yes);
      const int response = msgBox.exec();

      if (response == QMessageBox::Yes)
      {
        // Load markup and apply it to the scene
        markup_.t_read(backup_filename);
        markup_.set_filename(markup_filename);
        saveMarkup(); // replace old file with restored backup
        scene_.update_from_annotation();
        markup_loaded = true;
      }
    }
  }

  if (!markup_loaded && markup_exists)
  {
    // Load markup and apply it to the scene
    markup_.t_read(markup_filename);
    markup_.set_filename(markup_filename);
    scene_.update_from_annotation();
    markup_loaded = true;
  }

  // If markup still not loaded then give it a blank filename so that the user
  // will be prompted when they ask to 'save'.
  if (!markup_loaded)
    markup_.set_filename("");

  // This determines the backup filename
  markup_.set_image_filename(image_filename);

  // Restart the timer
  markup_.start_timing();

  // Start new image in 'Grade Image' mode
  applyAnnotationChange();
	scene_.set_edit_mode(ncm_qscene::ModeDisplayAll); // Warnings can be skipped here
  updateImageControls();
  updateStatusLeft();

	//Filter the markup list view to show markup for this image
	vcl_string file_name = "*" + vul_file::strip_extension(vul_file::strip_directory(image_filename)) + "*.txt";
	//QRegExp rx(QString(file_name.c_str()));
	ui.markupListView->setRootIndex(markupFilterModel_->mapFromSource(baseFileModel_->index(QString(dirname.c_str()))));
	//markupFilterModel_->setFilterWildcard(QString("*.txt"));
	markupFilterModel_->setFilterWildcard(QString(file_name.c_str()));
	//markupFilterModel_->setFilterRegExp(rx);
	
  
  return true;
}
bool nailfold_qfinal_gui::selectCaseDialog() 
{
	QString qdir_path = QFileDialog::getExistingDirectory(
			this, tr("Select case to analyse"), rootCasePath_.c_str(),
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

		if (!qdir_path.isEmpty())
			updateCasePath(qdir_path.toStdString());

		return true;

}

bool nailfold_qfinal_gui::loadImageDialog()
{
  // get filename from dialog box
	QString selfilter;
	QString fileName = QFileDialog::getOpenFileName(
    this, "Load image", currentCasePath_.c_str(),
    tr("All Images (*.jpg *.png *.bmp)" ), // filter
    &selfilter );

  // if cancelled, do nothing (return)
  if (fileName.isEmpty())
    return false;

  // load the image
  return loadImageFrom(fileName.toStdString());
}

//
//: Load an image from file
bool nailfold_qfinal_gui::loadMarkupFrom(vcl_string markup_filename)
{
  // Ignore if no image name is supplied
  if (markup_filename == "")
    return false;

	// Check markup file exists
  const bool markup_exists = vul_file::exists(markup_filename);

	if (!markup_exists)
		return false;

  // If the markup has been modified, ask user if they want to save
  // If they 'Cancel', return without further action
  if (!saveIfModified())
    return false;

  // Clear scene (including markup) before we load anything
  // Turn off snapshots temporarily otherwise we'll create a new backup file
  markup_.set_snapshots(false);
  scene_.clear_markup();
  markup_.set_snapshots(true);

  // Define backup filename check if it exists
  const vcl_string backup_filename = markup_filename + '~';
  const bool backup_exists = vul_file::exists(backup_filename);

  bool markup_loaded = false;

  // Check for an autosaved markup that still exists (indicative of a crash)
  if (backup_exists)
  {
    // if markup exists then compare timestamps to check if backup is newer
    bool backup_is_newer = false;

    QFileInfo fileinfo;

    fileinfo.setFile(QString(backup_filename.c_str()));
    QDateTime backup_timestamp = fileinfo.lastModified();
    fileinfo.setFile(QString(markup_filename.c_str()));
    QDateTime markup_timestamp = fileinfo.lastModified();

    backup_is_newer = (markup_timestamp < backup_timestamp);

    // If backup is there but not the 'non-backup', or if backup is newer than
    // the 'non-backup' then ask user if they want to restore from the backup
    // instead
    if (backup_is_newer)
    {
      QMessageBox msgBox;
      msgBox.setText("Restore From Backup?");
      msgBox.setInformativeText("Restore previous session from backup?");
      msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
      msgBox.setDefaultButton(QMessageBox::Yes);
      const int response = msgBox.exec();

      if (response == QMessageBox::Yes)
      {
        // Load markup and apply it to the scene
        markup_.t_read(backup_filename);
        markup_.set_filename(markup_filename);
        saveMarkup(); // replace old file with restored backup
        scene_.update_from_annotation();
        markup_loaded = true;
      }
    }
  }

  if (!markup_loaded)
  {
    // Load markup and apply it to the scene
    markup_.t_read(markup_filename);
    markup_.set_filename(markup_filename);
    scene_.update_from_annotation();
		updateMarkupPath(vul_file::dirname(markup_filename));
    markup_loaded = true;
  }

  // Restart the timer
  markup_.start_timing();

  // Start new image in 'Display' mode
  applyAnnotationChange();
  scene_.set_edit_mode(ncm_qscene::ModeDisplayAll); // Warnings can be skipped here
  updateImageControls();
  updateStatusLeft();
  
  return true;
}

bool nailfold_qfinal_gui::loadMarkupDialog()
{
  // get filename from dialog box
	QString selfilter;
	QString fileName = QFileDialog::getOpenFileName(
    this, "Load image", currentCasePath_.c_str(),
    tr("Markup files (*.txt)" ), // filter
    &selfilter );

  // if cancelled, do nothing (return)
  if (fileName.isEmpty())
    return false;

  // load the markup
  return loadMarkupFrom(fileName.toStdString());
}

//: Process signal that an item in the image tree view has been double clicked
void nailfold_qfinal_gui::on_caseTreeView_clicked( const QModelIndex & index )
{
}

//: Process signal that an item in the image tree view has been double clicked
void nailfold_qfinal_gui::on_caseTreeView_doubleClicked( const QModelIndex & index )
{
	QModelIndex srcIndex = imageFilterModel_->mapToSource(index);
	QString fileName = baseFileModel_->filePath( srcIndex );
	vcl_cout << fileName.toStdString() << vcl_endl;
	if (!baseFileModel_->isDir(srcIndex))
	{
		loadImageFrom(fileName.toStdString());
	}	
}

void nailfold_qfinal_gui::on_actionSelectCase_triggered()
{
  selectCaseDialog();
}

void nailfold_qfinal_gui::on_actionLoadImage_triggered()
{
  loadImageDialog();
}

void nailfold_qfinal_gui::on_actionLoadMarkup_triggered()
{
  loadMarkupDialog();
}

//
//: Save markup to text file
bool nailfold_qfinal_gui::saveMarkupDialog()
{
  // Get a filename from the user via a dialog box

  // Propose a filename based on the image name
  const vcl_string initial_filename = 
      vul_file::strip_extension(imageProcessor_.filename()) + "_markup.txt";

  // get filename from dialog box
	QString selfilter;
	QString fileName = QFileDialog::getSaveFileName(
    this, "Save markup", initial_filename.c_str(),
    tr("Text files (*.txt)" ), // filter
    &selfilter );

  // if cancelled, do nothing (return)
  if (fileName.isEmpty())
    return false;

  // set markup filename
  markup_.set_filename(fileName.toStdString());

  return true;
}

//
//: Return error code based on success of saving operation:
//    0   Success
//    1   Filename not defined
//    2   Save to Inbox failed
//    3   Delete backup failed
//    4   Save copy to Outbox failed
//    5   Transfer Outbox to remote server failed
int nailfold_qfinal_gui::saveMarkup()
{
  bool filename_defined = !(markup_.filename().empty());

  // If no filename defined then get one from a dialog
  if (!filename_defined)
  {
    //filename_defined = saveMarkupDialog();

    markup_.set_filename(
        vul_file::strip_extension(imageProcessor_.filename()) + "_markup.txt");
    filename_defined = true;
  }

  if (!filename_defined)
    return 1;

  // Set the observer name before we try saving
  markup_.set_observer_name(preferences_.username());

  if (!markup_.save())
    return 2;

  if (!deleteBackup())
    return 3;

  // Flag this image as marked once it has been saved locally
  fileServer_.set_current_marked();

  // Copy file to Outbox folder.
  if (!fileServer_.copy_to_outbox(markup_.filename(), 
                                   markup_.tagged_filename()))
    return 4;

  // Empty Outbox to remote server.
  if (preferences_.use_remote_server() &&
      fileServer_.push_outbox_to(*serverHandler_) > 0)
    return 5;

  // Otherwise, return 0 for complete success.
  return 0;
}

void nailfold_qfinal_gui::on_actionSaveMarkup_triggered()
{
  saveMarkup();
  closeIfFinished();
}

//
//: Give the opportunity to save if the file has been modified.
//  Returns true if the markup is considered up to date - not necessarily 
//  because it wasn't saved. It is considered up to date if:
//  - It wasn't modified to begin with
//  - The user chose to save it (selected 'Yes' when prompted)
//  - The user chose to continue without saving (selected 'No' when prompted)
//
//  Only if the user selects 'Cancel' when prompted does the function return
//  false (i.e. the markup is not up to date).
bool nailfold_qfinal_gui::saveIfModified()
{
  // If not modified then do nothing
  if (!markup_.is_modified())
    return true;

  QMessageBox msgBox;
  msgBox.setText("Markup has been modified");
  msgBox.setInformativeText("Do you want to save this markup?");
  msgBox.setStandardButtons(QMessageBox::Yes | 
                            QMessageBox::No | 
                            QMessageBox::Cancel );
  msgBox.setDefaultButton(QMessageBox::Yes);

  // if user chooses not to save, or chooses to save and does so, return true
  const int response = msgBox.exec();
  if (response == QMessageBox::Yes)
  {
    // 'Yes'
    return (saveMarkup() == 0);
  }
  else if (response == QMessageBox::No)
  {
    // 'No'
    deleteBackup();
    // Reset modified flag to indicate that the user has defined this as up to
    // date.
    markup_.set_modified(false);
    return true;
  }
  else
  {
    // 'Cancel'
    return false;
  }
}


//
//: User selected 'Preferences' so bring up dialog box
void nailfold_qfinal_gui::on_actionPreferences_triggered()
{
  ncm_qfinal_options optionsWindow(preferences_, this);
  optionsWindow.exec();

  updateVesselList();

  // Update grid spacing
  const double gridSpacingPixels = 
    preferences_.grid_pixels_per_mm() * preferences_.grid_spacing_mm();
  scene_.set_grid_spacing(gridSpacingPixels);

  // Update UI elements that are visible
  applyUserLevelChange();
}

//
//: Enter "What's This" mode
void nailfold_qfinal_gui::on_actionWhatsThis_triggered()
{
  QWhatsThis::enterWhatsThisMode();
}

void nailfold_qfinal_gui::on_actionHelp_triggered()
{
QDesktopServices::openUrl(QUrl::fromLocalFile(QString::fromStdString(
    "ncm_software.chm")));
}

//
//: Brightness and contrast controls
void nailfold_qfinal_gui::on_brightnessSlider_sliderMoved(int value)
{
  const int contrast = ui.contrastSlider->value();
  const int brightness = value;
  const int lower_limit = (255-brightness) - (255-contrast)/2;
  const int upper_limit = (255-brightness) + (255-contrast)/2;

  imageProcessor_.set_contrast(lower_limit, upper_limit);
}

void nailfold_qfinal_gui::on_contrastSlider_sliderMoved(int value)
{
  const int contrast = value;
  const int brightness = ui.brightnessSlider->value();
  const int lower_limit = (255-brightness) - (255-contrast)/2;
  const int upper_limit = (255-brightness) + (255-contrast)/2;

  imageProcessor_.set_contrast(lower_limit, upper_limit);
}

//
//: Update line colormap when slider changed
void nailfold_qfinal_gui::on_compSizeSlider_sliderMoved(int value)
{
  imageProcessor_.set_line_size_threshold(value);
}

//
//: Now that we have stopped sliding, find the pixels that are nonzero
void nailfold_qfinal_gui::on_compSizeSlider_sliderReleased()
{
  imageProcessor_.get_line_pixels();
}

void nailfold_qfinal_gui::on_deleteAllVessels_clicked()
{
  // Prompt before doing anything
  QMessageBox msgBox;
  msgBox.setText("Delete All Vessels");
  msgBox.setInformativeText("Are you sure you want to delete all vessels?");
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  msgBox.setDefaultButton(QMessageBox::No);

  // delete all vessels from the scene (which, in turn, delete vessel entries
  // from the markup)
  if (msgBox.exec() == QMessageBox::Yes)
    scene_.delete_all_vesselitems();
}
void nailfold_qfinal_gui::on_deleteAllHaemorrhages_clicked()
{
  // Prompt before doing anything
  QMessageBox msgBox;
  msgBox.setText("Delete All Haemorrhages");
  msgBox.setInformativeText("Are you sure you want to delete all haemorrhages?");
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  msgBox.setDefaultButton(QMessageBox::No);

  // delete all haemorrhages from the scene (which, in turn, delete haemorrhage
  // entries from the markup)
  if (msgBox.exec() == QMessageBox::Yes)
    scene_.delete_all_haemorrhageitems();
}

void nailfold_qfinal_gui::setGrade(int grade_index)
{
  markup_.image_grade().set_value(grade_index);
}

void nailfold_qfinal_gui::updateGridSpacing()
{
  // Update grid spacing
  const double gridSpacingPixels = 
    preferences_.grid_pixels_per_mm() * preferences_.grid_spacing_mm();

  scene_.set_grid_spacing(gridSpacingPixels);
}

void nailfold_qfinal_gui::on_graphicsView_viewModeChanged(int)
{
  updateStageControls();
  updateToolbar();
  updateHelpPanel();
}

void nailfold_qfinal_gui::onImageProcessor_imageLoaded(bool success)
{
  // Keep scene and sceneview correctly sized
  scene_.fit_to_image();
  scene_.refresh_grid();
  ui.actionSaveMarkup->setEnabled(success);
  ui.graphicsView->fit_both();
  updateZoomLimits();
}

//
//  Private methods
//

void nailfold_qfinal_gui::createModeActionGroup()
{
  // not ideal, but seems to be the only way to do this
  modeGroup_ = new QActionGroup(this);
  modeGroup_->addAction(ui.actionPan);
  modeGroup_->addAction(ui.actionLabelImage);
  modeGroup_->addAction(ui.actionAddVessels);
  modeGroup_->addAction(ui.actionDefineVesselSize);
  modeGroup_->addAction(ui.actionDefineVesselShape);
  modeGroup_->addAction(ui.actionLabelApices);
  modeGroup_->addAction(ui.actionDrawVesselPath);
  modeGroup_->addAction(ui.actionAddHaemorrhages);
  modeGroup_->setExclusive(true);
  ui.actionLabelImage->setChecked(true);
}

void nailfold_qfinal_gui::createGradeRadioGroup()
{
  // Add grade radio buttons to a group for convenience
  gradeRadioGroup_ = new QButtonGroup(this);
  gradeRadioGroup_->addButton(ui.gradeNormalRadio, 1);
  gradeRadioGroup_->addButton(ui.gradeEarlyRadio, 2);
  gradeRadioGroup_->addButton(ui.gradeActiveRadio, 3);
  gradeRadioGroup_->addButton(ui.gradeLateRadio, 4);
  gradeRadioGroup_->addButton(ui.gradeExtremeRadio, 5);
  gradeRadioGroup_->addButton(ui.gradePoorQualityRadio, 6);
  gradeRadioGroup_->addButton(ui.gradeNonspecificRadio, 7);
  gradeRadioGroup_->setExclusive(true);
  //ui.gradeUndefinedRadio->setChecked(true);
}
void nailfold_qfinal_gui::createRowRadioGroup()
{
  // Add vessel row buttons to a group for convenience
  rowRadioGroup_ = new QButtonGroup(this);
  rowRadioGroup_->addButton(ui.distalRadio, 0);
  rowRadioGroup_->addButton(ui.nondistalRadio, 1);
  rowRadioGroup_->setExclusive(true);
  ui.distalRadio->setChecked(true);
}
void nailfold_qfinal_gui::createSizeRadioGroup()
{
  // Add vessel size radio buttons to a group for convenience
  sizeRadioGroup_ = new QButtonGroup(this);
  sizeRadioGroup_->addButton(ui.sizeNormalRadio, 0);
  sizeRadioGroup_->addButton(ui.sizeEnlargedRadio, 1);
  sizeRadioGroup_->addButton(ui.sizeGiantRadio, 2);
  sizeRadioGroup_->addButton(ui.sizeIrregularRadio, 3);
  sizeRadioGroup_->setExclusive(true);
  ui.sizeNormalRadio->setChecked(true);
}
void nailfold_qfinal_gui::createShapeRadioGroup()
{
  // Add vessel shape radio buttons to a group for convenience
  shapeRadioGroup_ = new QButtonGroup(this);
  shapeRadioGroup_->addButton(ui.shapeNormalRadio, 0);
  shapeRadioGroup_->addButton(ui.shapeTortuousRadio, 1);
  shapeRadioGroup_->addButton(ui.shapeRamifiedRadio, 2);
  shapeRadioGroup_->setExclusive(true);
  ui.shapeNormalRadio->setChecked(true);
}


//
//: Enable/disable controls that are dependent on user level (basic or advanced)
void nailfold_qfinal_gui::applyUserLevelChange()
{
  const bool advanced_user = preferences_.is_advanced_user();

  // Disable anything associated with drawing the vessel path.
  ui.actionDrawVesselPath->setEnabled(advanced_user);
  ui.drawVesselPathPage->setEnabled(advanced_user);

  // Change enabled status of controls in apex labelling mode.
  ui.show_edges_checkBox->setVisible(advanced_user);

  updateProgressFrame();
}

//
//: Update radio buttons to reflect current image grade.
void nailfold_qfinal_gui::updateGradeRadios()
{
  if (markup_.image_grade().value() == ncm_image_grade::GradeUndefined)
    uncheckGradeRadios();
  else
    gradeRadioGroup_->button(markup_.image_grade().value())->toggle();
}

//
//: Uncheck all of the radio buttons corresponding to image grade (i.e. if the
//  image is 'Undefined').
void nailfold_qfinal_gui::uncheckGradeRadios()
{
  // Make the group temporarily non-exclusive, because exclusive groups aren't 
  // allowed to be all unchecked.
  gradeRadioGroup_->setExclusive(false);

  for (int i = 0; i < gradeRadioGroup_->buttons().count(); ++i)
  {
    gradeRadioGroup_->buttons()[i]->setChecked(false);
  }

  gradeRadioGroup_->setExclusive(true);
}

//
//: Initialize the registry accessor and default parameters if they don't
//  already exist.
void nailfold_qfinal_gui::initialize_registry_accessor()
{
  ra_.get_current_user_key();

  long result = 
      ra_.get_subkey("Software\\University of Manchester\\NCM QCapture");

  if (result != ERROR_SUCCESS)
    ra_.create_subkey("Software\\University of Manchester\\NCM QCapture");
}

//
//: Get image and markup directory paths from the registry
void nailfold_qfinal_gui::get_registry_paths()
{
	long result = ERROR_SUCCESS;
  vcl_string dir_path;
	
	/*
	//Get image directory
	result = ra_.read_string("ncmImageDir", dir_path);

	if (result == ERROR_SUCCESS) {
		updateImagePath(dir_path);
	}
	else {
		// get dir path from dialog box
		QString qdir_path = QFileDialog::getExistingDirectory(
			this, tr("Select directory containing images"), currentImagePath_.c_str(),
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

		if (!qdir_path.isEmpty())
			updateImagePath(qdir_path.toStdString());		
	}
	*/
	/*
	//Get markup dir
	result = ra_.read_string("ncmMarkupDir", dir_path);	
	
	if (result == ERROR_SUCCESS) {
		updateMarkupPath(dir_path);
	}
	else {
		// get dir path from dialog box
		QString qdir_path = QFileDialog::getExistingDirectory(
			this, tr("Select directory containing image markup"), currentMarkupPath_.c_str(),
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

		if (!qdir_path.isEmpty())
			updateMarkupPath(qdir_path.toStdString());
	}
	*/
	
	//Get root dir
	result = ra_.read_string("ncmRootCaseDir", dir_path);	

	if (result == ERROR_SUCCESS) {
		updateRootPath(dir_path);
	}
	else {
		// get dir path from dialog box
		QString qdir_path = QFileDialog::getExistingDirectory(
			this, tr("Select directory containing all cases"), rootCasePath_.c_str(),
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

		if (!qdir_path.isEmpty())
			updateRootPath(qdir_path.toStdString());
	}
}

//
//: Does what it says on the tin
void nailfold_qfinal_gui::connectSignalsToSlots()
{
  QObject::connect( ui.autoContrast, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setAutoContrast(bool)) );
  QObject::connect( ui.autoContrast, SIGNAL(toggled(bool)),
                    this, SLOT(updateContrastControls()) );

  // make connections between UI widgets and local member variables
  // we could possibly make these slots part of the graphicsview rather than the
  // scene so that the connections can be set up in QDesigner (but it's not a
  // big deal at the moment).
  QObject::connect( ui.show_edges_checkBox, SIGNAL(toggled(bool)),
                    &scene_, SLOT(set_edges_visible(bool)) );
  QObject::connect( ui.show_lines_checkBox, SIGNAL(toggled(bool)),
                    &scene_, SLOT(set_lines_visible(bool)) );
  QObject::connect( ui.show_lines_checkBox, SIGNAL(toggled(bool)),
                    ui.snap_checkBox, SLOT(setEnabled(bool)) );
  QObject::connect( ui.snap_checkBox, SIGNAL(toggled(bool)),
                    &scene_, SLOT(set_snapping(bool)) );
  QObject::connect( ui.show_vessels_checkBox, SIGNAL(toggled(bool)),
                    &scene_, SLOT(set_vessels_visible(bool)) );

  QObject::connect( ui.show_lines_checkBox, SIGNAL(toggled(bool)),
                    &imageProcessor_, SLOT(find_ridges()) );
  QObject::connect( ui.show_edges_checkBox, SIGNAL(toggled(bool)),
                    &imageProcessor_, SLOT(find_edges()) );

  QObject::connect( ui.actionToggleImageControls, SIGNAL(toggled(bool)),
                    ui.imageGroup, SLOT(setVisible(bool)) );
  QObject::connect( ui.actionToggleZoomControls, SIGNAL(toggled(bool)),
                    ui.zoomGroup, SLOT(setVisible(bool)) );
  QObject::connect( ui.actionToggleGrid, SIGNAL(toggled(bool)),
                    &scene_, SLOT(set_grid_visible(bool)) );

  QObject::connect( &imageProcessor_, SIGNAL(imageLoaded(bool)),
                    this, SLOT(onImageProcessor_imageLoaded(bool)) );

  // When vessel combobox changes selection, select that vessel and zoom in
  QObject::connect( &scene_, SIGNAL(selectionChanged(int)),
                    ui.vesselCombo, SLOT(setCurrentIndex(int)) );
  QObject::connect( &scene_, SIGNAL(selectionChanged(int)),
                    this, SLOT(updateVesselControls()) );
  QObject::connect( &scene_, SIGNAL(selectionChanged(int)),
                    this, SLOT(applyVesselZoom(int)) );

  // When a vessel is picked from the list, make it selected in the scene, zoom
  // in on it, and update the vessel controls
  QObject::connect( ui.vesselCombo, SIGNAL(currentIndexChanged(int)),
                    &scene_, SLOT(select_vessel(int)) );
  // Treat any click, regardless of whether it changes the selection, as a
  // reason to redraw and select
  QObject::connect( ui.vesselCombo, SIGNAL(activated(int)),
                    &scene_, SLOT(select_vessel(int)) );
  // This just works nicely - when sweeping over the list, the screen updates
  // in realtime.
  QObject::connect( ui.vesselCombo, SIGNAL(highlighted(int)),
                    &scene_, SLOT(select_vessel(int)) );

  // React to changes in scene
  QObject::connect( &scene_, SIGNAL(vesselAdded()),
                    this, SLOT(updateVesselList()) );
  QObject::connect( &scene_, SIGNAL(vesselDeleted()),
                    this, SLOT(updateVesselList()) );
  QObject::connect( &scene_, SIGNAL(vesselChanged()),
                    this, SLOT(updateVesselList()) );

  QObject::connect( ui.graphicsView, SIGNAL(changed()),
                    this, SLOT(updateZoomControls()) );
  QObject::connect( ui.graphicsView, SIGNAL(previous()),
                    this, SLOT(on_prevVesselButton_clicked()) );
  QObject::connect( ui.graphicsView, SIGNAL(next()),
                    this, SLOT(on_nextVesselButton_clicked()) );

  // Changing the image grade
  QObject::connect( gradeRadioGroup_, SIGNAL(buttonClicked(int)),
                    this, SLOT(setGrade(int)) );

  // Set vessel properties
  QObject::connect( rowRadioGroup_, SIGNAL(buttonClicked(int)),
                    this, SLOT(updateRowButtons(int)) );
  QObject::connect( ui.distalRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setDistal(bool)) );
  QObject::connect( ui.nondistalRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setNondistal(bool)) );

  QObject::connect( sizeRadioGroup_, SIGNAL(buttonClicked(int)),
                    this, SLOT(updateSizeButtons(int)) );
  QObject::connect( ui.sizeNormalRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setSizeNormal(bool)) );
  QObject::connect( ui.sizeEnlargedRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setSizeEnlarged(bool)) );
  QObject::connect( ui.sizeGiantRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setSizeGiant(bool)) );
  QObject::connect( ui.sizeIrregularRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setSizeIrregular(bool)) );
  
  QObject::connect( shapeRadioGroup_, SIGNAL(buttonClicked(int)),
                    this, SLOT(updateShapeButtons(int)) );
  QObject::connect( ui.shapeNormalRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setShapeNormal(bool)) );
  QObject::connect( ui.shapeTortuousRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setShapeTortuousOnly(bool)) );
  QObject::connect( ui.shapeRamifiedRadio, SIGNAL(toggled(bool)),
                    &scene_.new_vessel_properties(), SLOT(setShapeRamifiedOnly(bool)) );

  QObject::connect( &imageProcessor_, SIGNAL(rawImageChanged()),
                    ui.graphicsView, SLOT(repaint()) );
  QObject::connect( &imageProcessor_, SIGNAL(rawImageChanged()),
                    this, SLOT(updateContrastControls()) );
  QObject::connect( &imageProcessor_, SIGNAL(imageTypeChanged()),
                    this, SLOT(updateToolbar()) );
  QObject::connect( &imageProcessor_, SIGNAL(imageTypeChanged()),
                    this, SLOT(updateStageControls()) );

  // Signals that should update the status bar
  QObject::connect( &imageProcessor_, SIGNAL(statusChanged(const QString&, double)),
                    this, SLOT(updateStatus(const QString&, double)) );
  QObject::connect( &scene_, SIGNAL(apexLengthChanged(double)),
                    this, SLOT(onApexLengthChanged(double)) );

  // Zoom controls
  QObject::connect( ui.zoomFitButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(fit_both()));
  QObject::connect( ui.zoomFitWButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(fit_width()));
  QObject::connect( ui.zoomFitHButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(fit_height()));

  connectSetFocus();
}

//
//: Connect signals from all controls such that they return focus to the
//  graphicsView after every operation.
void nailfold_qfinal_gui::connectSetFocus()
{
  // After any of the following events, return the focus to the
  // graphics view to process keyboard events
  QObject::connect( ui.brightnessSlider, SIGNAL(sliderReleased()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.contrastSlider, SIGNAL(sliderReleased()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.autoContrast, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setFocus()));

  QObject::connect( ui.zoomSlider, SIGNAL(sliderReleased()),
                    ui.graphicsView, SLOT(setFocus()));
  //QObject::connect( ui.zoomEdit, SIGNAL(editFinished()),
  //                  ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.zoomInButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.zoomOutButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.zoomFitButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.zoomFitWButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.zoomFitHButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));

  QObject::connect( ui.vesselCombo, SIGNAL(currentIndexChanged(int)),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.firstVesselButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.prevVesselButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.nextVesselButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.lastVesselButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.deleteAllVessels, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.deleteAllHaemorrhages, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));

  QObject::connect( ui.stackedWidget, SIGNAL(currentChanged(int)),
                    ui.graphicsView, SLOT(setFocus()));

  QObject::connect( ui.sizeNormalRadio, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.sizeEnlargedRadio, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.sizeGiantRadio, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setFocus()));
  
  QObject::connect( ui.shapeNormalRadio, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.shapeTortuousRadio, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.shapeRamifiedRadio, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setFocus()));

  QObject::connect( ui.show_edges_checkBox, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setFocus()));

  QObject::connect( ui.show_lines_checkBox, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.compSizeSlider, SIGNAL(sliderReleased()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.snap_checkBox, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.show_vessels_checkBox, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setFocus()));

  QObject::connect( ui.prevStageButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));
  QObject::connect( ui.nextStageButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(setFocus()));
}


