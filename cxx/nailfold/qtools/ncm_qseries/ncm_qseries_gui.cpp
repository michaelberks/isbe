#include "ncm_qseries_gui.h"

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

#include "ncm_qseries_options.h"

//
//: Constructor
nailfold_qseries_gui::nailfold_qseries_gui(ncm_qseries_preferences& preferences,
                                           QWidget *parent, Qt::WFlags flags)
: QMainWindow(parent, flags),
  preferences_(preferences),
  constructed_(false),
  scene_(this),
	scene2_(this),
	currentStage_(Stage1)
{
  // Name objects to allow autoconnection of signals to slots
  scene_.setObjectName("scene");
	scene2_.setObjectName("scene");
  imageProcessor_.setObjectName("imageProcessor");
	imageProcessor2_.setObjectName("imageProcessor");

  // setup the UI
  ui.setupUi(this);

  setWindowTitle(windowTitle());

  style_ = new QWindowsXPStyle;

  // Get standard icons from the platform-dependent theme
  setIcons();

  // associate the scene with the markup and image processor
	ui.gridButton->setDown(true);
  scene_.set_image_processor(&imageProcessor_);
	scene_.set_grid_visible(ui.gridButton->isDown());

	scene2_.set_image_processor(&imageProcessor2_);
	scene2_.set_grid_visible(ui.gridButton->isDown());

  const double border = preferences_.normal_vessel_zoom() *
                        preferences_.enlarged_zoom_relative() * 
                        preferences_.giant_zoom_relative();
  scene_.set_border_size(border);
	scene2_.set_border_size(border);

  // use Pan mode by default, with reverse zoom (wheel up = zoom in)
	ui.graphicsView->usePanMode();
  ui.graphicsView->reverse_zoom();
  ui.graphicsView->setFocus();

	ui.graphicsViewBottom->usePanMode();
  ui.graphicsViewBottom->reverse_zoom();
  ui.graphicsViewBottom->setFocus();

  // tell graphics view to use this scene
  ui.graphicsView->setScene(&scene_);
	ui.graphicsViewBottom->setScene(&scene2_);

  // Make first item selected on startup
  ui.stackedWidget->setCurrentIndex(0);

	//Hide the help to start with
	ui.tabControl->hide();
  
  // Create and initialize various objects
  createGradeRadioGroup();
  createReasonsRadioGroup();

	ui.stackedWidget->setEnabled( false );

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

/*#else
  // hide the Debug button in the menubar for the Release version
  ui.actionDebug->setVisible(false);
  ui.actionCrash->setVisible(false);
#endif*/ 
  
  // Flag to show that constructor ended successfully
  constructed_ = true;
}

//
//: Destructor
nailfold_qseries_gui::~nailfold_qseries_gui()
{

  imageThread_.quit();

  // Stop taking snapshots now that we're finished.
  // Otherwise, when the scene cleans itself up by deleting the ncm_vesselitems,
  // the subsequent deletion of ncm_vessels causes a new snapshot to be taken
  // such that the backup file is created again.
  grade_.set_snapshots(false);

  // If backup file exists then delete it
  vcl_string backup_filename = grade_.filename() + '~';
  QFile(backup_filename.c_str()).remove();
}

bool nailfold_qseries_gui::initialize()
{

	//The images should all be local on the machine so just get them from the image list
	getSiteName();
	imageManager_.read_image_list();
	statusBar_main_.setText("");

  // Go to first image in Inbox and load it
  bool first_found = imageManager_.first();
  if (first_found)
    loadImage();
  else
    hideScene();

  return true;
}

QString nailfold_qseries_gui::windowTitle() const
{
  QString version_string = QString::number(preferences_.version());
  return QString("NCM Markup v" + version_string);
}

bool nailfold_qseries_gui::wasConstructed() const
{
  return constructed_;
}

void nailfold_qseries_gui::hideScene()
{
  ui.graphicsView->setScene(NULL);
	ui.graphicsViewBottom->setScene(NULL);

  ui.imageGroup->setEnabled(false);
  ui.zoomGroup->setEnabled(false);

  ui.stackedWidget->setEnabled(false);
}

void nailfold_qseries_gui::showScene()
{
  ui.graphicsView->setScene(&scene_);
  ui.graphicsView->fit_both();

	ui.graphicsViewBottom->setScene(&scene2_);
  ui.graphicsViewBottom->fit_both();

  ui.imageGroup->setEnabled(true);
  ui.zoomGroup->setEnabled(true);

  updateStageControls();
}

bool nailfold_qseries_gui::getSiteName()
{
  int error_code = -1;

  error_code = imageManager_.read_site_name();

  if (error_code == 0)
    return true;

  if (error_code == 1)
  {
    // File not found - create a new one using the site selector
    ncm_qsite_selector siteSelector(this);
    siteSelector.exec();

    // Try once more - if it worked, great
    error_code = imageManager_.read_site_name();
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


//
//: Optionally close the application if there are no more images to mark up
bool nailfold_qseries_gui::closeIfFinished()
{
  // If that was the last image then ask user if they want to close
  if (imageManager_.is_empty())
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

void nailfold_qseries_gui::setIcons()
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

void nailfold_qseries_gui::setupImageThread()
{
  // Assign all image processing to a separate thread to avoid the GUI hanging
  // up on large images
  imageProcessor_.moveToThread(&imageThread_);
	imageProcessor2_.moveToThread(&imageThread_);

  // Have the thread running at all times to handle brightness and contrast
  // changes
  imageThread_.start();
}

//
//: Create the status bar's labels and progress bar
void nailfold_qseries_gui::setupStatusbar()
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
void nailfold_qseries_gui::updateStatus(const QString& status, 
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

void nailfold_qseries_gui::updateStatusLeft()
{
  QString statusString = QString::number(imageManager_.n_unmarked()) + 
                         " images left to mark";
  statusBar_LPanel_.setText(statusString);
}

//
//: Apply a specified image zooming to the graphics view
void nailfold_qseries_gui::applyImageZoom(ncm_qseries_preferences::ImageZoom zoom)
{
  switch (zoom)
  {
    case ncm_qseries_preferences::ZoomNoResize:
      // do nothing
      break;
    case ncm_qseries_preferences::ZoomFit:
      ui.graphicsView->fit_both();
			ui.graphicsViewBottom->fit_both();
      break;
    case ncm_qseries_preferences::ZoomFitWidth:
      ui.graphicsView->fit_width();
			ui.graphicsViewBottom->fit_width();
      break;
    case ncm_qseries_preferences::ZoomFitHeight:
      ui.graphicsView->fit_height();
			ui.graphicsViewBottom->fit_height();
      break;
    default:
      // throw error?
      break;
  }
}

void nailfold_qseries_gui::applyUserZoom()
{
	applyImageZoom(preferences_.label_image_zoom());
}

//
//: Update toolbar at top of screen
void nailfold_qseries_gui::updateToolbar()
{
  // Enable/disable buttons depending on whether the image is a dermatogram or
  // capillarogram
  const bool image_is_capillarogram = 
      (imageProcessor_.image_type() == ncm_qimagehandler::TypeCapillarogram);

  const bool user_is_advanced = preferences_.is_advanced_user();
}

//
//: Update brightness and contrast that may change as a result of autocontrast
void nailfold_qseries_gui::updateContrastControls()
{
  const int min_contrast = imageProcessor_.min_contrast();
  const int max_contrast = imageProcessor_.max_contrast();
  const int brightness = 255 - (min_contrast + max_contrast)/2;
  const int contrast = 255 - max_contrast + min_contrast;

  ui.brightnessSlider->setValue(brightness);
  ui.contrastSlider->setValue(contrast);

  ui.brightnessSlider->setEnabled(!ui.autoContrast->isChecked());
  ui.contrastSlider->setEnabled(!ui.autoContrast->isChecked());

	const int min_contrast2 = imageProcessor2_.min_contrast();
  const int max_contrast2 = imageProcessor2_.max_contrast();
  const int brightness2 = 255 - (min_contrast2 + max_contrast2)/2;
  const int contrast2 = 255 - max_contrast2 + min_contrast2;

  ui.brightnessSlider_2->setValue(brightness2);
  ui.contrastSlider_2->setValue(contrast2);

  ui.brightnessSlider_2->setEnabled(!ui.autoContrast_2->isChecked());
  ui.contrastSlider_2->setEnabled(!ui.autoContrast_2->isChecked());
}
 
//: Update vessel selection controls, depending on how many items are in the
//  list and which one is selected
//void nailfold_qseries_gui::updateVesselControls()


//
//: Update limits on zoom controls
void nailfold_qseries_gui::updateZoomLimits()
{
  const double minimum_scale = ui.graphicsView->fitHeightScale();
  ui.zoomSlider->setMinimum(static_cast<int>(100*minimum_scale + 1));

  const double maximum_scale = 10 * minimum_scale;
  ui.zoomSlider->setMaximum(static_cast<int>(100*maximum_scale - 1));

  ui.graphicsView->setMaxScale(maximum_scale);
}

//
//: Update zoom controls
void nailfold_qseries_gui::updateZoomControls()
{
  const double scale = ui.graphicsView->transform().m11();
  ui.zoomSlider->setValue(static_cast<int>(100*scale + 1));

  QString zoom_string = QString::number(static_cast<int>(100*scale + 1));
  ui.zoomEdit->setText(zoom_string+"%");
}

//
//: Make sure help panel is up to date
void nailfold_qseries_gui::updateHelpPanel()
{
  // Display correct help page for this mode
  /*if (ui.graphicsView->isPanning())
  {
    ui.helpStackedWidget->setCurrentWidget(ui.panHelpPage);
    return;
  }*/

  switch (currentStage_)
  {
  case Stage1:
		ui.helpStackedWidget->setCurrentWidget(ui.stage1HelpPage);
		break;
  case Stage2:
    ui.helpStackedWidget->setCurrentWidget(ui.stage2HelpPage);
    break;

  default:
    assert(false);
  }
}

//
//  Private slots
//

//
//: Update stage controls
void nailfold_qseries_gui::updateStageControls()
{
  // Display correct set of controls for this mode
  ui.stackedWidget->setCurrentIndex(currentStage_);

  // If there are no images then disable 
  // widget altogether and return
  if (imageManager_.is_empty())
  {
    ui.stackedWidget->setEnabled(false);
		return;
  }
	
	// Otherwise, stacked widget should be enabled
  ui.stackedWidget->setEnabled(true);


  //If we are in mode one 
	if (currentStage_ == Stage1) 
	{
		//If this is the first image disable "Go Back"
		if (imageManager_.is_first()) //MB isfirst image
			ui.goBackButton->setEnabled(false);
		else
			ui.goBackButton->setEnabled(true);

		//If a normal/abnormal grade hasn't been given disable "Continue"
		if (!grade_.is_graded())
		{
			//Disable the slider controls
			ui.progressionSlider->setValue(ncm_qseries_grade::None);
			ui.progressionLabel->setEnabled(false);
			ui.progressionBox->setEnabled(false);
			ui.continueButton->setEnabled(false);
			
		}
		else
		{
			//We can continue as we have a grade
			ui.continueButton->setEnabled(true);

			if (grade_.is_normal())
			{
				//Don't mark the progression slider
				ui.progressionSlider->setValue(ncm_qseries_grade::None);
				ui.progressionLabel->setEnabled(false);
				ui.progressionBox->setEnabled(false);
				
			}
			else //At least one image is abnormal, enable the progression slider
			{
				ui.progressionLabel->setEnabled(true);
				ui.progressionBox->setEnabled(true);
				ui.progressionSlider->setValue(grade_.progression_level());
			}
		}
	}
	else //Must be in mode 2 (currentStage_ == Stage1)
  {
		//Set text appropriately
		switch (grade_.progression_level())
		{
		case ncm_qseries_grade::MajorTop:
			{
				ui.reasonsStatus->setText("You selected the top image \nas showing significant signs of \ndisease progression");
				break;
			}
		case ncm_qseries_grade::MinorTop:
			{
				ui.reasonsStatus->setText("You selected the top image \nas showing some signs of \ndisease progression");
				break;
			}
		case ncm_qseries_grade::MinorBottom:
			{
				ui.reasonsStatus->setText("You selected the bottom image \nas showing some signs of \ndisease progression");
				break;
			}
		case ncm_qseries_grade::MajorBottom:
			{
				ui.reasonsStatus->setText("You selected the bottom image \nas showing significant signs of \ndisease progression");
				break;
			}
		case ncm_qseries_grade::None:
			{
				//Something has gone wrong!
				assert(false);
			}

			default:
				assert(false);
				//Something has gone wrong!
		}

		//Both go back and go forward enabled
		ui.goBackButton->setEnabled(true);
		ui.continueButton->setEnabled(true);
  }

}

void nailfold_qseries_gui::on_goBackButton_clicked()
{
  const int first_stage = 0;

	//If first stage go back to previous image
  if (ui.stackedWidget->currentIndex() == first_stage)
  {
		requestPreviousImage();
	}
	else
	{
		//Go back to first stage
    requestStageChange(Stage1);
  }
}
void nailfold_qseries_gui::on_continueButton_clicked()
{
		//If first stage
	if (ui.stackedWidget->currentIndex() == Stage1)
	{
		//If images normal, continue to next image
		if (grade_.is_graded())
		{
			if (grade_.is_normal() || (grade_.progression_level()==ncm_qseries_grade::None))
				requestNextImage();
			else //Images are abnormal, got to stage 2
				requestStageChange(Stage2);
		}
		//else: Do nothing, we shouldn't ever get here if continue has been enabled correctly
	}
	else //In second stage, continue to next image
	{
		requestNextImage();		
	}
}

void nailfold_qseries_gui::requestPreviousImage()
{
  saveIfModified();
  bool found_previous = imageManager_.previous();
  if (found_previous)
    loadImage();
		
  else
  {
    //We should never get to this, but just in case
    ui.goBackButton->setEnabled(false);
  }
}
void nailfold_qseries_gui::requestNextImage()
{
  // Scene edit mode is changed (without warnings) if the image is
  // successfully loaded. Apply warnings here, before the image is loaded.
  const bool dismiss_warning = warnOnStageChange();
  if (!dismiss_warning)
    return;

  if (!imageManager_.is_last())
  {
    // There are more images to come, so save data and continue
		saveIfModified();
    bool found_next = imageManager_.next();
    loadImage();
  }
  else
  {
    // This was the last image. Tell the user well done and ask if they want to close
		saveIfModified();
		closeIfFinished();
  }
}

void nailfold_qseries_gui::on_unmarkedCheckbox_toggled(bool checked)
{
  imageManager_.use_unmarked_only(checked);

	if (imageManager_.current_image().size() != 2)
  {
    // If no image was loaded then load the first one.
    // This can happen when all images have been marked then we switch to show
    // all images.
    imageManager_.first();
    loadImage();
  }
  else if (imageManager_.is_using_unmarked_only() &&
           imageManager_.is_current_marked())
  {
    // If we switch to unmarked images and the current image is marked then
    // look for an unmarked image to show. If there are no such images then
    // clear the image altogether.
    if (imageManager_.is_empty())
    {
      updateStageControls();
    }
    else
    {
      // Look forward.
      bool found = imageManager_.next();

      // If nothing there then look back.
      if (!found)
        imageManager_.previous();

      // Load the image found (there must be one as we've already checked that
      // the image queue is not empty
      loadImage();
    }
  }
  else
  {
    // Just update the stage and image controls
    updateStageControls();
  }
}

//
//: Toolbar actions

void nailfold_qseries_gui::on_gridButton_toggled()
{

	if (ui.gridButton->isChecked())
	{
		scene_.set_grid_visible(true);
		scene2_.set_grid_visible(true);
	}
	else
	{
		scene_.set_grid_visible(false);
		scene2_.set_grid_visible(false); 
	}
}

//
//  Zoom controls
//
void nailfold_qseries_gui::setZoomFrom(int value)
{
  const double minimum_scale_top = ui.graphicsView->fitScale();
	const double minimum_scale_bottom = ui.graphicsViewBottom->fitScale();
  double new_scale = static_cast<double>(value-1) / 100;

  if (new_scale < minimum_scale_top)
    new_scale = minimum_scale_top;
	if (new_scale < minimum_scale_bottom)
    new_scale = minimum_scale_bottom;

  ui.graphicsView->setScale(new_scale);
	ui.graphicsViewBottom->setScale(new_scale);
}
void nailfold_qseries_gui::on_zoomSlider_sliderMoved(int value)
{
  setZoomFrom(value);
}

void nailfold_qseries_gui::on_zoomEdit_editingFinished()
{
  // Remove last character (%) from string and convert
  QString value_string = ui.zoomEdit->text();
  value_string.chop(1);
  const int value = value_string.toInt();

  setZoomFrom(value);
}

void nailfold_qseries_gui::on_zoomOutButton_clicked()
{
  // Avoid calling functions that generate more signals, such as this:
  //ui.zoomSlider->setValue(ui.zoomSlider->value()-ui.zoomSlider->pageStep());

  // Instead, set the graphicsView zoom directly, and let the resulting signal
  // trigger an update of the GUI controls
  int new_zoom = ui.zoomSlider->value() - ui.zoomSlider->pageStep();
  setZoomFrom(new_zoom);
}
void nailfold_qseries_gui::on_zoomInButton_clicked()
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
void nailfold_qseries_gui::on_hidePanelButton_clicked()
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
//: Sound warning when no reasons have been given in edit mode 2
//  Returns true unless user elects to heed the warning.
bool nailfold_qseries_gui::warnOnStageChange()
{
	if (currentStage_ == Stage1)
		return true;
	
	//We must be in stage 2
	if (preferences_.warn_if_no_reasons)
  {
    // Warn if no reasons have been given
		// Alias for convenience
		QButtonGroup*& radioGroup = reasonsRadioGroup_;

		// selected should never be zero, as the button should be disabled in this
		// instance
		int selected = radioGroup->checkedId();
    if (selected == -1)
    {
			QMessageBox msgBox;
			msgBox.setIcon(QMessageBox::Warning);
			msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
			msgBox.setDefaultButton(QMessageBox::No);
			msgBox.setText("No reasons given");
			msgBox.setInformativeText("Continue?");
      
      return (msgBox.exec() == QMessageBox::Yes);
    }
  }
	return true;
}

//
//: Process a request to change edit mode.
void nailfold_qseries_gui::requestStageChange(int newStage)
{
  // If the new mode is the same as the current one then do nothing.
  if (static_cast<GradeStage>(newStage) == currentStage_)
    return;

  const bool dismiss_warning = warnOnStageChange();
  if (dismiss_warning)
    goToStage(newStage);
}

//
//: Move to the first stage of grading the image pairs
void nailfold_qseries_gui::goToStage(int stage)
{
	currentStage_ = static_cast<GradeStage>(stage);
	ui.stackedWidget->setCurrentIndex(stage);
	if (currentStage_ == Stage1)
		updateGradeRadios();
	else
		updateReasonsRadios();
	updateStageControls();
}

//
// Events
//
void nailfold_qseries_gui::closeEvent(QCloseEvent *ev)
{
  // only close if (a) markup is not modified or (b) user declines to save or 
  // (c) user successfully saves to a file.
  if (saveIfModified())
    ev->accept();
  else
    ev->ignore();
}

//
//  File I/O
//

//
//: Delete the temporary backup file
bool nailfold_qseries_gui::deleteBackup() const
{
  vcl_string backup_filename = grade_.snapshot_filename();
  return QFile(backup_filename.c_str()).remove();
}

//
//: Load an image from file
bool nailfold_qseries_gui::loadImage()
{

	vcl_vector<vcl_string> im_filenames = imageManager_.current_image();

  // Ignore if no image name is supplied
	if (im_filenames.size() != 2)
    return false;

  // Update the window title
  setWindowTitle(windowTitle() + " - " + im_filenames[0].c_str());

	// Clear scene (including markup) before we load anything
  // Turn off snapshots temporarily otherwise we'll create a new backup file
  grade_.set_snapshots(false);
  scene_.clear();
	scene2_.clear();
  grade_.set_snapshots(true);

  // Tell the image processor which file to load then load it.
  imageProcessor_.set_filename(im_filenames[0]);
  imageProcessor_.load_image();

	imageProcessor2_.set_filename(im_filenames[1]);
  imageProcessor2_.load_image();

  imageManager_.set_current(im_filenames[0]);

	vcl_cout << "Image " << im_filenames[0] << " loaded" << vcl_endl;

	vcl_string grade_filename = imageManager_.current_grade();
	const bool grade_exists = imageManager_.is_current_marked();

  // Define backup filename check if it exists
  const vcl_string backup_filename = grade_filename + '~';
  const bool backup_exists = vul_file::exists(backup_filename);

	//Start with a clear grade
	grade_.clear();
	grade_.set_filename(grade_filename);

  bool grade_loaded = false;

  // Check for an autosaved markup that still exists (indicative of a crash)
  if (backup_exists)
  {
    // if markup exists then compare timestamps to check if backup is newer
    bool backup_is_newer = false;
    if (grade_exists)
    {
      QFileInfo fileinfo;

      fileinfo.setFile(QString(backup_filename.c_str()));
      QDateTime backup_timestamp = fileinfo.lastModified();
      fileinfo.setFile(QString(grade_filename.c_str()));
      QDateTime grade_timestamp = fileinfo.lastModified();

      backup_is_newer = (grade_timestamp < backup_timestamp);
    }

    // If backup is there but not the 'non-backup', or if backup is newer than
    // the 'non-backup' then ask user if they want to restore from the backup
    // instead
    if (!grade_exists || backup_is_newer)
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
        grade_.t_read(backup_filename);
        saveMarkup(); // replace old file with restored backup
        grade_loaded = true;
      }
    }
  }

  if (!grade_loaded && grade_exists)
  {
    // Load markup and apply it to the scene
    grade_.t_read(grade_filename);
    grade_loaded = true;
  }		

  // Restart the timer
  grade_.start_timing();

  // Start new image in 'Grade Image' mode
  goToStage(0);
  updateStatusLeft();
  
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
int nailfold_qseries_gui::saveMarkup()
{
  //Check filename is defined
  if (grade_.filename().empty())
    return 1;

  // Set the observer name before we try saving
  grade_.set_observer_name(preferences_.username());

	//Save the grade here, should return 0 if ok
  if (!grade_.save())
    return 2;

	//Try deleting backup
  if (!deleteBackup())
    return 3;

  // Flag this image as marked once it has been saved locally
  imageManager_.set_current_marked();

  // Otherwise, return 0 for complete success.
  return 0;
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
bool nailfold_qseries_gui::saveIfModified()
{
  // If not modified then do nothing
  if (!grade_.is_modified())
    return true;

	//Otherwise try and save
	return (saveMarkup() == 0);

  /*QMessageBox msgBox;
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
    grade_.set_modified(false);
    return true;
  }
  else
  {
    // 'Cancel'
    return false;
  }
	*/
}

//
// User clicked to save markup
void nailfold_qseries_gui::on_actionSaveMarkup_triggered()
{
  saveMarkup();
  closeIfFinished();
}

//
//: User selected 'Preferences' so bring up dialog box
void nailfold_qseries_gui::on_actionPreferences_triggered()
{
  ncm_qseries_options optionsWindow(preferences_, this);
  optionsWindow.exec();

  // Update grid spacing
  const double gridSpacingPixels = 
    preferences_.grid_pixels_per_mm() * preferences_.grid_spacing_mm();
  scene_.set_grid_spacing(gridSpacingPixels);
	scene2_.set_grid_spacing(gridSpacingPixels);

  // Update UI elements that are visible
  applyUserLevelChange();
}

//
//: Brightness and contrast controls
void nailfold_qseries_gui::on_brightnessSlider_sliderMoved(int value)
{
  const int contrast = ui.contrastSlider->value();
  const int brightness = value;
  const int lower_limit = (255-brightness) - (255-contrast)/2;
  const int upper_limit = (255-brightness) + (255-contrast)/2;

  imageProcessor_.set_contrast(lower_limit, upper_limit);
}

void nailfold_qseries_gui::on_brightnessSlider_2_sliderMoved(int value)
{
  const int contrast = ui.contrastSlider_2->value();
  const int brightness = value;
  const int lower_limit = (255-brightness) - (255-contrast)/2;
  const int upper_limit = (255-brightness) + (255-contrast)/2;

	imageProcessor2_.set_contrast(lower_limit, upper_limit);
}

void nailfold_qseries_gui::on_contrastSlider_sliderMoved(int value)
{
  const int contrast = value;
  const int brightness = ui.brightnessSlider->value();
  const int lower_limit = (255-brightness) - (255-contrast)/2;
  const int upper_limit = (255-brightness) + (255-contrast)/2;

  imageProcessor_.set_contrast(lower_limit, upper_limit);
}

void nailfold_qseries_gui::on_contrastSlider_2_sliderMoved(int value)
{
  const int contrast = value;
  const int brightness = ui.brightnessSlider_2->value();
  const int lower_limit = (255-brightness) - (255-contrast)/2;
  const int upper_limit = (255-brightness) + (255-contrast)/2;

	imageProcessor2_.set_contrast(lower_limit, upper_limit);
}

//
//
// Are the images normal or abnormal
/*void nailfold_qseries_gui::on_imagesNormalRadio_toggled(bool checked)
{
	if (checked)
	{
		//Disable the slider controls
		ui.progressionLabel->setEnabled(false);
		ui.progressionBox->setEnabled(false);

		//Enable continue
		ui.continueButton->setEnabled(true);
	}

	//Set the grade
	grade_.set_normal(checked);

	//Make sure group is now exclusive
	gradeRadioGroup_->setExclusive(true);

	updateStageControls();
}

void nailfold_qseries_gui::on_imagesAbnormalRadio_toggled(bool checked)
{
	if (checked)
	{
		//Disable the slider controls
		ui.progressionLabel->setEnabled(true);
		ui.progressionBox->setEnabled(true);

		//Enable continue
		ui.continueButton->setEnabled(true);
	}

	//Set the grade
	grade_.set_normal(!checked);

	//Make sure group is now exclusive
	gradeRadioGroup_->setExclusive(true);

	updateStageControls();
}*/

//
//Choosing selecting reasons
void nailfold_qseries_gui::setGrade(int imageGrade)
{

	switch (static_cast<ncm_qseries_grade::ImageGrade>(imageGrade))
	{
	case ncm_qseries_grade::Normal:
		{
			grade_.set_normal(true);
			break;
		}
	case ncm_qseries_grade::Abnormal:
		{
			grade_.set_normal(false);
			break;
		}
	case ncm_qseries_grade::Ungradeable:
		{
			grade_.set_normal(true);
			break;
		}
		default:
			assert(false);
	}

	updateStageControls();
}

//
//Choosing selecting reasons
void nailfold_qseries_gui::selectReasons(int reason_idx)
{
	bool selected = reasonsRadioGroup_->button(reason_idx)->isChecked();
	grade_.set_reason(selected, static_cast<ncm_qseries_grade::ProgressReason>(reason_idx));
}

//
//Moving the progression slider
void nailfold_qseries_gui::on_progressionSlider_sliderMoved(int value)
{
	grade_.set_progression(static_cast<ncm_qseries_grade::ProgressionLevel>(value));
}


void nailfold_qseries_gui::updateGridSpacing()
{
  // Update grid spacing
  const double gridSpacingPixels = 
    preferences_.grid_pixels_per_mm() * preferences_.grid_spacing_mm();

  scene_.set_grid_spacing(gridSpacingPixels);
	scene2_.set_grid_spacing(gridSpacingPixels);
}

void nailfold_qseries_gui::on_graphicsView_viewModeChanged(int)
{
  updateStageControls();
  updateToolbar();
  updateHelpPanel();
}

void nailfold_qseries_gui::onImageProcessor_imageLoaded(bool success)
{
  // Keep scene and sceneview correctly sized
  scene_.fit_to_image();
  scene_.refresh_grid();

  ui.graphicsView->fit_both();
  updateZoomLimits();
}

void nailfold_qseries_gui::onImageProcessor2_imageLoaded(bool success)
{
  // Keep scene and sceneview correctly sized
	scene2_.fit_to_image();
  scene2_.refresh_grid();

  ui.graphicsViewBottom->fit_both();

	//After loading the bottom image, make sure scales are equal
	QTransform tt = ui.graphicsView->transform();
	QTransform tb = ui.graphicsViewBottom->transform();
	if (tt.m11() > tb.m11())
		ui.graphicsView->setScale(tb.m11());
	else
		ui.graphicsViewBottom->setScale(tt.m11());

  updateZoomLimits();
}

//
//  Private methods
//

void nailfold_qseries_gui::createGradeRadioGroup()
{
  // Add grade radio buttons to a group for convenience
  gradeRadioGroup_ = new QButtonGroup(this);
	gradeRadioGroup_->addButton(ui.imagesNormalRadio, ncm_qseries_grade::Normal);
  gradeRadioGroup_->addButton(ui.imagesAbnormalRadio, ncm_qseries_grade::Abnormal);
  ui.imagesNormalRadio->setChecked(false);
	ui.imagesAbnormalRadio->setChecked(false);
	ui.imagesUngradeableRadio->setChecked(false);
	gradeRadioGroup_->setExclusive(false);
	//gradeRadioGroup_->setExclusive(true); Only make exclusive when a selection has been made
}
void nailfold_qseries_gui::createReasonsRadioGroup()
{
  // Add vessel row buttons to a group for convenience
  reasonsRadioGroup_ = new QButtonGroup(this);
	reasonsRadioGroup_->addButton(ui.reason1Radio, ncm_qseries_grade::Reason0);
  reasonsRadioGroup_->addButton(ui.reason2Radio, ncm_qseries_grade::Reason1);
	reasonsRadioGroup_->addButton(ui.reason3Radio, ncm_qseries_grade::Reason2);
	reasonsRadioGroup_->addButton(ui.reason4Radio, ncm_qseries_grade::Reason3);
	reasonsRadioGroup_->addButton(ui.reason5Radio, ncm_qseries_grade::Reason4);
	reasonsRadioGroup_->addButton(ui.reason6Radio, ncm_qseries_grade::Reason5);
	reasonsRadioGroup_->addButton(ui.reason7Radio, ncm_qseries_grade::Reason6);
	reasonsRadioGroup_->addButton(ui.reason8Radio, ncm_qseries_grade::Reason7);
	reasonsRadioGroup_->addButton(ui.reason9Radio, ncm_qseries_grade::Reason8);
	reasonsRadioGroup_->addButton(ui.reason10Radio, ncm_qseries_grade::Reason9);
	reasonsRadioGroup_->addButton(ui.reason11Radio, ncm_qseries_grade::Reason10);
	reasonsRadioGroup_->addButton(ui.reason12Radio, ncm_qseries_grade::Reason11);
	reasonsRadioGroup_->addButton(ui.reason13Radio, ncm_qseries_grade::Reason12);
  reasonsRadioGroup_->setExclusive(false);
}

//
//: Enable/disable controls that are dependent on user level (basic or advanced)
void nailfold_qseries_gui::applyUserLevelChange()
{
  const bool advanced_user = preferences_.is_advanced_user();
}

//
//: Update radio buttons to reflect current image grade.
void nailfold_qseries_gui::updateGradeRadios()
{
	if (grade_.is_graded())
	{
		gradeRadioGroup_->setExclusive(true);

		if (grade_.is_normal())
			ui.imagesNormalRadio->setChecked(true);
		else
			ui.imagesAbnormalRadio->setChecked(true);
	}
	else
		uncheckGradeRadios();

}

//
//: Uncheck all of the radio buttons corresponding to image grade
void nailfold_qseries_gui::uncheckGradeRadios()
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
//: Uncheck all of the radio reasons buttons
void nailfold_qseries_gui::uncheckReasonsRadios()
{
  for (int i = 0; i < reasonsRadioGroup_->buttons().count(); ++i)
  {
    reasonsRadioGroup_->buttons()[i]->setChecked(false);
  }
}

//
//: Update radio buttons to reflect currently selected reasons
void nailfold_qseries_gui::updateReasonsRadios()
{
	for (int i = 0; i < reasonsRadioGroup_->buttons().count(); ++i)
  {
		reasonsRadioGroup_->buttons()[i]->setChecked( grade_.progress_reasons(static_cast<ncm_qseries_grade::ProgressReason>(i) ) );
  }
}

//
//: Does what it says on the tin
void nailfold_qseries_gui::connectSignalsToSlots()
{
  QObject::connect( ui.autoContrast, SIGNAL(toggled(bool)),
                    ui.graphicsView, SLOT(setAutoContrast(bool)) );
	QObject::connect( ui.autoContrast_2, SIGNAL(toggled(bool)),
                    ui.graphicsViewBottom, SLOT(setAutoContrast(bool)) );
  QObject::connect( ui.autoContrast, SIGNAL(toggled(bool)),
                    this, SLOT(updateContrastControls()) );
	QObject::connect( ui.autoContrast_2, SIGNAL(toggled(bool)),
                    this, SLOT(updateContrastControls()) );

  // make connections between UI widgets and local member variables
  // we could possibly make these slots part of the graphicsview rather than the
  // scene so that the connections can be set up in QDesigner (but it's not a
  // big deal at the moment).

  QObject::connect( &imageProcessor_, SIGNAL(imageLoaded(bool)),
                    this, SLOT(onImageProcessor_imageLoaded(bool)) );
	QObject::connect( &imageProcessor2_, SIGNAL(imageLoaded(bool)),
                    this, SLOT(onImageProcessor2_imageLoaded(bool)) );

  QObject::connect( ui.graphicsView, SIGNAL(changed()),
                    this, SLOT(updateZoomControls()) );

	// Set grade
  QObject::connect( gradeRadioGroup_, SIGNAL(buttonClicked(int)),
                    this, SLOT(setGrade(int)) );

  // Set reasons for grade progression
  QObject::connect( reasonsRadioGroup_, SIGNAL(buttonClicked(int)),
                    this, SLOT(selectReasons(int)) );

	// Contrast/zoom adjustment
  QObject::connect( &imageProcessor_, SIGNAL(rawImageChanged()),
                    ui.graphicsView, SLOT(repaint()) );
  QObject::connect( &imageProcessor_, SIGNAL(rawImageChanged()),
                    this, SLOT(updateContrastControls()) );

	QObject::connect( &imageProcessor2_, SIGNAL(rawImageChanged()),
                    ui.graphicsViewBottom, SLOT(repaint()) );
  QObject::connect( &imageProcessor2_, SIGNAL(rawImageChanged()),
                    this, SLOT(updateContrastControls()) );


  // Zoom controls
  QObject::connect( ui.zoomFitButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(fit_both()));
  QObject::connect( ui.zoomFitWButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(fit_width()));
  QObject::connect( ui.zoomFitHButton, SIGNAL(clicked()),
                    ui.graphicsView, SLOT(fit_height()));

	QObject::connect( ui.zoomFitButton, SIGNAL(clicked()),
                    ui.graphicsViewBottom, SLOT(fit_both()));
  QObject::connect( ui.zoomFitWButton, SIGNAL(clicked()),
                    ui.graphicsViewBottom, SLOT(fit_width()));
  QObject::connect( ui.zoomFitHButton, SIGNAL(clicked()),
                    ui.graphicsViewBottom, SLOT(fit_height()));

  connectSetFocus();
}

//
//: Connect signals from all controls such that they return focus to the
//  graphicsView after every operation.
void nailfold_qseries_gui::connectSetFocus()
{
  // After any of the following events, return the focus to the
  // graphics view to process keyboard events
  /*QObject::connect( ui.brightnessSlider, SIGNAL(sliderReleased()),
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

  QObject::connect( ui.stackedWidget, SIGNAL(currentChanged(int)),
                    ui.graphicsView, SLOT(setFocus()));*/
}


