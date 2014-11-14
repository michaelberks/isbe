#include "ncm_qapt_util_gui.h"

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

#include <vcl_cmath.h>
#include <vcl_iostream.h>

//#include <vnl/vnl_random.h>

#include <vul/vul_file.h>

#include <nailfold/ncm_apt_server.h>
#include <nailfold/ncm_apt_controller.h>

//
//: Constructor
ncm_qapt_util_gui::ncm_qapt_util_gui(QWidget *parent, Qt::WFlags flags)
: QMainWindow(parent, flags),
  blob_size_(1.0),
  motors_are_live_(true),
  x_precision_(0.1f),
  y_precision_(0.1f),
  z_precision_(0.01f)
{
  // Name objects to allow autoconnection of signals to slots
  scene_.setObjectName("scene");

  // setup the UI
  ui.setupUi(this);

  // tell graphics view to use this scene
  ui.graphicsView->setScene(&scene_);

  apt_.set_event_dialogs(true);

  connect_signals_to_slots();

  initialize_timer();

  // Flag to show that constructor ended successfully
  constructed_ = true;
}

//
//: Destructor
ncm_qapt_util_gui::~ncm_qapt_util_gui()
{
  if (motors_are_live_)
    apt_.move_all_home(true);
}

bool ncm_qapt_util_gui::initialize()
{
  // Set up sliders and spinboxes
  update_controlsX();
  update_controlsY();
  update_controlsZ();

  initialize_scene();
  update_scene();
  update_comboboxes();

  return true;
}


//
// Public slots
//

//
// Private slots
//

void ncm_qapt_util_gui::on_controllerX_activated(int index)
{
  // If combobox is empty (i.e. the list has just been cleared) then do nothing.
  switch (index)
  {
    case -1:
      ui.homeX->setEnabled(false);
      return;
    case 0:
      ui.homeX->setEnabled(false);
      apt_.set_X_id(-1);
      break;
    default:
      ui.homeX->setEnabled(true);
      apt_.set_X_id(ui.controllerX->currentText().toLong());
  }

  update_comboboxes();
  set_targetX(apt_.X()->position());
}
void ncm_qapt_util_gui::on_controllerY_activated(int index)
{
  // If combobox is empty (i.e. the list has just been cleared) then do nothing.
  switch (index)
  {
    case -1:
      ui.homeY->setEnabled(false);
      return;
    case 0:
      ui.homeY->setEnabled(false);
      apt_.set_Y_id(-1);
      break;
    default:
      ui.homeY->setEnabled(true);
      apt_.set_Y_id(ui.controllerY->currentText().toLong());
  }

  update_comboboxes();
  set_targetY(apt_.Y()->position());
}
void ncm_qapt_util_gui::on_controllerZ_activated(int index)
{
  // If combobox is empty (i.e. the list has just been cleared) then do nothing.
  switch (index)
  {
    case -1:
      ui.homeZ->setEnabled(false);
      return;
    case 0:
      ui.homeZ->setEnabled(false);
      apt_.set_Z_id(-1);
      break;
    default:
      ui.homeZ->setEnabled(true);
      apt_.set_Z_id(ui.controllerZ->currentText().toLong());
  }

  update_comboboxes();
  set_targetZ(apt_.Z()->position());
}

//
//:
void ncm_qapt_util_gui::on_homeX_clicked()
{
  if (motors_are_live_)
  {
    apt_.X()->set_home_offset(apt_.X()->mid_position());
    apt_.X()->move_home(false);
  }

  update_controlsX();
  update_scene();
}
void ncm_qapt_util_gui::on_homeY_clicked()
{
  apt_.Y()->set_home_offset(apt_.Y()->mid_position());
  apt_.Y()->move_home();
  update_controlsY();
  update_scene();
}
void ncm_qapt_util_gui::on_homeZ_clicked()
{
  apt_.Z()->set_home_offset(apt_.Z()->mid_position());
  apt_.Z()->move_home();
  update_controlsZ();
  update_scene();
}

//
//:
void ncm_qapt_util_gui::on_sliderX_valueChanged(int value)
{
  set_targetX(static_cast<double>(value) * x_precision_);
}
void ncm_qapt_util_gui::on_sliderY_valueChanged(int value)
{
  set_targetY(static_cast<double>(value) * y_precision_);
}
void ncm_qapt_util_gui::on_dialZ_valueChanged(int value)
{
  set_targetZ(static_cast<double>(value) * z_precision_);
}

//
//:
void ncm_qapt_util_gui::on_spinX_valueChanged(double value)
{
  //set_targetX(value);
}
void ncm_qapt_util_gui::on_spinY_valueChanged(double value)
{
  //set_targetY(value);
}
void ncm_qapt_util_gui::on_spinZ_valueChanged(double value)
{
  //set_targetZ(value);
}

//
//:
void ncm_qapt_util_gui::on_graphicsView_clicked(double x, double y)
{
  set_targetX(x);
  set_targetY(y);
}



//
//: Update the scene at regular intervals determined by a timer.
void ncm_qapt_util_gui::onPollTimer()
{
  update_scene();
}

//
// Private member functions
//

//
//: Update controls
void ncm_qapt_util_gui::update_controlsX()
{
  // Set up sliders and spinboxes
  ncm_apt_controller_base* c = apt_.X();
  
  ui.sliderX->setMinimum(static_cast<int>(c->min_position() / x_precision_));
  ui.sliderX->setMaximum(static_cast<int>(c->max_position() / x_precision_));
  ui.sliderX->setValue(static_cast<int>(c->position() / x_precision_));
  
  QDoubleSpinBox* sp = ui.spinX;
  sp->setMinimum(c->min_position());
  sp->setMaximum(c->max_position());
  sp->setValue(c->position());
  sp->setSingleStep(x_precision_);
  sp->setDecimals(-vcl_log10(x_precision_));
}

void ncm_qapt_util_gui::update_controlsY()
{
  // Set up sliders and spinboxes
  ncm_apt_controller_base* c = apt_.Y();

  ui.sliderY->setMinimum(static_cast<int>(c->min_position() / y_precision_));
  ui.sliderY->setMaximum(static_cast<int>(c->max_position() / y_precision_));
  ui.sliderY->setValue(static_cast<int>(c->position() / y_precision_));
  
  QDoubleSpinBox* sp = ui.spinY;
  sp->setMinimum(c->min_position());
  sp->setMaximum(c->max_position());
  sp->setValue(c->position());
  sp->setSingleStep(y_precision_);
  sp->setDecimals(-vcl_log10(y_precision_));
}

void ncm_qapt_util_gui::update_controlsZ()
{
  // Set up sliders and spinboxes
  ncm_apt_controller_base* c = apt_.Z();

  ui.dialZ->setMinimum(static_cast<int>(c->min_position() / z_precision_));
  ui.dialZ->setMaximum(static_cast<int>(c->max_position() / z_precision_));
  ui.dialZ->setValue(static_cast<int>(c->position() / z_precision_));
  
  QDoubleSpinBox* sp = ui.spinZ;
  sp->setMinimum(c->min_position());
  sp->setMaximum(c->max_position());
  sp->setValue(c->position());
  sp->setSingleStep(z_precision_);
  sp->setDecimals(-vcl_log10(z_precision_));
}

//
//: Create scene elements that we will reuse later.
void ncm_qapt_util_gui::update_scene_limits()
{
  const qreal left = ui.spinX->minimum();
  const qreal top = ui.spinY->minimum();
  const qreal width = ui.spinX->maximum() - left;
  const qreal height = ui.spinY->maximum() - top;

  scene_.setSceneRect(left, top, width, height);
  ui.graphicsView->fitInView(ui.graphicsView->sceneRect());
}
void ncm_qapt_util_gui::initialize_scene()
{
  const QRectF target_rect(-0.5*blob_size_, -0.5*blob_size_, 
                           blob_size_, blob_size_);
  QPen target_pen(QColor(0,0,0));
  QBrush target_brush(QColor(0,0,255));
  target_pos_item_ = scene_.addRect(target_rect, target_pen, target_brush);

  const QRectF actual_rect(-0.25*blob_size_, -0.25*blob_size_, 
                           0.5*blob_size_, 0.5*blob_size_);
  QPen actual_pen(QColor(0,0,0));
  QBrush actual_brush(QColor(255,0,0));
  actual_pos_item_ = scene_.addEllipse(actual_rect, actual_pen, actual_brush);

  update_scene_limits();
}

void ncm_qapt_util_gui::update_scene()
{
  update_scene_limits();

  const double target_x = ui.spinX->value();
  const double target_y = ui.spinY->value();
  target_pos_item_->setPos(target_x, target_y);

  const double actual_x = apt_.X()->position();
  const double actual_y = apt_.Y()->position();
  actual_pos_item_->setPos(actual_x, actual_y);
}

void ncm_qapt_util_gui::initialize_timer()
{
  poll_timer_.setInterval(50);
  poll_timer_.start();
}


//
//: Update controller comboboxes that enable you to assign controllers to axes.
void ncm_qapt_util_gui::update_comboboxes()
{
  // Get list of IDs that are connected
  vcl_vector<long> id_vector = apt_.ids();

  // For each axis, remove the IDs of any controller currently assigned to 
  // another axis and set the current item index to that matching the currently
  // assigned controller (if any).

  ui.controllerX->clear();
  ui.controllerX->addItem("None");
  for (unsigned i = 0; i < id_vector.size(); ++i)
  {
    // Add controller to list if not already in use by another axis.
    if ((apt_.Y()->id() != id_vector[i]) && 
        (apt_.Z()->id() != id_vector[i]))
      ui.controllerX->addItem(QString::number(id_vector[i]));

    // Store the index of the currently assigned controller.
    if (apt_.X()->id() == id_vector[i])
      ui.controllerX->setCurrentIndex(ui.controllerX->count()-1);
  }
    
  ui.controllerY->clear();
  ui.controllerY->addItem("None");
  for (unsigned i = 0; i < id_vector.size(); ++i)
  {
    // Add controller to list if not already in use by another axis.
    if ((apt_.X()->id() != id_vector[i]) && 
        (apt_.Z()->id() != id_vector[i]))
      ui.controllerY->addItem(QString::number(id_vector[i]));

    // Store the index of the currently assigned controller.
    if (apt_.Y()->id() == id_vector[i])
      ui.controllerY->setCurrentIndex(ui.controllerY->count()-1);
  }

  ui.controllerZ->clear();
  ui.controllerZ->addItem("None");
  for (unsigned i = 0; i < id_vector.size(); ++i)
  {
    // Add controller to list if not already in use by another axis.
    if ((apt_.X()->id() != id_vector[i]) && 
        (apt_.Y()->id() != id_vector[i]))
      ui.controllerZ->addItem(QString::number(id_vector[i]));

    // Store the index of the currently assigned controller.
    if (apt_.Z()->id() == id_vector[i])
      ui.controllerZ->setCurrentIndex(ui.controllerZ->count()-1);
  }
}

//
//: Set the target position for X axis.
void ncm_qapt_util_gui::set_targetX(double x)
{
  int int_kx = static_cast<int>(x / x_precision_);
  ui.sliderX->setValue(int_kx);
  ui.spinX->setValue(x);
  if (motors_are_live_)
    apt_.X()->move_to(x);
}

void ncm_qapt_util_gui::set_targetY(double y)
{
  int int_ky = static_cast<int>(y / y_precision_);
  ui.sliderY->setValue(int_ky);
  ui.spinY->setValue(y);
  if (motors_are_live_)
    apt_.Y()->move_to(y);
}

void ncm_qapt_util_gui::set_targetZ(double z)
{
  int int_kz = static_cast<int>(z / z_precision_);
  ui.dialZ->setValue(int_kz);
  ui.spinZ->setValue(z);
  if (motors_are_live_)
    apt_.Z()->move_to(z);
}

//
//: Connect signals to slots (not surprisingly...)
void ncm_qapt_util_gui::connect_signals_to_slots()
{
  QObject::connect( &poll_timer_, SIGNAL(timeout()),
                    this, SLOT(onPollTimer()) );
}
