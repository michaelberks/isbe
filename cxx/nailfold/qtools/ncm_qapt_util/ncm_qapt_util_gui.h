#ifndef NCM_QAPT_UTIL_GUI_H
#define NCM_QAPT_UTIL_GUI_H

#include "ui_ncm_qapt_util_gui.h"

#include <QMainWindow>
#include <QWheelEvent>
#include <QProgressBar>
#include <QThread>
#include <QDir>
#include <QTimer>

#include <nailfold/ncm_apt_server.h>

class ncm_qapt_util_gui : public QMainWindow
{
  Q_OBJECT

//
// INTERFACE
//

public:
  ncm_qapt_util_gui(QWidget *parent = 0, Qt::WFlags flags = 0);

  ~ncm_qapt_util_gui();

  bool initialize();

public slots:

//
// IMPLEMENTATION
//

protected:

private slots:

  void on_controllerX_activated(int index);
  void on_controllerY_activated(int index);
  void on_controllerZ_activated(int index);

  void on_homeX_clicked();
  void on_homeY_clicked();
  void on_homeZ_clicked();

  void on_sliderX_valueChanged(int value);
  void on_sliderY_valueChanged(int value);
  void on_dialZ_valueChanged(int value);

  void on_spinX_valueChanged(double value);
  void on_spinY_valueChanged(double value);
  void on_spinZ_valueChanged(double value);

  void on_graphicsView_clicked(double x, double y);

  void onPollTimer();

private:
  //
  // Private member functions
  //

  void update_controlsX();
  void update_controlsY();
  void update_controlsZ();

  void update_scene_limits();
  void initialize_scene();
  void update_scene();
  void initialize_timer();
  void update_comboboxes();

  void set_targetX(double x);
  void set_targetY(double y);
  void set_targetZ(double z);

  void connect_signals_to_slots();

  //
  // Private member variables
  //

  //: User interface object
  Ui::MainWindow ui;

  //: Flag set to true if the window was constructed completely
  bool constructed_;

  //: canvas on which to draw stuff
  QGraphicsScene scene_;
  QGraphicsItem* target_pos_item_;
  QGraphicsItem* actual_pos_item_;

  QTimer poll_timer_;

  ncm_apt_server apt_;

  const qreal blob_size_;
  const bool motors_are_live_;

  //: Precision of target positions in millimetres.
  const float x_precision_;
  const float y_precision_;
  const float z_precision_;
};

#endif // NCM_QAPT_UTIL_GUI_H
