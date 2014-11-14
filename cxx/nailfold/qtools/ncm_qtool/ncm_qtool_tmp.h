#ifndef NCM_QTOOL_TMP_H
#define NCM_QTOOL_TMP_H

#include <QMainWindow>
#include <QWheelEvent>
#include "ui_ncm_qtool_tmp.h"

#include <vil/vil_image_view.h>

class ncm_qtool_tmp : public QMainWindow
{
  Q_OBJECT

public:
  ncm_qtool_tmp(QWidget *parent = 0, Qt::WFlags flags = 0);
  ~ncm_qtool_tmp();

  //: Non UI members manipulated by the UI
  void get_frame();

private slots:
  void on_actionLoad_image_activated();
  void on_actionReset_activated();
  void on_actionSmooth_activated();
  void on_actionFind_lines_activated();
  void on_actionGaussian_deriv_activated();
  void on_actionMaxima_activated();
  void on_actionHarris_activated();
  void on_actionDoG_activated();

  void on_actionPan_activated();
  void on_actionMarkup_activated();

  void on_actionDebug_activated();

  void on_pushButton_clicked();

  void on_brightnessSlider_valueChanged(int value);
  void on_contrastSlider_valueChanged(int value);
  void on_lowerLimitSlider_valueChanged(int value);
  void on_upperLimitSlider_valueChanged(int value);
  void on_compSizeSlider_valueChanged(int value);
  void on_lowerLimitSlider_sliderPressed();
  void on_upperLimitSlider_sliderPressed();
  void on_lowerLimitSlider_sliderReleased();
  void on_upperLimitSlider_sliderReleased();

private:
  Ui::NailfoldMainWindow ui;

  QGraphicsScene scene_;

  //: brightness change (relative to original)
  float brightness_;

  //: contrast change (relative to original)
  float contrast_;

  int lower_limit_;
  int upper_limit_;

  vil_image_view<vxl_byte> image_;
  vil_image_view<vxl_byte> small_image_;
  vil_image_view<vxl_byte> normalized_;
  vil_image_view<vxl_byte> processed_;

  vil_image_view<int> component_label_image_;
  vcl_vector<unsigned> component_sizes_;
  vil_image_view<vxl_byte> component_size_image_;

  QVector<QRgb> raw_colour_table_;
  QVector<QRgb> component_colour_table_;

  QPointF preview_scene_centre_;

  void redraw_raw();
  void quick_preview_raw();

  void update_raw_colour_table();
  void update_component_colour_table(int value);

protected:

  //void timerEvent( QTimerEvent * );
};

#endif // ncm_qtool_TMP_H
