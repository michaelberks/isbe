#ifndef NCM_QCAMERA_GUI_H
#define NCM_QCAMERA_GUI_H

#include <QBitmap>
#include <QGraphicsPixmapItem>
#include <QPushButton>
#include <QGroupBox>
#include <QSlider>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QGraphicsScene>

#include <vcl_iostream.h>
#include <vcl_cmath.h>
#include <vcl_iomanip.h>
#include <vcl_string.h>
#include <vcl_sstream.h>
#include <vcl_iosfwd.h>
#include <vcl_vector.h>
#include <vcl_cstddef.h>
#include <vil/vil_image_view.h>

#include <nailfold/qtools/ncm_cameras/ncm_camera.h>
#include <nailfold/qtools/ncm_exposure_spinbox.h>
#include <nailfold/qtools/nailfold_graphics_view.h>

class ncm_qcamera_gui
{

public:
	ncm_qcamera_gui();

	void set_ui_controls(
									QPushButton *_loadCameraButton,
									QGroupBox *_cameraControls,
									QSlider *_exposureSlider,
									QExposureSpinBox *_exposureSpinBox,
									QCheckBox *_autoExposure,
									QSlider *_gainSlider,
									QDoubleSpinBox *_gainSpinBox,
									QCheckBox *_autoGain, 
									QComboBox *_frameRateSelection,
									QGroupBox *_displayControls,
									QSlider *_brightnessSlider,
									QSlider *_contrastSlider,
									QPushButton *_stopLiveButton );

	void set_camera( ncm_camera *_camera );

	void set_scene_controls(
									nailfold_graphics_view *_graphicsView,
									QGraphicsPixmapItem *_raw_pixmap_item,
									QVector<QRgb> *_raw_colour_table,
									QGraphicsScene *_scene,
									QImage *_qimage );

	void loadCameraButton_clicked();
	void stopLiveButton_clicked();

	//Display controls
	void brightnessSliderMoved(int value);
	void contrastSliderMoved(int value);

	//Camera controls
	void frameRateSelectionActivated (const QString &text);
	void gainSpinBoxChanged(double value);
	void gainSliderMoved(int value);
	void autoGainChanged(int state);
	void exposureSpinBoxChanged(double value);
	void exposureSliderMoved(int value);
	void autoExposureChanged(int state);
	void redraw_scene();

private:

	void update_display(int & contrast_min, int & contrast_max, int contrast, int brightness);
	

  //: Update raw colormap to enhance contrast
  void update_raw_colour_table();

	//: pixmap item for raw image in scene
  QGraphicsPixmapItem *raw_pixmap_item_;

  //: Colormap for raw images
  QVector<QRgb> *raw_colour_table_;

  //: canvas on which to draw stuff
  QGraphicsScene *scene_;

	//: image that gets drawn in canvas
  QImage *qimage_;

	//Camera object
	ncm_camera *camera_;

  //: grey value to map to 0 (black)
  int contrast_min_;

  //: grey value to map to 255 (white)
  int contrast_max_;
	
	//Flags:
	bool camera_connected_; //Are cameras connected/live?
	bool camera_live_;
	

	//UI widgets
	nailfold_graphics_view *graphicsView_;
	QPushButton *loadCameraButton_;
	QGroupBox *cameraControls_;

	QSlider *exposureSlider_;
	QExposureSpinBox *exposureSpinBox_;
  QCheckBox *autoExposure_;
  
	QSlider *gainSlider_;
	QDoubleSpinBox *gainSpinBox_;
  QCheckBox *autoGain_;
  
  QComboBox *frameRateSelection_;

  QGroupBox *displayControls_;
  QSlider *brightnessSlider_;
  QSlider *contrastSlider_;
  QPushButton *stopLiveButton_;
};


#endif // NCM_QCAMERA_GUI_H