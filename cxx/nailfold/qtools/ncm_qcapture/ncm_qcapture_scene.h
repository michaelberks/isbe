#ifndef ___NCM_QCAPTURE_SCENE_H__
#define ___NCM_QCAPTURE_SCENE_H__

#include <QGraphicsScene>

class QImage;
class QGraphicsItem;
class QGraphicsPixmapItem;

class ncm_qcapture_scene : public QGraphicsScene
{
  Q_OBJECT

public: // methods

  ncm_qcapture_scene();

  void set_image(QImage* image);

  void move_pixmap_to(int x, int y);

  void show_crosshairs();
  void hide_crosshairs();

  void show_joystick();
  void hide_joystick();

	void show_x_left_limit();
	void hide_x_left_limit();

	void show_x_right_limit();
	void hide_x_right_limit();

	void show_y_top_limit();
	void hide_y_top_limit();

	void show_y_bottom_limit();
	void hide_y_bottom_limit();

	void show_z_upper_limit();
	void hide_z_upper_limit();

	void show_z_lower_limit();
	void hide_z_lower_limit();

	void show_saving_msg();
	void hide_saving_msg();

	void show_mosaic_msg();
	void hide_mosaic_msg();
	
	void rotate_text();

public slots:

  void set_contrast_from(int contrast, int brightness);

signals:

  void joystickMoved(double, double);
  void joystickReleased();
  void doubleClicked(double, double, Qt::MouseButton);
  void wheelMoved(int);

protected: // events

  void mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent);
  void mouseMoveEvent(QGraphicsSceneMouseEvent* mouseEvent);
  void mouseReleaseEvent(QGraphicsSceneMouseEvent* mouseEvent);
  void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* mouseEvent);
  void wheelEvent(QGraphicsSceneWheelEvent* wheelEvent);

private: // methods

  bool image_is_valid();

  void update_raw_colour_table(
    int contrast_min, 
    int contrast_max);

  void update_pixmap();

  void create_crosshairs();
  void update_crosshairs();

  void create_joystick();
  void move_joystick_to(QPointF pos);

  double joystick_limit() const;

	void create_limits();
	void create_text();

private: // variables

  //: Pointer to the image to be displayed at the back of the scene
  QImage* image_;

  //: Pixmap item for raw image in scene
  QGraphicsPixmapItem* raw_pixmap_item_;

  //: Colormap for raw images
  QVector<QRgb> raw_colour_table_;

  //: Crosshair lines indicating the centre of the image.
  QGraphicsLineItem* horizontal_crosshair_;
  QGraphicsLineItem* vertical_crosshair_;

  //: 'Joystick' item that controls panning
  //ncm_qcapture_joystick joystick_;

  QGraphicsLineItem* joystick_line_;
  QGraphicsEllipseItem* joystick_endpoint_;
  QGraphicsEllipseItem* joystick_limit_line_;

	//: Boundary borders of the motor grid
	QGraphicsLineItem* left_border_;
	QGraphicsLineItem* right_border_;
	QGraphicsLineItem* top_border_;
	QGraphicsLineItem* bottom_border_;
	QGraphicsSimpleTextItem* z_upper_msg_;
	QGraphicsSimpleTextItem* z_lower_msg_;

	QGraphicsSimpleTextItem* saving_frames_msg_;
	QGraphicsSimpleTextItem* making_mosaic_msg_;

  int endpoint_size_;

  bool is_dragging_joystick_;

};

#endif ___NCM_QCAPTURE_SCENE_H__
