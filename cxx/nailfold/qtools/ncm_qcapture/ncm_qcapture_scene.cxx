#include "ncm_qcapture_scene.h"

#include <QPen>
#include <QColor>
#include <QImage>
#include <QGraphicsItem>
#include <QGraphicsPixmapItem>
#include <QGraphicsSceneMouseEvent>
#include <QWidget>
#include <vcl_iostream.h>
#include <vcl_cmath.h>

//  Forward declarations of helper functions
int limit(int input, int lower_limit, int upper_limit);

//
//  Public methods
//

ncm_qcapture_scene::ncm_qcapture_scene()
: image_(new QImage()),
  endpoint_size_(8),
  is_dragging_joystick_(false)
{
  raw_pixmap_item_ = addPixmap( QPixmap::fromImage(*image_) );

  // Initialize colour table
  update_raw_colour_table(0, 255);
	
	create_limits();
  create_crosshairs();
  create_joystick();
	create_text();

  //hide_joystick();
}
 
void ncm_qcapture_scene::set_image(QImage* image)
{
  image_ = image;

  update_pixmap();
}
 
void ncm_qcapture_scene::move_pixmap_to(int x, int y)
{
  raw_pixmap_item_->setPos(x, y);
}
 
void ncm_qcapture_scene::show_crosshairs()
{
  horizontal_crosshair_->show();
  vertical_crosshair_->show();
}
 
void ncm_qcapture_scene::hide_crosshairs()
{
  horizontal_crosshair_->hide();
  vertical_crosshair_->hide();
}
 
void ncm_qcapture_scene::show_joystick()
{
  joystick_line_->show();
  joystick_endpoint_->show();
  joystick_limit_line_->show();
}
 
void ncm_qcapture_scene::hide_joystick()
{
  joystick_line_->hide();
  joystick_endpoint_->hide();
  joystick_limit_line_->hide();
}

void ncm_qcapture_scene::show_x_left_limit()
{
	left_border_->show();
}
void ncm_qcapture_scene::hide_x_left_limit()
{
	left_border_->hide();
}

void ncm_qcapture_scene::show_x_right_limit()
{
	right_border_->show();
}
void ncm_qcapture_scene::hide_x_right_limit()
{
	right_border_->hide();
}

void ncm_qcapture_scene::show_y_top_limit()
{
	top_border_->show();
}
void ncm_qcapture_scene::hide_y_top_limit()
{
	top_border_->hide();
}

void ncm_qcapture_scene::show_y_bottom_limit()
{
	bottom_border_->show();
}
void ncm_qcapture_scene::hide_y_bottom_limit()
{
	bottom_border_->hide();
}

void ncm_qcapture_scene::show_z_upper_limit()
{
	z_upper_msg_->show();
}
void ncm_qcapture_scene::hide_z_upper_limit()
{
	z_upper_msg_->hide();
}

void ncm_qcapture_scene::show_z_lower_limit()
{
	z_lower_msg_->show();
}
void ncm_qcapture_scene::hide_z_lower_limit()
{
	z_lower_msg_->hide();
}

void ncm_qcapture_scene::show_saving_msg()
{
	saving_frames_msg_->show();
}
void ncm_qcapture_scene::hide_saving_msg()
{
	saving_frames_msg_->hide();
}

void ncm_qcapture_scene::show_mosaic_msg()
{
	making_mosaic_msg_->show();
}
void ncm_qcapture_scene::hide_mosaic_msg()
{
	making_mosaic_msg_->hide();
}

//
//  Public slots
//

void ncm_qcapture_scene::set_contrast_from(int contrast, int brightness)
{
  int contrast_min = (255-brightness) - (255-contrast)/2;
  int contrast_max = (255-brightness) + (255-contrast)/2;

  contrast_min = limit(contrast_min, 0, 255);
  contrast_max = limit(contrast_max, 0, 255);

  update_raw_colour_table(contrast_min, contrast_max);

  update_pixmap();
}
 
//
//  Protected events
//

void ncm_qcapture_scene::mousePressEvent(
  QGraphicsSceneMouseEvent* mouseEvent)
{
  // Ignore anything not close to the centre of the image.
  if (                   (mouseEvent->scenePos().x() < endpoint_size_) && 
       (-endpoint_size_ < mouseEvent->scenePos().x()) && 
                         (mouseEvent->scenePos().y() < endpoint_size_) && 
       (-endpoint_size_ < mouseEvent->scenePos().y()) )
  {
    move_joystick_to(mouseEvent->scenePos());
    joystick_limit_line_->show();

    is_dragging_joystick_ = true;

    mouseEvent->accept();
  }
}

void ncm_qcapture_scene::mouseMoveEvent(
  QGraphicsSceneMouseEvent* mouseEvent)
{
  if (is_dragging_joystick_)
  {
    double max_length = joystick_limit();

    // Normalize (x,y) coordinates.
    double x = mouseEvent->scenePos().x() / max_length;
    double y = mouseEvent->scenePos().y() / max_length;
    double length = vcl_sqrt(x*x + y*y);

    // Limit to within a fixed circle.
    if (length > 1.0)
    {
      x /= length;
      y /= length;
    }

    emit joystickMoved(x, y);

    move_joystick_to( QPointF(x*max_length, y*max_length) );
    
    mouseEvent->accept();
  }
}

void ncm_qcapture_scene::mouseReleaseEvent(
  QGraphicsSceneMouseEvent* mouseEvent)
{
  if (is_dragging_joystick_)
  {
    emit joystickReleased();

    move_joystick_to(QPointF(0,0));
    joystick_limit_line_->hide();

    is_dragging_joystick_ = false;

    mouseEvent->accept();
  }
}
 
void ncm_qcapture_scene::mouseDoubleClickEvent(
  QGraphicsSceneMouseEvent* mouseEvent)
{
  emit doubleClicked( mouseEvent->scenePos().x(),
                      mouseEvent->scenePos().y(),
                      mouseEvent->button() );
}

void ncm_qcapture_scene::wheelEvent(QGraphicsSceneWheelEvent* wheelEvent)
{
	emit wheelMoved( wheelEvent->delta() );
}
 
//
//  Private methods
//

bool ncm_qcapture_scene::image_is_valid()
{
  return (image_ != NULL) && !(image_->isNull());
}

void ncm_qcapture_scene::update_raw_colour_table(
  int contrast_min, 
  int contrast_max)
{
  raw_colour_table_.resize(256);

  // fill everything from first to constrast_min_ with 0 (black)
  for (int i = 0; i < contrast_min; i++)
    raw_colour_table_[i] = qRgb(0,0,0);

  // fill everything from constrast_max_ to end with 255 (white)
  for (int i = contrast_max; i < 256; i++)
    raw_colour_table_[i] = qRgb(255,255,255);

  // fill values in between with a linear gradient
  int denom = contrast_max - contrast_min;
  for (int i = contrast_min; i < contrast_max; i++)
  {
    int grey = (255 * (i - contrast_min)) / denom;
    raw_colour_table_[i] = qRgb(grey, grey, grey);
  }

  update_pixmap();
}

//: Add the qimage to the canvas.
void ncm_qcapture_scene::update_pixmap()
{
  if ( image_is_valid() )
  {
	  image_->setColorTable(raw_colour_table_);

	  raw_pixmap_item_->setPixmap( QPixmap::fromImage(*image_) );
    raw_pixmap_item_->setPos(-image_->width()/2, -image_->height()/2);

	  setSceneRect(QRectF(-image_->width()/2, -image_->height()/2, 
                        image_->width(), image_->height() ) );

    update_crosshairs();

    update();
  }
}

void ncm_qcapture_scene::create_crosshairs()
{
  QPen crosshair_pen(QColor(0, 255, 0, /* alpha = */ 128));

  const qreal halfWidth = width() * 0.5;
  const qreal halfHeight = height() * 0.5;

  horizontal_crosshair_ = addLine(-halfWidth, 0, halfWidth, 0, crosshair_pen);
  vertical_crosshair_ = addLine(0, -halfHeight, 0, halfHeight, crosshair_pen);
}
 
void ncm_qcapture_scene::update_crosshairs()
{
  const qreal halfWidth = width() * 0.5;
  const qreal halfHeight = height() * 0.5;

  if (horizontal_crosshair_ != NULL)
    horizontal_crosshair_->setLine(-halfWidth, 0, halfWidth, 0);

  if (vertical_crosshair_ != NULL)
    vertical_crosshair_->setLine(0, -halfHeight, 0, halfHeight);
}

void ncm_qcapture_scene::create_joystick()
{
  QColor joystick_colour = QColor(255, 255, 0, /* alpha = */ 128);

  joystick_line_ = addLine(0, 0, 0, 0, 
                           QPen(joystick_colour));

  joystick_endpoint_ = addEllipse(-0.5*endpoint_size_, -0.5*endpoint_size_, 
                                  endpoint_size_, endpoint_size_, 
                                  QPen(joystick_colour), 
                                  QBrush(joystick_colour));

  double lim = joystick_limit();
  joystick_colour = QColor(255, 255, 0, /* alpha = */ 32);
  joystick_limit_line_ = addEllipse(-lim, -lim, 2*lim, 2*lim, 
                                    QPen(joystick_colour));

}
 
void ncm_qcapture_scene::move_joystick_to(QPointF pos)
{
  if (joystick_line_ != NULL) 
    joystick_line_->setLine(0, 0, pos.x(), pos.y());

  if (joystick_endpoint_ != NULL) 
    joystick_endpoint_->setRect(pos.x() - 0.5*endpoint_size_, 
                                pos.y() - 0.5*endpoint_size_, 
                                endpoint_size_, endpoint_size_);

  if (joystick_limit_line_ != NULL)
  {
    double lim = joystick_limit();
    joystick_limit_line_->setRect(-lim, -lim, 2*lim, 2*lim);
  }
}
 
double ncm_qcapture_scene::joystick_limit() const
{
  if (sceneRect().height() < sceneRect().width())
    return sceneRect().height() / 2;
  else
    return sceneRect().width() / 2;
}

void ncm_qcapture_scene::create_limits()
{
	const qreal halfWidth = 320;// width() * 0.5;
  const qreal halfHeight = 240;//height() * 0.5;

	QPen p;
	p.setColor(QColor(255,0,0));
	p.setWidth(10);
	
	left_border_ =		addLine(-halfWidth, -halfHeight, -halfWidth, halfHeight, p);
	right_border_ =		addLine(halfWidth, -halfHeight, halfWidth, halfHeight, p);
	top_border_ =			addLine(-halfWidth, -halfHeight, halfWidth, -halfHeight, p);
	bottom_border_ =	addLine(-halfWidth, halfHeight, halfWidth, halfHeight, p);

	left_border_->setOpacity(0.5);
	right_border_->setOpacity(0.5);
	top_border_->setOpacity(0.5);
	bottom_border_->setOpacity(0.5);

	left_border_->hide();
	right_border_->hide();
	top_border_->hide();
	bottom_border_->hide();
}

void ncm_qcapture_scene::create_text()
{
	const qreal halfWidth = 320;// width() * 0.5;
  const qreal halfHeight = 240;//height() * 0.5;

	QBrush b(QColor(255,0,0,128)); //Semi-opaque red lettering

	z_upper_msg_ = addSimpleText("Z upper limit reached", QFont("Arial", 24));
	QRectF r = z_upper_msg_->boundingRect();
	z_upper_msg_->setPos(-r.center().x(),-r.center().y());
	z_upper_msg_->setBrush(b);
	z_upper_msg_->hide();

	z_lower_msg_ = addSimpleText("Z lower limit reached", QFont("Arial", 24));
	r = z_lower_msg_->boundingRect();

	z_lower_msg_->setPos(-r.center().x(),-r.center().y());
	z_lower_msg_->setBrush(b);
	z_lower_msg_->hide();

	//Green opaque lettering in top left corner
	b.setColor(QColor(0,255,0,255));
	saving_frames_msg_ = addSimpleText("Saving frames...", QFont("Arial", 24));
	saving_frames_msg_->setPos(-halfWidth,-halfHeight);
	saving_frames_msg_->setBrush(b);
	saving_frames_msg_->hide();

	making_mosaic_msg_ = addSimpleText("Making mosaic...", QFont("Arial", 24));
	making_mosaic_msg_->setPos(-halfWidth,-halfHeight);
	making_mosaic_msg_->setBrush(b);
	making_mosaic_msg_->hide();

}

void ncm_qcapture_scene::rotate_text()
{
	QPointF p;
	p = z_upper_msg_->pos();
	z_upper_msg_->setPos(0,0);
	z_upper_msg_->rotate(180);
	z_upper_msg_->setPos(-p);

	p = z_lower_msg_->pos();
	z_lower_msg_->setPos(0,0);
	z_lower_msg_->rotate(180);
	z_lower_msg_->setPos(-p);

	p = saving_frames_msg_->pos();
	saving_frames_msg_->setPos(0,0);
	saving_frames_msg_->rotate(180);
	saving_frames_msg_->setPos(-p);

	p = making_mosaic_msg_->pos();
	making_mosaic_msg_->setPos(0,0);
	making_mosaic_msg_->rotate(180);
	making_mosaic_msg_->setPos(-p);
}

//
//  Helper functions
//

int limit(int input, int lower_limit, int upper_limit)
{
	if (input < lower_limit)
    return lower_limit;
  else if (input > upper_limit)
    return upper_limit;
  else
    return input;
}
