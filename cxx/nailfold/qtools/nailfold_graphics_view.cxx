//:
// \file
// \brief View including ability to zoom in and out
//        Based on Tim's qvcr_zoom_view
// \author Phil Tresadern

#include "nailfold_graphics_view.h"

#include <vcl_cmath.h>
#include <vcl_algorithm.h>
#include <vcl_iostream.h>

#include <QFileDialog>
#include <QMessageBox>
#include <QtGui>

nailfold_graphics_view::nailfold_graphics_view(QWidget * parent)
  : QGraphicsView(parent)
{
  set_defaults();
}

nailfold_graphics_view::nailfold_graphics_view(QGraphicsScene* scene,
                                 QWidget * parent)
  : QGraphicsView(scene,parent)
{
  set_defaults();
}

//: Destructor
nailfold_graphics_view::~nailfold_graphics_view()
{
  // delete all objects in the trash
  for (unsigned i = 0; i < trash_.size(); ++i)
    delete trash_[i];
}

//
//: Sets various properties
void nailfold_graphics_view::set_defaults()
{
  // Define general background colour from palette
  setBackgroundRole(QPalette::Mid);

  // Allow movement of viewport by dragging background
  setDragMode(QGraphicsView::ScrollHandDrag);

  // Default background is greenish
  // Don't use as occasionally get bits left behind
//  setBackgroundBrush(QBrush(QColor(10,100,10)));

  // Define default zoom direction to increase zoom when wheel
  // is rotated toward you
  zoom_dir_ = -1;

  // use snapping by default
  use_snapping_ = true;

  // Trigger mouseMoveEvent even if no buttons pressed
  setMouseTracking(true);

  vessel_pen_ = QPen(QColor(255,0,0));
  vessel_pen_.setWidth(2);

  selected_vessel_pen_ = QPen(QColor(255,255,0));
  selected_vessel_pen_.setWidth(2);

  // data is not modified at start
  is_modified_ = false;

  trash_.clear();
}

void nailfold_graphics_view::get_vessel_pixels(QGraphicsPixmapItem* pixmap_item)
{
  // if no pixmap supplied then do nothing (return)
  if (pixmap_item == 0)
    return;

  // clear current list
  nonzero_pixels_.resize(0);

  QImage image = pixmap_item->pixmap().toImage();
  for (int y = 0; y < image.height(); ++y)
    for (int x = 0; x < image.width(); ++x)
    {
      if (qGray(image.pixel(x,y)) > 0)
        nonzero_pixels_.push_back(QPoint(x,y));
    }
}

QPoint nailfold_graphics_view::pixel_nearest_to(QPointF pos)
{
  // search through list for nearest point
  // could be done more effectively by storing nonzero pixels as a
  // vector< vector<QPoint> > with one vector per row (or column) of the image
  // the list could then be searched starting from the current row/col,
  // propagating the search raster outward

  QPoint nearest_point(0,0);
  int min_distance = -1;

  for (unsigned i = 0; i < nonzero_pixels_.size(); ++i)
  {
    int d = (pos - nonzero_pixels_[i]).manhattanLength();

    if ((min_distance == -1) || (d < min_distance))
    {
      nearest_point = nonzero_pixels_[i];
      min_distance = d;
    }
  }

  return nearest_point;
}


void nailfold_graphics_view::clear_nearby_lines()
{
  // reset lines that had previously changed colour back to red
  for (unsigned i = 0; i < nearby_lines_.size(); ++i)
    nearby_lines_[i]->setPen(vessel_pen_);
  nearby_lines_.resize(0);
}

//
//: Add a point to the current vessel, snapping where needed and avoiding
//  duplicated points
void nailfold_graphics_view::addVesselPoint(QPoint pos)
{
  if (use_snapping_)
  {
    QPoint new_p = pixel_nearest_to(mapToScene(pos));

    if (polyline_.empty())
    {
      // if this is the first then add
      polyline_ << new_p;
    }
    else
    {
      // avoid duplicating points
      QPoint last_p = polyline_[polyline_.size()-1].toPoint();
      if ((new_p.x() != last_p.x()) ||
          (new_p.y() != last_p.y()))
        polyline_ << new_p;
    }
  }
  else
    polyline_ << mapToScene(pos);
}

// Event handlers
//
//: React to key presses
void nailfold_graphics_view::keyPressEvent(QKeyEvent *event)
{
  switch (event->key()) {
    case Qt::Key_Up:
//         translate(0, -20);
        break;
    case Qt::Key_Down:
//         translate(0, 20);
        break;
    case Qt::Key_Left:
//         translate(-20, 0);
        break;
    case Qt::Key_Right:
//         translate(20, 0);
        break;
    case Qt::Key_Plus:
        scaleView(1.2);
        break;
    case Qt::Key_Minus:
        scaleView(1 / 1.2);
        break;
    default:
        QGraphicsView::keyPressEvent(event);
  }
}

//
//: React to mouse button press
void nailfold_graphics_view::mousePressEvent(QMouseEvent *event)
{
  // let QGraphicsView handle Pan mode
  if (dragMode() == QGraphicsView::ScrollHandDrag)
  {
    event->ignore();
    QGraphicsView::mousePressEvent(event);
    return;
  }

  switch (event->button())
  {
  case Qt::LeftButton:
    // begin drawing line
    polyline_.clear();
    addVesselPoint(event->pos());

    if (scene() != 0)
    {
      QPainterPath qpath;
      qpath.addPolygon(polyline_);
      current_polyline_ = scene()->addPath(qpath, vessel_pen_);
      is_modified_ = true;
    }
    break;
  }
}
//
//: React to mouse movement
void nailfold_graphics_view::mouseMoveEvent(QMouseEvent *event)
{
  // let QGraphicsView handle Pan mode
  if (dragMode() == QGraphicsView::ScrollHandDrag)
  {
    event->ignore();
    QGraphicsView::mouseMoveEvent(event);
    return;
  }

  clear_nearby_lines();

  QGraphicsPathItem* nearby_path = 0;

  switch (event->buttons())
  {
  case Qt::NoButton:
    // if there is more than just the pixmap then look for nearby lines
    nearby_path = vessel_nearest_to(event->pos());
    if (nearby_path != 0)
    {
      nearby_path->setPen(selected_vessel_pen_);
      nearby_lines_.push_back(nearby_path);
    }
    break;

  case Qt::LeftButton:
    // drawing line
    addVesselPoint(event->pos());

    if (scene() != 0)
    {
      // get rid of old polyline and delete object created by previous call
      // to scene()->addPath()
      scene()->removeItem(current_polyline_);
      delete current_polyline_;

      QPainterPath qpath;
      qpath.addPolygon(polyline_);
      current_polyline_ = scene()->addPath(qpath, vessel_pen_);

      is_modified_ = true;
    }
  }
}

//
//: React to mouse button release
void nailfold_graphics_view::mouseReleaseEvent(QMouseEvent *event)
{
  // let QGraphicsView handle Pan mode
  if (dragMode() == QGraphicsView::ScrollHandDrag)
  {
    event->ignore();
    QGraphicsView::mouseReleaseEvent(event);
    return;
  }

  QGraphicsPathItem* nearby_path = 0; 
  
  switch (event->button())
  {
  case Qt::LeftButton:
    // end line drawing
    addVesselPoint(event->pos());

    if (scene() != 0)
    {
      // get rid of old polyline and delete object created by previous call
      // to scene()->addPath()
      scene()->removeItem(current_polyline_);
      delete current_polyline_;
      current_polyline_ = 0;

      QPainterPath qpath;
      qpath.addPolygon(polyline_);
      QGraphicsPathItem* path_item = scene()->addPath(qpath, vessel_pen_);

      // add a tooltip to indicate number of points
      int length = path_item->path().elementCount();
      path_item->setToolTip("#points: "+QString::number(length));

      is_modified_ = true;
    }
    break;
 
  case Qt::RightButton:
    // delete nearest line
    nearby_path = vessel_nearest_to(event->pos());
    if (nearby_path != 0)
    {
      scene()->removeItem(nearby_path);
      trash_.push_back(nearby_path);
      nearby_path = 0;

      is_modified_ = true;
    }
    break;
  }
}

//
//: Use mouse-wheel movement to scale view
void nailfold_graphics_view::wheelEvent(QWheelEvent *event)
{
  if (dragMode() == QGraphicsView::ScrollHandDrag)
    scaleView(vcl_pow((double)2, zoom_dir_*event->delta() / 240.0));
}

//
//: When widget resized, update view transform
void nailfold_graphics_view::resizeEvent ( QResizeEvent * event )
{
  if (event->oldSize().isEmpty()) 
  {
    show_all_scene();
    return;
  }

  double x_scale = double(event->size().width())
                          /event->oldSize().width();
  double y_scale = double(event->size().height())
                          /event->oldSize().height();

  // Find the smallest scaling, so as to ensure that 
  // the whole of the current view remains in view
  double s = vcl_min(x_scale,y_scale);

  scaleView(s);
}



// Private functions
//
//: Scale the view by given factor
void nailfold_graphics_view::scaleView(qreal scaleFactor)
{
  qreal factor = matrix().scale(scaleFactor, scaleFactor).mapRect(QRectF(0, 0, 1, 1)).width();
  if (factor < 0.07 || factor > 100)
      return;

  scale(scaleFactor, scaleFactor);

  if (scene())
  {
    // Prevent zooming out if whole scene already visible
    QPolygon scene_region = mapFromScene(scene()->sceneRect());
    QRect view_rect = rect();
    if (view_rect.contains(scene_region.point(0)) &&
        view_rect.contains(scene_region.point(1)) &&
        view_rect.contains(scene_region.point(2)) &&
        view_rect.contains(scene_region.point(3))    )
    {
      // View is already fully visible
      fitInView(scene()->sceneRect(),Qt::KeepAspectRatio);
    }
  }

}

//
//: Set view to display all of scene rectangle
void nailfold_graphics_view::show_all_scene()
{
  if (scene())
  {
    fitInView(scene()->sceneRect(),Qt::KeepAspectRatio);
  }
}


//
//: Find marked vessel nearest to pos
QGraphicsPathItem* nailfold_graphics_view::vessel_nearest_to(QPoint pos)
{
  QGraphicsPathItem* nearby_path = 0;

  // if there is more than just the pixmap then look for nearby lines
  if ((scene() != 0) && (scene()->items().size() > 1))
  {
    // add invisible disc of radius 5
    QGraphicsItem* disc = scene()->addEllipse(mapToScene(pos).x()-2.5, 
                                              mapToScene(pos).y()-2.5,
                                              5, 5);

    // find nearby items (i.e. that overlap with the disc)
    QList<QGraphicsItem*> nearby = scene()->collidingItems(disc);

    // find highest nearby line on the pile
    QList<QGraphicsItem*>::iterator item_it = nearby.begin();
    while ( (nearby_path == 0) && (item_it != nearby.end()) )
    {
      nearby_path = qgraphicsitem_cast<QGraphicsPathItem*>(*item_it);
      ++item_it;
    }

    // if not found then return zero later
    if ( item_it == nearby.end() )
      nearby_path = 0;

    // get rid of the disc
    scene()->removeItem(disc);
    delete disc;
    disc = 0;
  }

  return nearby_path;
}

//
//: Invoke file dialog and save screenshot
void nailfold_graphics_view::save_screenshot()
{
  QString filename = QFileDialog::getSaveFileName(this);
  if (!filename.isEmpty())
    save_screenshot(filename);
}

//
//: Save screenshot to named file
void nailfold_graphics_view::save_screenshot(const QString& image_path)
{
  // By writing to a QImage we get anti-aliasing, which isn't
  // always present on screen.

  QImage buffer_image( size(), QImage::Format_RGB32 );
  QPainter image_painter(&buffer_image);
  // Ensure background is drawn:
  image_painter.setBackgroundMode(Qt::OpaqueMode);  
  render(&image_painter);
  image_painter.end();
  if (!buffer_image.save(image_path))
    QMessageBox::warning(this,
        tr("Screenshot"),
        tr("Failed to save screenshot to\n")
        +image_path);

/*
  QPixmap pixmap=QPixmap::grabWidget(this);
  if (!pixmap.save(image_path))
    QMessageBox::warning(this,
        tr("Screenshot"),
        tr("Failed to save screenshot to\n")
        +image_path);
*/
}

