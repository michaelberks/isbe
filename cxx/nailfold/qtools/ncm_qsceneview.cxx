//:
// \file
// \brief View including ability to zoom in and out
//        Based on Tim's qvcr_zoom_view
// \author Phil Tresadern

#include "ncm_qsceneview.h"

#include <vcl_cmath.h>
#include <vcl_algorithm.h>
#include <vcl_iostream.h>

#include <vnl/vnl_random.h>

#include <vgl/vgl_point_2d.h>

#include <nailfold/ncm_annotation.h>
#include <nailfold/ncm_vessel.h>

#include <nailfold/qtools/ncm_qscene.h>
#include <nailfold/qtools/ncm_qimagehandler.h>

#include <QFileDialog>
#include <QMessageBox>
#include <QtGui>

//
//: Constructor
ncm_qsceneview::ncm_qsceneview(QWidget * parent)
  : QGraphicsView(parent)
{
  set_defaults();
}

//
//: Constructor
ncm_qsceneview::ncm_qsceneview(QGraphicsScene* scene,
                                 QWidget * parent)
  : QGraphicsView(scene,parent)
{
  set_defaults();
}

//
//: Destructor
ncm_qsceneview::~ncm_qsceneview()
{
}

//
//: Sets various properties
void ncm_qsceneview::set_defaults()
{
  // Define general background colour from palette
  setBackgroundRole(QPalette::Mid);

  // Load cursor for marking up
  QPixmap cursor_pixmap;
  cursor_pixmap.load(":/cursor_round");
  markup_cursor_ = QCursor(cursor_pixmap);

  // Define default zoom direction to increase zoom when wheel
  // is rotated toward you
  zoom_dir_ = -1;

  // Trigger mouseMoveEvent even if no buttons pressed
  setMouseTracking(true);

  // Use Pan mode to begin with
  useMarkupMode();

  //// Have scrollbars always off
  //setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  //setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

  auto_contrast_ = true;

  // Very large value - effectively infinite
  max_scale_ = 1.0e9;
}

//: Return handle to ncm_qscene rather than default QGraphicsScene
ncm_qscene* ncm_qsceneview::ncm_scene()
{
  return dynamic_cast<ncm_qscene*>(scene());
}

//
//: Centre on a particular vessel
void ncm_qsceneview::centre_on_vessel(int vessel_index)
{
  focus_on_vessel(vessel_index);
}

//
//: Focus on a particular vessel
void ncm_qsceneview::focus_on_vessel(int vessel_index,
                                     double zoom /* = 40.0 */,
                                     double translate_y /* = 0.0 */)
{
  const ncm_vessel* vessel = ncm_scene()->annotation()->vessel(vessel_index);
  if (vessel == NULL)
    return;

  // get coordinates of vessel anchor point
  const vgl_point_2d<double>& anchor = vessel->anchor();

  // define rectangle that we want to be in the view
  const double left = anchor.x() - zoom;
  const double top = anchor.y() - (1-translate_y)*zoom;
  const double width = 2*zoom;
  const double height = 2*zoom;
  QRectF zoom_rect(left, top, width, height);

  // fit rectangle in view
  fitInView(zoom_rect, Qt::KeepAspectRatio);
  emit changed();
}

//
//: Change direction of zoom with respect to mouse wheel movement
//  Default zooms in when wheel rolled toward you - this switches to the (I
//  think) more intuitive system where the view zooms in when the wheel is 
//  rolled away from you (as in Google Maps, for example).
void ncm_qsceneview::reverse_zoom()
{ 
  zoom_dir_ = -zoom_dir_; 
}

//
//: Get current mode
ncm_qsceneview::ViewMode ncm_qsceneview::current_mode()
{
  return current_mode_;
}

//
//: Set to Pan mode
void ncm_qsceneview::usePanMode()
{
  // If already in this mode then quit, changing nothing
  if (current_mode_ == ModePan)
    return;

  current_mode_ = ModePan;

  // do not forward mouse events to scene
  setInteractive(false);

  setDragMode(QGraphicsView::ScrollHandDrag);

  emit viewModeChanged(static_cast<int>(current_mode_));
}

//
//: Set to Pan/Markup mode
void ncm_qsceneview::usePanMarkupMode()
{
  // If already in this mode then quit, changing nothing
  if (current_mode_ == ModePanMarkup)
    return;

  current_mode_ = ModePanMarkup;

  // do not forward mouse events to scene
  setInteractive(false);

  setDragMode(QGraphicsView::ScrollHandDrag);

  emit viewModeChanged(static_cast<int>(current_mode_));
}

//
//: Set to Markup mode
void ncm_qsceneview::useMarkupMode()
{
  // If already in this mode then quit, changing nothing
  if (current_mode_ == ModeMarkup)
    return;

  current_mode_ = ModeMarkup;

  // do not forward mouse events to scene
  setInteractive(true);

  setDragMode(QGraphicsView::NoDrag);
  setCursor(markup_cursor_);

  emit viewModeChanged(static_cast<int>(current_mode_));
}

bool ncm_qsceneview::isPanning() const
{
  return (current_mode_ == ModePan ||
          current_mode_ == ModePanMarkup);
}

//
//: Automatically set contrast
void ncm_qsceneview::setAutoContrast(bool auto_contrast /* = true */)
{
  auto_contrast_ = auto_contrast;

  if (auto_contrast_)
  {
    QRect view = mapToScene(rect()).boundingRect().toRect();
    ncm_scene()->image_processor()->set_contrast_from(view);
  }

  update();
}

//
//  Event handlers
//

//: React to key presses
void ncm_qsceneview::keyPressEvent(QKeyEvent *event)
{
  switch (event->key()) {
    case Qt::Key_Comma: // (same key as '<')
      if ((ncm_scene()->edit_mode() == ncm_qscene::ModeLabelApices) ||
          (ncm_scene()->edit_mode() == ncm_qscene::ModeDrawVesselPath))
      {
        emit previous();
        event->accept();
      }
      break;

    case Qt::Key_Period: // (same key as '>')
      if ((ncm_scene()->edit_mode() == ncm_qscene::ModeLabelApices) ||
          (ncm_scene()->edit_mode() == ncm_qscene::ModeDrawVesselPath))
      {
        emit next();
        event->accept();
      }
      break;

    case Qt::Key_Plus:
      scaleView(1.2);
      emit changed();
      break;

    case Qt::Key_Minus:
      scaleView(1 / 1.2);
      emit changed();
      break;

    case Qt::Key_Control:
      if (current_mode_ == ModeMarkup)
        usePanMarkupMode();
      break;

    default:
      QGraphicsView::keyPressEvent(event);
      emit changed();
  }
}

void ncm_qsceneview::keyReleaseEvent(QKeyEvent *event)
{
  switch (event->key()) {
    case Qt::Key_Control:
      if (current_mode_ == ModePanMarkup)
        useMarkupMode();
      break;

    default:
      QGraphicsView::keyPressEvent(event);
      emit changed();
  }
}

//
//: Use mouse-wheel movement to scale view
void ncm_qsceneview::wheelEvent(QWheelEvent *event)
{
  scaleViewAbout(event->pos(),
                 vcl_pow((double)2, zoom_dir_*event->delta() / 240.0));
}

//
//: Scale at which the image width fits in the window
double ncm_qsceneview::fitWidthScale() const
{
  const int sbo = 23; // to account for scrollbar offset
  const QRect view_rect = rect();
  return (view_rect.width() - sbo) / sceneRect().width();
}

//
//: Scale at which the image height fits in the window
double ncm_qsceneview::fitHeightScale() const
{
  const int sbo = 23; // to account for scrollbar offset
  const QRect view_rect = rect();
  return (view_rect.height() - sbo) / sceneRect().height();
}

//
//: Scale at which the image fits in the window
double ncm_qsceneview::fitScale() const
{
  return vcl_min(fitWidthScale(), fitHeightScale());
}

double ncm_qsceneview::maxScale() const
{
  return max_scale_;
}

//
//:
void ncm_qsceneview::setMaxScale(double max_scale)
{
  max_scale_ = max_scale;

  // Limit zoom to a maximum value
  const double current_scale = transform().m11();
  if (current_scale > max_scale_)
    setScale(max_scale_);
}

//
//:
void ncm_qsceneview::setScale(double new_scale)
{
  const double current_scale = transform().m11();

  const double ratio = new_scale / current_scale;
  QTransform new_transform = transform();
  new_transform.scale(ratio, ratio);

  setTransform(new_transform);

  // Emit signal only if scale has changed appreciably
  if (vcl_abs(ratio - 1.0) > 1e-12)
    emit changed();
}

/*
//
//: Set zoom factor in real terms (1.0 = fit in view)
void ncm_qsceneview::setZoomFactor(double zoomFactor)
{
  //xscale = (view_rect.width() - 23) / sceneRect().width();
  //yscale = (view_rect.height() - 23) / sceneRect().height();
  //m11 = vcl_min(xscale,yscale);

  //QTransform new_transform = transform();
  //new_transform.scale(amount);

  //setTransform(new_transform);
}

double ncm_qsceneview::zoomFactor() const
{
  // Get scene dimensions in screen coordinates
  const QPolygon scene_region = mapFromScene(scene()->sceneRect());
  const QPoint diagonal = scene_region.point(2) - scene_region.point(0);
  const qreal scene_width = diagonal.x();
  const qreal scene_height = diagonal.y();

  // Including the scrollbars, there is a 23 (!) pixel border between the view
  // of the scene and the frame of the GraphicsView. Accounting for this, we
  // can compute the zoom factor in real terms
  const QPoint scrollbar_offset(23,23);
  const QRect view_rect = rect();
  const qreal xscale = scene_width / (view_rect.width() - scrollbar_offset.x());
  const qreal yscale = scene_height / (view_rect.height() - scrollbar_offset.y());

  qreal xscale2 = (view_rect.width() - 23) / sceneRect().width();
  qreal yscale2 = (view_rect.height() - 23) / sceneRect().height();
  qreal m11 = vcl_min(xscale2,yscale2);
  vcl_cout << m11 << " vs " << transform().m11() << vcl_endl;

  // The dimension in which the scene fills the view will have the higher ratio
  // in which we are interested.
  return vcl_max(xscale, yscale);
}
*/

//
//: When widget resized, update view transform
void ncm_qsceneview::resizeEvent ( QResizeEvent * event )
{
  if (event->oldSize().isEmpty()) 
  {
    fit_both();
    return;
  }

  double x_scale = double(event->size().width()) / event->oldSize().width();
  double y_scale = double(event->size().height()) / event->oldSize().height();

  // Find the smallest scaling, to ensure that 
  // the whole of the current view remains in view
  double s = vcl_min(x_scale,y_scale);

  scaleView(s);
}

void ncm_qsceneview::paintEvent(QPaintEvent* event)
{
  //static vnl_random rnd;
  //int r = rnd.lrand32();
  //vcl_cout << r << vcl_endl;
  //rnd.reseed(r);

  if (ncm_scene() != NULL && 
      ncm_scene()->image_processor() != NULL)
  {
    if (auto_contrast_)
    {
      QRect view = mapToScene(rect()).boundingRect().toRect();
      ncm_scene()->image_processor()->set_contrast_from(view);
    }

    // only repaint if the brightness/contrast are up-to-date
    if (ncm_scene()->image_processor()->isValid())
      QGraphicsView::paintEvent(event);
  }
}

//
// Private functions
//

//: Scale the view by given factor
void ncm_qsceneview::scaleView(qreal scaleFactor)
{
  scaleViewAbout(rect().center(), scaleFactor);
}

//: Scale the view by given factor about a given fixed point
void ncm_qsceneview::scaleViewAbout(const QPoint& fixed_point, qreal scaleFactor)
{
  // width of one image pixel in screen pixels
  QMatrix scaled_matrix = matrix().scale(scaleFactor, scaleFactor);
  qreal factor = scaled_matrix.mapRect(QRectF(0, 0, 1, 1)).width();
  
  // Limit zoom to a maximum value
  const double current_scale = transform().m11();
  const double new_scale = scaleFactor * current_scale;
  if (new_scale > max_scale_)
    scaleFactor *= (max_scale_ / new_scale);

  // Set new centre such that fixed_point remains static in the view
  // With scrollbars, we need a small offset to give the true centre of the
  // scene
  const QPoint scrollbar_offset(8,8);
  QPointF fixed_point_scene = mapToScene(fixed_point - scrollbar_offset);
  QPointF current_centre = mapToScene(rect().center() - scrollbar_offset);
  QPointF new_centre = fixed_point_scene + 
                          (current_centre - fixed_point_scene) / scaleFactor;

  scale(scaleFactor, scaleFactor);
  centerOn(new_centre);

  if (scene())
  {
    // Prevent zooming out if whole scene already visible
    QPolygon scene_region = mapFromScene(scene()->sceneRect());
    QRect view_rect = rect();

    if (view_rect.contains(scene_region.point(0)) &&
        view_rect.contains(scene_region.point(1)) &&
        view_rect.contains(scene_region.point(2)) &&
        view_rect.contains(scene_region.point(3)))
    {
      // View is already fully visible
      fitInView(scene()->sceneRect(), Qt::KeepAspectRatio);
    }
  }

  emit changed();
}

//
//: Set view to display all of scene rectangle
void ncm_qsceneview::fit_both()
{
  if (ncm_scene())
  {
    fitInView(ncm_scene()->imageRect(), Qt::KeepAspectRatio);
    emit changed();
  }
}

void ncm_qsceneview::fit_width()
{
  if (ncm_scene())
  {
    QRectF image_rect = ncm_scene()->imageRect();
    image_rect.setTop(0.5 * (image_rect.top()+image_rect.bottom()));
    image_rect.setHeight(0.00001); // must be nonzero
    fitInView(image_rect, Qt::KeepAspectRatio);
    emit changed();
  }
}

void ncm_qsceneview::fit_height()
{
  if (ncm_scene())
  {
    QRectF image_rect = ncm_scene()->imageRect();
    image_rect.setLeft(0.5 * (image_rect.left()+image_rect.right()));
    image_rect.setWidth(0.00001); // must be nonzero
    fitInView(image_rect, Qt::KeepAspectRatio);
    emit changed();
  }
}

//
//: Invoke file dialog and save screenshot
void ncm_qsceneview::save_screenshot()
{
  QString filename = QFileDialog::getSaveFileName(this);
  if (!filename.isEmpty())
    save_screenshot(filename);
}

//
//: Save screenshot to named file
void ncm_qsceneview::save_screenshot(const QString& image_path)
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

