//:
// \file
// \brief QGraphicsItem-derived class for creating a visualization of the two
//        points either side of the apex (during annotation).
// \author Phil Tresadern

#include "ncm_qapexitem.h"
#include "ncm_qscene.h"
#include "ncm_qvesselitem.h"

#include <vcl_cmath.h>
#include <vcl_algorithm.h>
#include <vcl_iostream.h>

#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_vector_2d.h>

#include <QFileDialog>
#include <QMessageBox>
#include <QtGui>
#include <QPolygonF>
#include <QGraphicsItem>

#include <nailfold/ncm_annotation.h>
#include <nailfold/ncm_apex.h>

// This uses the source (rather than header) file so that QGraphicsApexItem's
// interface can remain hidden
#include "ncm_qapexitem_p.cxx"

/* static */ ncm_qapexitem_appearance& QGraphicsApexItem::appearance()
{
  return QGraphicsApexItemPrivate::appearance_;
}


/*!
    Constructs a QGraphicsPath. \a parent is passed to
    QAbstractGraphicsShapeItem's constructor.

    \sa QGraphicsScene::addItem()
*/
QGraphicsApexItem::QGraphicsApexItem(ncm_apex* const apex_data,
                                     bool interactive /* = false */,
                                     QGraphicsItem *parent /* = 0 */
#ifndef Q_QDOC
                                     // obsolete argument
                                     , QGraphicsScene *scene
#endif
)
: QGraphicsItem(*new QGraphicsApexItemPrivate, parent, scene)
{
  Q_D(QGraphicsApexItem);

  // FIXME: we should check that vessel_data is not zero
  d->set_apex_data(apex_data);

  if (!interactive)
    d->updateBoundingRect();
  
  // call hoverMoveEvent when the mouse hovers over the vessel
  setAcceptHoverEvents(true);

  // If the apex was created by a user clicking on the scene (interactive mode)
  // then use mousemove events to define the outer point
  if (interactive)
  {
    d->is_hovered_over_ = true;
    d->is_moving_outer_point_ = true;
  }
}

//
//: Destructor
//  Needs to take care of deleting the apex information from the markup, too
QGraphicsApexItem::~QGraphicsApexItem()
{
  Q_D(QGraphicsApexItem);

  assert( d->ncm_scene() != NULL );
  
  // Disconnect all signals
  disconnect();

  // Delete apex data via the vessel (to ensure that the vessel's list of 
  // apices is up to date)
  ncm_apex* const apex = d->apex_data();
  apex->parent_vessel().delete_apex(apex);

  // Remove apexitem from list stored within parent vesselitem
  d->vesselitem()->remove_apexitem(this);

  // Detach this vesselitem from the scene before it is destroyed
  d->ncm_scene()->removeItem(this);
}

//
//: What to do when the mouse is moved while no button is pressed
void QGraphicsApexItem::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
  Q_D(QGraphicsApexItem);
  d->cache_cursor(event->scenePos());
  d->is_hovered_over_ = true;
  d->updateBoundingRect();
}

void QGraphicsApexItem::hoverMoveEvent(QGraphicsSceneHoverEvent *event)
{
  Q_D(QGraphicsApexItem);
  d->cache_cursor(event->scenePos());
  d->is_hovered_over_ = true;
  d->updateBoundingRect();
}

void QGraphicsApexItem::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
  Q_D(QGraphicsApexItem);
  d->cache_cursor(event->scenePos());
  d->is_hovered_over_ = false;
  d->updateBoundingRect();
}

//
//: What to do when a mouse button is pressed
void QGraphicsApexItem::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  Q_D(QGraphicsApexItem);

  // unset accepted flag
  event->ignore();

  d->cache_cursor(event->scenePos());

  switch (event->button())
  {
    case Qt::LeftButton:
    {
      if (d->can_move_inner_point())
      {
        d->apex_data()->set_inner_point(d->cursor_pos_.x(), 
                                        d->cursor_pos_.y());
        d->is_moving_inner_point_ = true;
        event->accept();
      }
      else if (d->can_move_outer_point())
      {
        d->apex_data()->set_outer_point(d->cursor_pos_.x(), 
                                        d->cursor_pos_.y());
        d->is_moving_outer_point_ = true;
        event->accept();
      }
      
      break;
    }

    case Qt::RightButton:
    {
      if (d->can_be_deleted())
      {
        event->accept();
        delete this;
        // need to update list in vesselitem, too
        return;
      }
    }
  }

  d->updateBoundingRect();
}

//
//: What to do when the mouse is moved while a mouse button is pressed
void QGraphicsApexItem::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  Q_D(QGraphicsApexItem);

  // unset accepted flag
  event->ignore();

  d->cache_cursor(event->scenePos());

  if (d->is_moving_inner_point_)
  {
    d->apex_data()->set_inner_point(d->cursor_pos_.x(), 
                                    d->cursor_pos_.y());
    emit lengthChanged(d->apex_data()->width());
    event->accept();
  }
  else if (d->is_moving_outer_point_)
  {
    d->apex_data()->set_outer_point(d->cursor_pos_.x(), 
                                    d->cursor_pos_.y());
    emit lengthChanged(d->apex_data()->width());
    event->accept();
  }

  d->updateBoundingRect();
}

//
//: What to do when a mouse button is released
void QGraphicsApexItem::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
  Q_D(QGraphicsApexItem);

  // unset accepted flag
  event->ignore();

  d->cache_cursor(event->scenePos());

  if (d->is_moving_inner_point_)
  {
    d->apex_data()->set_inner_point(d->cursor_pos_.x(), 
                                    d->cursor_pos_.y());
    d->is_moving_inner_point_ = false;
    emit lengthChanged(-1.0); // Signal end of length changes
    event->accept();
  }
  else if (d->is_moving_outer_point_)
  {
    d->apex_data()->set_outer_point(d->cursor_pos_.x(), 
                                    d->cursor_pos_.y());
    d->is_moving_outer_point_ = false;
    emit lengthChanged(-1.0); // Signal end of length changes
    event->accept();
  }

  d->updateBoundingRect();

  // because we'll have made this item the grabber manually (when it was first 
  // created) we have to ungrab from it manually, too.
  ungrabMouse();

  // If apex is too small (usually as a result of not dragging) then delete it
  const double width_threshold = 5.0;
  if (d->apex_data()->width() < width_threshold)
    delete this;
}

/*!
    Returns the item's path as a QPainterPath. If no item has been set, an
    empty QPainterPath is returned.

    \sa setPath()
*/

QPainterPath QGraphicsApexItem::path() const
{
  Q_D(const QGraphicsApexItem);

  QPainterPath path;
  
  // connect inner and outer points
  path.moveTo(QPointF(d->apex_data()->inner_point().x(),
                      d->apex_data()->inner_point().y()));
  path.lineTo(QPointF(d->apex_data()->outer_point().x(),
                      d->apex_data()->outer_point().y()));

  // get normal to line joining inner and outer points
  vgl_vector_2d<double> connector = d->apex_data()->outer_point() - 
                                    d->apex_data()->inner_point();
  vgl_vector_2d<double> normal = 
      vgl_vector_2d<double>(-connector.y(), connector.x());
  normal *= 0.1;

  vgl_point_2d<double> p1, p2;

  // draw cap at inner point
  p1 = d->apex_data()->inner_point() + normal;
  p2 = d->apex_data()->inner_point() - normal;
  path.moveTo(QPointF(p1.x(), p1.y()));
  path.lineTo(QPointF(p2.x(), p2.y()));
                                      
  // draw cap at outer point
  p1 = d->apex_data()->outer_point() + normal;
  p2 = d->apex_data()->outer_point() - normal;
  path.moveTo(QPointF(p1.x(), p1.y()));
  path.lineTo(QPointF(p2.x(), p2.y()));

  return path;
}

//
//: Return the outline of the vessel so that the convex hull can be computed
//  for collision detection
QPainterPath QGraphicsApexItem::shape() const
{
  Q_D(const QGraphicsApexItem);

  //  We need this hack as QPainterPathStroker will set a width of 1.0
  //  if we pass a value of 0.0 to QPainterPathStroker::setWidth()
  const qreal penWidthZero = qreal(0.00000001);

  // if blank then return blank
  if (path() == QPainterPath())
    return path();

  QPen pen(d->pen());

  QPainterPathStroker ps;
  ps.setCapStyle(pen.capStyle());
  if (pen.widthF() > 0.0)
      ps.setWidth(pen.widthF());
  else
      ps.setWidth(penWidthZero);
  ps.setJoinStyle(pen.joinStyle());
  ps.setMiterLimit(pen.miterLimit());

  // add outline of vessel stroke to painterpath
  QPainterPath p = ps.createStroke(path());
  p.addPath(path());

  return p;
}

/*!
    \reimp
*/
QRectF QGraphicsApexItem::boundingRect() const
{
  Q_D(const QGraphicsApexItem);

  if (d->boundingRect.isNull()) 
  {
    qreal pw = d->pen().widthF();
    if (pw == 0.0)
      d->boundingRect = path().controlPointRect();
    else 
      d->boundingRect = shape().controlPointRect();

    // expand bounding box by 2 units
    d->boundingRect.adjust(-pw, -pw, pw, pw);
  }

  return d->boundingRect;
}


//
//: How the vessel should be drawn in the scene
void QGraphicsApexItem::paint(QPainter *painter, 
                              const QStyleOptionGraphicsItem *option,
                              QWidget *widget)
{
  Q_D(const QGraphicsApexItem);

  Q_UNUSED(widget);

  // Return if apex should not be visible
  if (!appearance().apex_is_visible())
    return;

  painter->setPen(d->pen());
  painter->setBrush(QBrush(Qt::NoBrush));
  painter->drawPath(path());
}

//: Reimplementation
int QGraphicsApexItem::type() const
{
  return Type;
}

bool QGraphicsApexItem::supportsExtension(Extension extension) const
{
  Q_UNUSED(extension);
  return false;
}

void QGraphicsApexItem::setExtension(Extension extension, const QVariant &variant)
{
  Q_UNUSED(extension);
  Q_UNUSED(variant);
}

QVariant QGraphicsApexItem::extension(const QVariant &variant) const
{
  Q_UNUSED(variant);
  return QVariant();
}
