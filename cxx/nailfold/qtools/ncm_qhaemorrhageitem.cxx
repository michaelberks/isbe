//:
// \file
// \brief QGraphicsItem-derived class for creating a visualization of the two
//        points either side of the apex (during annotation).
// \author Phil Tresadern

#include "ncm_qhaemorrhageitem.h"
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
#include <nailfold/ncm_haemorrhage.h>

// This uses the source (rather than header) file so that QGraphicsApexItem's
// interface can remain hidden
#include "ncm_qhaemorrhageitem_p.cxx"

//
//: Access to appearance properties
/* static */ ncm_qhaemorrhageitem_appearance& QGraphicsHaemorrhageItem::appearance()
{
  return QGraphicsHaemorrhageItemPrivate::appearance_;
}

/*!
    Constructs a QGraphicsPath. \a parent is passed to
    QAbstractGraphicsShapeItem's constructor.

    \sa QGraphicsScene::addItem()
*/
QGraphicsHaemorrhageItem::QGraphicsHaemorrhageItem(
    ncm_haemorrhage* const haemorrhage_data,
    QGraphicsItem *parent /* = 0 */
#ifndef Q_QDOC
    // obsolete argument
    , QGraphicsScene *scene
#endif
)
: QGraphicsItem(*new QGraphicsHaemorrhageItemPrivate, parent, scene)
{
  Q_D(QGraphicsHaemorrhageItem);

  // FIXME: we should check that vessel_data is not zero
  d->set_haemorrhage_data(haemorrhage_data);
  
  // call hoverMoveEvent when the mouse hovers over the vesselnap 
  setAcceptHoverEvents(true);
}

//
//: Destructor
//  Needs to take care of deleting the apex information from the markup, too
QGraphicsHaemorrhageItem::~QGraphicsHaemorrhageItem()
{
  Q_D(QGraphicsHaemorrhageItem);

  assert( d->ncm_scene() != NULL );
  assert( d->ncm_scene()->annotation() != NULL );

  // Remove the vessel entry from the corresponding markup
  d->ncm_scene()->annotation()->delete_haemorrhage(d->haemorrhage_data());

  // Detach this vesselitem from the scene before it is destroyed
  d->ncm_scene()->removeItem(this);
}

//
//: What to do when the mouse is moved while no button is pressed
void QGraphicsHaemorrhageItem::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
  Q_D(QGraphicsHaemorrhageItem);
  d->cache_cursor(event->scenePos());
  d->is_hovered_over_ = true;
  d->updateBoundingRect();
}

void QGraphicsHaemorrhageItem::hoverMoveEvent(QGraphicsSceneHoverEvent *event)
{
  Q_D(QGraphicsHaemorrhageItem);
  d->cache_cursor(event->scenePos());
  d->is_hovered_over_ = true;
  d->updateBoundingRect();
}

void QGraphicsHaemorrhageItem::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
  Q_D(QGraphicsHaemorrhageItem);
  d->cache_cursor(event->scenePos());
  d->is_hovered_over_ = false;
  d->updateBoundingRect();
}

//
//: What to do when a mouse button is pressed
void QGraphicsHaemorrhageItem::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  Q_D(QGraphicsHaemorrhageItem);

  // unset accepted flag
  event->ignore();

  d->cache_cursor(event->scenePos());

  switch (event->button())
  {
    case Qt::LeftButton:
    {
      if (d->can_be_dragged())
      {
      }
      else if (d->can_draw_outline())
      {
      }

      break;
    }

    case Qt::RightButton:
    {
      if (d->can_delete_outline())
      {
        // delete the boundary definition
      }
      if (d->can_be_deleted())
      {
        delete this;
        // emit signal from scene?
        event->accept();
        return;
      }
    }
  }

  d->updateBoundingRect();
}

//
//: What to do when the mouse is moved while a mouse button is pressed
void QGraphicsHaemorrhageItem::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  Q_D(QGraphicsHaemorrhageItem);

  // unset accepted flag
  event->ignore();

  d->cache_cursor(event->scenePos());

  if (d->is_drawing_outline_)
  {
    //d->apex_data()->set_inner_point(d->cursor_pos_.x(), 
    //                                d->cursor_pos_.y());
    //event->accept();
  }

  d->updateBoundingRect();
}

//
//: What to do when a mouse button is released
void QGraphicsHaemorrhageItem::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
  Q_D(QGraphicsHaemorrhageItem);

  // unset accepted flag
  event->ignore();

  d->cache_cursor(event->scenePos());

  if (d->is_drawing_outline_)
  {
    //d->apex_data()->set_inner_point(d->cursor_pos_.x(), 
    //                                d->cursor_pos_.y());
    //d->is_moving_inner_point_ = false;
    //event->accept();
  }

  d->updateBoundingRect();
}

/*!
    Returns the item's path as a QPainterPath. If no item has been set, an
    empty QPainterPath is returned.

    \sa setPath()
*/

QPainterPath QGraphicsHaemorrhageItem::path() const
{
  Q_D(const QGraphicsHaemorrhageItem);

  QPainterPath path;

  if (appearance().placeholder_is_visible() ||
      (!d->haemorrhage_data()->outline_is_defined()))
    path.addPath(d->placeholder_path());

  if (appearance().outline_is_visible() &&
      d->haemorrhage_data()->outline_is_defined())
    path.addPath(d->outline_path());

  return path;
}

//
//: Return the outline of the vessel so that the convex hull can be computed
//  for collision detection
QPainterPath QGraphicsHaemorrhageItem::shape() const
{
  Q_D(const QGraphicsHaemorrhageItem);

  //  We need this hack as QPainterPathStroker will set a width of 1.0
  //  if we pass a value of 0.0 to QPainterPathStroker::setWidth()
  const qreal penWidthZero = qreal(0.00000001);

  // if blank then return blank
  if (path() == QPainterPath())
      return path();

  QPen pen(d->pen());

  QPainterPath p;
  QPainterPathStroker ps;
  ps.setCapStyle(pen.capStyle());
  ps.setJoinStyle(pen.joinStyle());
  ps.setMiterLimit(pen.miterLimit());


  // Add placeholder path (if needed)
  if (appearance().placeholder_is_visible() ||
      (!d->haemorrhage_data()->outline_is_defined()))
  {
    p.addEllipse(d->placeholder_path().boundingRect());
  }

  // Add vessel path (if needed)
  if (appearance().outline_is_visible() &&
      d->haemorrhage_data()->outline_is_defined())
  {
    pen = d->pen();

    if (pen.widthF() > 0.0)
    {
      ps.setWidth(pen.widthF());
      p.addPath(ps.createStroke(d->outline_path())); // add outline
    }
    else
      p.addPath(d->outline_path()); // add centreline
  }

  return p;
}

/*!
    \reimp
*/
QRectF QGraphicsHaemorrhageItem::boundingRect() const
{
    Q_D(const QGraphicsHaemorrhageItem);

    if (d->boundingRect.isNull()) {
        qreal pw = d->pen().widthF();
        if (pw == 0.0)
            d->boundingRect = path().controlPointRect();
        else {
            d->boundingRect = shape().controlPointRect();
        }
    }

    // Expand bounding box by the width of the pen that is used to define
    // shape()
    const double lw = d->line_width();
    d->boundingRect.adjust(-lw, -lw, lw, lw);

    return d->boundingRect;
}


//
//: How the vessel should be drawn in the scene
void QGraphicsHaemorrhageItem::paint(QPainter *painter, 
                                const QStyleOptionGraphicsItem *option,
                                QWidget *widget)
{
  Q_D(const QGraphicsHaemorrhageItem);

  Q_UNUSED(widget);
 
  painter->setPen(d->pen());
  painter->setBrush(QBrush(Qt::NoBrush));
  painter->drawPath(path());
}

//: Reimplementation
int QGraphicsHaemorrhageItem::type() const
{
  return Type;
}

bool QGraphicsHaemorrhageItem::supportsExtension(Extension extension) const
{
  Q_UNUSED(extension);
  return false;
}

void QGraphicsHaemorrhageItem::setExtension(Extension extension, const QVariant &variant)
{
  Q_UNUSED(extension);
  Q_UNUSED(variant);
}

QVariant QGraphicsHaemorrhageItem::extension(const QVariant &variant) const
{
  Q_UNUSED(variant);
  return QVariant();
}


