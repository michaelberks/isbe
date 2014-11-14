//:
// \file
// \brief Visualization of a blood vessel. This is based on the 
//        QGraphicsPathItem, hence why it doesn't quite fit in the style of the
//        other code
// \author Phil Tresadern

#include "ncm_qvesselitem.h"
#include "ncm_qscene.h"
#include "ncm_qapexitem.h"
#include "ncm_qglobal.h"

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
//#include "QtGui/private/qgraphicsitem_p.h"

#include <nailfold/ncm_annotation.h>

// This uses the source (rather than header) file so that QGraphicsVesselItem's
// declaration can remain hidden
#include "ncm_qvesselitem_p.cxx"

//
//  Static member functions
//

//
//: Access to appearance properties
/* static */ ncm_qvesselitem_appearance& QGraphicsVesselItem::appearance()
{
  return QGraphicsVesselItemPrivate::appearance_;
}


/*!
    Constructs a QGraphicsPath. \a parent is passed to
    QAbstractGraphicsShapeItem's constructor.

    \sa QGraphicsScene::addItem()
*/
QGraphicsVesselItem::QGraphicsVesselItem(QGraphicsItem *parent
#ifndef Q_QDOC
                                     // obsolete argument
                                     , QGraphicsScene *scene
#endif
)
: QGraphicsItem(*new QGraphicsVesselItemPrivate, parent, scene)
{
}

/*!
    Constructs a QGraphicsPath item using \a path as the default path. \a
    parent is passed to QAbstractGraphicsShapeItem's constructor.

    \sa QGraphicsScene::addItem()
*/
QGraphicsVesselItem::QGraphicsVesselItem(ncm_vessel* const vessel_data,
                                         QGraphicsItem *parent
#ifndef Q_QDOC
                                     // obsolete argument
                                     , QGraphicsScene *scene
#endif
)
: QGraphicsItem(*new QGraphicsVesselItemPrivate, parent, scene)
{
  // Although scene is marked as an obsolete argument, it's still used (and 
  // useful) here because the apexitems have signals that are connected to the
  // scene which needs to be known at construction (which is not the case if
  // the vesselitem is added to the scene only after construction). This may
  // need changing if a future version of Qt deletes the scene argument 
  // altogether.

  Q_D(QGraphicsVesselItem);

  d->set_vessel_data(vessel_data);

  // call hoverMoveEvent when the mouse hovers over the vessel
  setAcceptHoverEvents(true);
}

/*!
    Destroys the QGraphicsVesselItem.
*/
QGraphicsVesselItem::~QGraphicsVesselItem()
{
  Q_D(QGraphicsVesselItem);

  assert( d->ncm_scene() != NULL );
  assert( d->ncm_scene()->annotation() != NULL );

  // Repeatedly remove the first apexitem in the list until it is empty
  while (!d->ApexItems_.empty())
    delete *(d->ApexItems_.begin());

  // Remove the vessel entry from the corresponding markup
  d->ncm_scene()->annotation()->delete_vessel(d->vessel_data());

  // Detach this vesselitem from the scene before it is destroyed
  d->ncm_scene()->removeItem(this);
}

void QGraphicsVesselItem::remove_apexitem(const QGraphicsApexItem* apexitem)
{
  Q_D(QGraphicsVesselItem);

  vcl_vector<QGraphicsApexItem*>::iterator it = 
      vcl_find(d->ApexItems_.begin(), d->ApexItems_.end(), apexitem);

  if (it != d->ApexItems_.end())
    d->ApexItems_.erase(it);
}

//
//: Whether this item matches the given vessel data
bool QGraphicsVesselItem::is_representing(const ncm_vessel* const vessel) const
{
  Q_D(const QGraphicsVesselItem);

  if (vessel != NULL)
    return (d->vessel_data() == vessel);
  else
    return false;
}

//
//: Is vesselitem selected in list?
bool QGraphicsVesselItem::is_selected() const
{
  Q_D(const QGraphicsVesselItem);
  return d->is_selected();
}

bool QGraphicsVesselItem::is_adding_points() const
{
  Q_D(const QGraphicsVesselItem);
  return (d->is_adding_venous_points() ||
          d->is_adding_arterial_points());
}

//
//: Redraw everything within the bounding rectangle
void QGraphicsVesselItem::updateBoundingRect() 
{
  Q_D(QGraphicsVesselItem);

  // recompute bounding rectangle and redraw
  prepareGeometryChange();
  d->boundingRect = QRect();
  update(boundingRect());
}

//
//: What to do when the mouse is moved while no button is pressed
void QGraphicsVesselItem::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
//#if _DEBUG
//  vcl_cout << __FUNCTION__ << vcl_endl;
//#endif

  Q_D(QGraphicsVesselItem);
  d->is_hovered_over_ = true;
  updateBoundingRect();

  // update snapped points for later use
  d->update_nearest_points_to(event->scenePos());

  // Keep tips up to date
  vcl_stringstream index_stream;
  index_stream << d->vessel_data()->index()+1; // 1,2,3,... instead of 0,1,2,...
  vcl_string index_string = "Vessel: " + index_stream.str();
  vcl_string size_string = "Size: " + d->vessel_data()->properties().size_string();
  vcl_string shape_string = "Shape: " + d->vessel_data()->properties().shape_string();
  vcl_string tooltip_string = index_string + '\n' + 
                              size_string + '\n' + 
                              shape_string;
  setToolTip(tooltip_string.c_str());
}

void QGraphicsVesselItem::hoverMoveEvent(QGraphicsSceneHoverEvent *event)
{
//#if _DEBUG
//  vcl_cout << __FUNCTION__ << vcl_endl;
//#endif

  Q_D(QGraphicsVesselItem);
  d->is_hovered_over_ = true;
  updateBoundingRect();

  // update snapped points for later use
  d->update_nearest_points_to(event->scenePos());
}

void QGraphicsVesselItem::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
//#if _DEBUG
//  vcl_cout << __FUNCTION__ << vcl_endl;
//#endif

  Q_D(QGraphicsVesselItem);
  d->is_hovered_over_ = false;
  updateBoundingRect();

  // update snapped points for later use
  d->update_nearest_points_to(event->scenePos());
}

//
//: What to do when a mouse button is pressed
void QGraphicsVesselItem::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
//#if _DEBUG
//  vcl_cout << __FUNCTION__ << vcl_endl;
//#endif

  // We can arrive here through two routes:
  // - The user has clicked on an existing vesselitem (e.g. to add an apex)
  // - The user has clicked on empty space while this vesselitem was selected
  //   in the list

  // This may trigger one of the following actions
  // - Set vessel properties
  // - Add vessel path
  // - Start appending arterial/venous points
  // - Create an apex
  // - Delete vessel
  // - Delete vessel path
  // - Trim vessel
  // - Associate a haemorrhage with a vessel

  Q_D(QGraphicsVesselItem);

  // unset accepted flag
  event->ignore();

  // update snapped points for later use
  d->update_nearest_points_to(event->scenePos());

  const bool shift_is_pressed = (event->modifiers() & Qt::ShiftModifier);

  switch (event->button())
  {
    case Qt::LeftButton:
    {
      if (d->can_set_size())
      {
        // Warn if vessel is already labelled
        if (  appearance().warn_on_overwrite_size &&
            !(d->vessel_data()->properties().is_size_undefined()) )
        {
          QMessageBox msgBox;
          msgBox.setIcon(QMessageBox::Warning);
          msgBox.setText("Vessel size already defined");
          msgBox.setInformativeText("Overwrite?");
          msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
          msgBox.setDefaultButton(QMessageBox::No);

          const int response = msgBox.exec();
          if (response == QMessageBox::No)
            break;
        }

        // copy properties from somewhere to this vessel
        ncm_vessel_properties& p = d->vessel_data()->properties();
        p.copy_size_from(d->ncm_scene()->new_vessel_properties());

        // A QGraphicsItem cannot emit signals of its own so we must tell the
        // scene to do it on our behalf
        d->ncm_scene()->emit_vessel_changed();
        event->accept();
        break;
      }
      else if (d->can_set_shape())
      {
        // Warn if vessel is already labelled
        if (  appearance().warn_on_overwrite_shape &&
            !(d->vessel_data()->properties().is_shape_undefined()) )
        {
          QMessageBox msgBox;
          msgBox.setIcon(QMessageBox::Warning);
          msgBox.setText("Vessel shape already defined");
          msgBox.setInformativeText("Overwrite?");
          msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
          msgBox.setDefaultButton(QMessageBox::No);

          const int response = msgBox.exec();
          if (response == QMessageBox::No)
            break;
        }

        // copy properties from somewhere to this vessel
        ncm_vessel_properties& p = d->vessel_data()->properties();
        p.copy_shape_from(d->ncm_scene()->new_vessel_properties());

        // A QGraphicsItem cannot emit signals of its own so we must tell the
        // scene to do it on our behalf
        d->ncm_scene()->emit_vessel_changed();
        event->accept();
        break;
      }
      else if (d->can_add_path())
      {
        // add first point to the path
        const QPointF pos = d->cursor_pos_; // see note in MouseMoveEvent below
        d->vessel_data()->add_venous_point(pos.x(), pos.y());
        d->is_adding_venous_points_ = true;
        event->accept();
        break;
      }
      else if (d->can_add_venous_points())
      {
        // add venous points
        d->is_adding_venous_points_ = true;
        event->accept();
        break;
      }
      else if (d->can_add_arterial_points())
      {
        // add arterial points
        d->is_adding_arterial_points_ = true;
        event->accept();
        break;
      }
      else if (d->can_add_apex())
      {
        // Create a new apex at given location
        d->create_apexitem_at(d->cursor_pos_);
        event->accept();
        break;
      }

      break;
    }

    case Qt::RightButton:
    {
      if (d->can_be_deleted())
      {
        // Shift not pressed => delete whole vessel
        // In this case, the vesselitem deletes itself and nothing
        // should try to access its members (data or functions). Therefore, we
        // take a copy of the pointer to the scene so that we can emit the
        // vessel_deleted() signal after the vessel has actually been deleted.
        // Otherwise, the list is updated first, therefore containing the now
        // deleted vessel.
        ncm_qscene* scene = d->ncm_scene();
        delete this;
        scene->emit_vessel_deleted();
        event->accept();
        return;
      }
      else if (d->can_delete_path() && !shift_is_pressed)
      {
        // delete path only
        d->vessel_data()->delete_path();
        event->accept();
        break;
      }
      else if (d->can_trim_path() && shift_is_pressed)
      {
        d->vessel_data()->trim_at(d->cursor_pos_.x(),
                                  d->cursor_pos_.y());
        event->accept();
        break;
      }

      break;
    }

    default:
    {}
  }

  updateBoundingRect();
}

void QGraphicsVesselItem::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
//#if _DEBUG
//  vcl_cout << __FUNCTION__ << vcl_endl;
//#endif

  Q_D(QGraphicsVesselItem);

  // unset accepted flag
  event->ignore();

  // update snapped points for later use
  d->update_nearest_points_to(event->scenePos());

  // copy new point to local variable
  QPointF pos;
  if (d->ncm_scene()->is_snapping())
    pos = d->nearest_line_point_;
  else
    pos = d->cursor_pos_;

  // Note: While writing this code I noticed a very peculiar 'feature' in debug
  //       mode. If you try to add d->cursor_pos_ directly, i.e.,
  //
  //         add_venous_point(d->cursor_pos_.x(), d->cursor_pos_.y())
  //
  //       then it has an adverse effect on the refresh rate when tracing
  //       out a line. I never figured out why exactly (it isn't just because
  //       you're accessing d->cursor_pos_ twice), but copying to a local
  //       variable is neater anyway. Very odd.

  if (d->is_adding_venous_points())
  {
    // add venous points
    d->vessel_data()->add_venous_point(pos.x(), pos.y());
    event->accept();
  }
  else if (d->is_adding_arterial_points())
  {
    // add arterial points
    d->vessel_data()->add_arterial_point(pos.x(), pos.y());
    event->accept();
  }

  updateBoundingRect();
}

void QGraphicsVesselItem::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
//#if _DEBUG
//  vcl_cout << __FUNCTION__ << vcl_endl;
//#endif

  Q_D(QGraphicsVesselItem);

  // unset accepted flag
  event->ignore();

  // update snapped points for later use
  d->update_nearest_points_to(event->scenePos());

  if (is_adding_points() && 
      d->vessel_data()->path_is_suspicious())
  {
    // Create a warning dialog
    QMessageBox msgBox;
    msgBox.setText("Path is closer to an unselected vessel");
    msgBox.setInformativeText("Keep path anyway?");
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::No);

    // Show the dialog and react accordingly
    if (msgBox.exec() == QMessageBox::No)
      d->vessel_data()->delete_path();

    // NOTE: this last statement will flag the annotation as modified when that
    // may not necessarily be the case (at least as far as we are concerned)
    // since the path will have been added, then deleted immediately.
    // In itself, this isn't a major problem - just worth being aware of.
  }

  if (d->is_adding_venous_points())
  {
    d->is_adding_venous_points_ = false;
    d->vessel_data()->set_modified();
    event->accept();
  }
  else if (d->is_adding_arterial_points())
  {
    d->is_adding_arterial_points_ = false;
    d->vessel_data()->set_modified();
    event->accept();
  }

  updateBoundingRect();

  // because we sometimes grab the mouse explicitly (i.e. when creating a new
  // vessel) we have to explicitly ungrab it, too.
  ungrabMouse();
}

/*!
    Returns the item's path as a QPainterPath. If no item has been set, an
    empty QPainterPath is returned.

    \sa setPath()
*/

QPainterPath QGraphicsVesselItem::path() const
{
  Q_D(const QGraphicsVesselItem);

  QPainterPath path;
  if (appearance().placeholder_is_visible() ||
      (!d->vessel_data()->path_is_defined()))
    path.addPath(d->placeholder_path());

  if (appearance().path_is_visible() &&
      d->vessel_data()->path_is_defined())
    path.addPath(d->vessel_path());

  return path;
}

//
//: Return the outline of the vessel so that the convex hull can be computed
//  for collision detection
QPainterPath QGraphicsVesselItem::shape() const
{
  Q_D(const QGraphicsVesselItem);

  //  We need this hack as QPainterPathStroker will set a width of 1.0
  //  if we pass a value of 0.0 to QPainterPathStroker::setWidth()
  const qreal penWidthZero = qreal(0.00000001);

  // if blank then return blank
  if (path() == QPainterPath())
      return path();

  QPen pen; // Default pen settings

  QPainterPath p;
  QPainterPathStroker ps;
  ps.setCapStyle(pen.capStyle());
  ps.setJoinStyle(pen.joinStyle());
  ps.setMiterLimit(pen.miterLimit());

  // Add placeholder path (if needed)
  if (appearance().placeholder_is_visible() ||
      (!d->vessel_data()->path_is_defined()))
  {
    p.addEllipse(d->placeholder_path().boundingRect());
  }

  // Add vessel path (if needed)
  if (appearance().path_is_visible() &&
      d->vessel_data()->path_is_defined())
  {
    pen = d->path_pen();

    if (pen.widthF() > 0.0)
    {
      ps.setWidth(pen.widthF());
      p.addPath(ps.createStroke(d->vessel_path())); // add outline
    }
    else
    {
      //ps.setWidth(penWidthZero);
      p.addPath(d->vessel_path()); // add centreline
    }

    // add ellipse(s) representing endpoint to painterpath
    if (d->is_adding_venous_points_ ||
        d->can_add_venous_points())
    {
      QPointF centre(d->vessel_data()->venous_endpoint()->x(),
                     d->vessel_data()->venous_endpoint()->y());
      p.addEllipse(centre, d->line_width(), d->line_width());
    }
    else if (d->is_adding_arterial_points_ ||
             d->can_add_arterial_points())
    {
      QPointF centre(d->vessel_data()->arterial_endpoint()->x(),
                     d->vessel_data()->arterial_endpoint()->y());
      p.addEllipse(centre, d->line_width(), d->line_width());
    }
  }

  return p;
}

/*!
    \reimp
*/
QRectF QGraphicsVesselItem::boundingRect() const
{
  Q_D(const QGraphicsVesselItem);

  if (d->boundingRect.isNull()) {
    //const qreal pen_width = d->pen().widthF();
    //const qreal pen_width = appearance().line_width();
    const qreal pen_width = d->line_width();

    //if (pen_width == 0.0)
    //  d->boundingRect = path().controlPointRect();
    //else
      d->boundingRect = shape().controlPointRect();

    // Expand bounding box by the width of the pen that is used to define
    // shape()
    d->boundingRect.adjust(-pen_width, -pen_width, 
                           pen_width, pen_width);
  }

  return d->boundingRect;
}


//
//: How the vessel should be drawn in the scene
void QGraphicsVesselItem::paint(QPainter *painter, 
                                const QStyleOptionGraphicsItem *option,
                                QWidget *widget)
{
  Q_D(const QGraphicsVesselItem);

  Q_UNUSED(widget);

  // Don't draw anything at all if the vessel is non-distal, unless we're in 
  // Tag Vessels mode
  if (!(d->vessel_data()->properties().is_distal()) &&
      d->ncm_scene()->edit_mode() != ncm_qscene::ModeAddVessels &&
			d->ncm_scene()->edit_mode() != ncm_qscene::ModeDisplayAll)
  {
    return;
  }
  
  painter->setBrush(QBrush(Qt::NoBrush));

  const bool path_is_expected = appearance().path_is_visible();
  const bool path_is_defined = (d->vessel_data()->path_is_defined());

  // Draw placeholder if specified or in the absence of a path that should be
  // shown.
  if (appearance().placeholder_is_visible() || 
      (path_is_expected && !path_is_defined))
  {
    painter->setPen(d->placeholder_pen());
    painter->drawPath(d->placeholder_path());
  }

  // Draw the path if it should be shown and has been defined.
  if (path_is_expected && path_is_defined)
  {
    painter->setPen(d->path_pen());
    painter->drawPath(d->vessel_path());

    // draw a blob at each end if the vessel is active to enable the user to add
    // more points to either end
    painter->setPen(QPen(Qt::NoPen));

    if (d->is_adding_venous_points_ ||
        d->can_add_venous_points())
    {
      painter->setBrush(d->brush());
      QPointF centre(d->vessel_data()->venous_endpoint()->x(),
                     d->vessel_data()->venous_endpoint()->y());
      painter->drawEllipse(centre, d->line_width(), d->line_width());
    }
    else if (d->is_adding_arterial_points_ ||
             d->can_add_arterial_points())
    {
      // draw arterial path
      painter->setBrush(d->brush());
      QPointF centre(d->vessel_data()->arterial_endpoint()->x(),
                     d->vessel_data()->arterial_endpoint()->y());
      painter->drawEllipse(centre, d->line_width(), d->line_width());
    }
  }
}

//: Reimplementation
int QGraphicsVesselItem::type() const
{
    return Type;
}

bool QGraphicsVesselItem::supportsExtension(Extension extension) const
{
    Q_UNUSED(extension);
    return false;
}

void QGraphicsVesselItem::setExtension(Extension extension, const QVariant &variant)
{
    Q_UNUSED(extension);
    Q_UNUSED(variant);
}

QVariant QGraphicsVesselItem::extension(const QVariant &variant) const
{
    Q_UNUSED(variant);
    return QVariant();
}


