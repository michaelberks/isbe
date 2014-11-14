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
#include "QtGui/private/qgraphicsitem_p.h"

#include <nailfold/ncm_annotation.h>
#include <nailfold/ncm_apex.h>

#include <nailfold/qtools/ncm_qapexitem_appearance.h>
#include <nailfold/qtools/ncm_qvesselitem_appearance.h>

//
//  Private data class
//

class QGraphicsApexItemPrivate : public QGraphicsItemPrivate
{
    Q_DECLARE_PUBLIC(QGraphicsApexItem)

public:
    //  Methods

    QGraphicsApexItemPrivate();

    void updateBoundingRect();

    //: Scene that contains this item
    ncm_qscene* ncm_scene();
    ncm_qscene const* ncm_scene() const;

    //: Parent vesselitem
    QGraphicsVesselItem* vesselitem();
    QGraphicsVesselItem const* vesselitem() const;

    //: Apex that this item represents
    void set_apex_data(ncm_apex* apex_data);
    ncm_apex* apex_data();
    ncm_apex const* apex_data() const;

    //: Appearance parameters specific to this vessel
    int opacity() const;
    int line_width() const;
    QBrush brush() const;
    QPen pen() const;

    //: Status indicators
    bool can_be_deleted() const;
    bool can_move_inner_point() const;
    bool can_move_outer_point() const;


    //  Variables

    //: Appearance options for *all* apiced (hence it is static)
    static ncm_qapexitem_appearance appearance_;

    mutable QRectF boundingRect;

    bool is_hovered_over_;
    bool is_moving_inner_point_;
    bool is_moving_outer_point_;

    void cache_cursor(const QPointF& cursor_pos);
    QPointF cursor_pos_;

private:

    bool can_move_point(vgl_point_2d<double> point) const;

    // pointer to the class containing the geometry and properties of the vessel
    // This is private so that you are forced to use vessel_data() that checks
    // the validity of the pointer
    ncm_apex* apex_data_;
};

// Define the static member variable
ncm_qapexitem_appearance QGraphicsApexItemPrivate::appearance_;


QGraphicsApexItemPrivate::QGraphicsApexItemPrivate()
: is_hovered_over_(false),
  is_moving_inner_point_(false),
  is_moving_outer_point_(false)
{
}

//
//: The scene that holds the apexitem
ncm_qscene* QGraphicsApexItemPrivate::ncm_scene()
{
  Q_Q(QGraphicsApexItem);
  return dynamic_cast<ncm_qscene*>(q->scene());
}
ncm_qscene const* QGraphicsApexItemPrivate::ncm_scene() const
{
  Q_Q(const QGraphicsApexItem);
  return dynamic_cast<ncm_qscene const*>(q->scene());
}

//
//: The vesselitem to which the apexitem belongs
QGraphicsVesselItem* QGraphicsApexItemPrivate::vesselitem()
{
  return qgraphicsitem_cast<QGraphicsVesselItem*>(parent);
}

QGraphicsVesselItem const* QGraphicsApexItemPrivate::vesselitem() const
{
  return qgraphicsitem_cast<QGraphicsVesselItem const*>(parent);
}

//
//: The ncm_apex that this item represents
void QGraphicsApexItemPrivate::set_apex_data(ncm_apex* apex_data)
{
  apex_data_ = apex_data;
}
ncm_apex* QGraphicsApexItemPrivate::apex_data()
{
  assert(apex_data_ != 0);
  return apex_data_;
}
ncm_apex const* QGraphicsApexItemPrivate::apex_data() const
{
  assert(apex_data_ != 0);
  return apex_data_;
}

//
//: Cache mouse cursor position
void QGraphicsApexItemPrivate::cache_cursor(const QPointF& cursor_pos)
{
  cursor_pos_ = cursor_pos;
}

//
//: Return opacity (dependent on active/inactive status)
int QGraphicsApexItemPrivate::opacity() const
{
  // set opacity
  //if (d->edit_mode_ != QGraphicsApexItemPrivate::ModeInactive)
    return appearance_.opacity_selected();
  //else
  //  return d->opacity_inactive_;
}

//
//: Return line width (dependent on vessel properties)
int QGraphicsApexItemPrivate::line_width() const
{
  int baseline = appearance_.line_width();

  if (appearance_.width_is_relative_to_image_size())
    baseline *= ncm_scene()->imageRect().height();

  // Modulate line width based on thickness of the vessel
  if (apex_data()->parent_vessel().properties().is_size_enlarged())
    return baseline * vesselitem()->appearance().scale_enlarged();
  if (apex_data()->parent_vessel().properties().is_size_giant())
    return baseline * vesselitem()->appearance().scale_enlarged()
                    * vesselitem()->appearance().scale_giant();
  else
    return baseline;
}

//
//: Set colour and opacity for arterial limb
QBrush QGraphicsApexItemPrivate::brush() const
{
  // set colour
  QColor rgba;
  
  if (can_move_inner_point() || can_move_outer_point())
    rgba = appearance_.colour_can_move();
  else if (can_be_deleted())
    rgba = appearance_.colour_can_delete();
  else
    rgba = appearance_.colour();

  rgba.setAlpha(opacity());

  return QBrush(rgba);
}

//
//: Set additional pen properties (i.e. width) for arterial limb
QPen QGraphicsApexItemPrivate::pen() const
{
  return QPen(brush(), line_width());
}

//
//: Redraw everything within the bounding rectangle
void QGraphicsApexItemPrivate::updateBoundingRect() 
{
  Q_Q(QGraphicsApexItem);

  // recompute bounding rectangle and redraw
  q->prepareGeometryChange();
  boundingRect = QRect();
  q->update(q->boundingRect());
}

//: Status indicators
bool QGraphicsApexItemPrivate::can_be_deleted() const
{
  return (ncm_scene()->edit_mode() == ncm_qscene::ModeLabelApices &&
          is_hovered_over_ &&
          vesselitem()->is_selected());
}

bool QGraphicsApexItemPrivate::can_move_inner_point() const
{
  return can_move_point(apex_data()->inner_point());
}

bool QGraphicsApexItemPrivate::can_move_outer_point() const
{
  return can_move_point(apex_data()->outer_point());
}

//
//  Private methods
//

bool QGraphicsApexItemPrivate::can_move_point(vgl_point_2d<double> point) const
{
  // must be in right mode
  const bool correct_mode = 
      (ncm_scene()->edit_mode() == ncm_qscene::ModeLabelApices);

  if (!correct_mode || !is_hovered_over_)
    return false;

  // get current mouse location as a vgl_point_2d
  vgl_point_2d<double> cursor_pos_vgl(cursor_pos_.x(), cursor_pos_.y());

  // if very close to the inner endpoint
  const vgl_vector_2d<double> pos_to_point = point - cursor_pos_vgl;

  const double dist_to_point = pos_to_point.length();

  // TODO: Create a user-defined variable to set this sensitivity
  const double threshold = 0.1 * apex_data()->width();
  return (dist_to_point < threshold);
}
