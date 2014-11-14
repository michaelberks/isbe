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
#include "QtGui/private/qgraphicsitem_p.h"

#include <nailfold/ncm_annotation.h>
#include <nailfold/ncm_haemorrhage.h>

#include <nailfold/qtools/ncm_qhaemorrhageitem_appearance.h>

//
//  Private data class
//

class QGraphicsHaemorrhageItemPrivate : public QGraphicsItemPrivate
{
    Q_DECLARE_PUBLIC(QGraphicsHaemorrhageItem)

public:
    //  Methods

    QGraphicsHaemorrhageItemPrivate();

    //: Scene that contains this item
    ncm_qscene* ncm_scene();
    ncm_qscene const* ncm_scene() const;

    //: Get parent vesselitem
    QGraphicsVesselItem* vesselitem();

    //: Haemorrhage that this item represents
    ncm_haemorrhage* haemorrhage_data();
    ncm_haemorrhage const* haemorrhage_data() const;
    void set_haemorrhage_data(ncm_haemorrhage* haemorrhage_data);

    void cache_cursor(const QPointF& cursor_pos);

    QPainterPath placeholder_path() const;
    QPainterPath outline_path() const;

    //: Redraw everything within bounding rectangle
    void updateBoundingRect();

    int opacity() const;
    int line_width() const;
    QBrush brush() const;
    QPen pen() const;

    // status indicators
    bool is_selected() const;
    bool can_be_dragged() const;
    bool can_be_deleted() const;
    bool can_draw_outline() const;
    bool can_delete_outline() const;

    
    //  Variables

    //: Appearance options for *all* apiced (hence it is static)
    static ncm_qhaemorrhageitem_appearance appearance_;

    mutable QRectF boundingRect;

    bool is_hovered_over_;
    bool is_drawing_outline_;

    QPointF cursor_pos_;


private:

    // pointer to the class containing the geometry and properties of the vessel
    // This is private so that you are forced to use vessel_data() that checks
    // the validity of the pointer
    ncm_haemorrhage* haemorrhage_data_;
};

// Define the static member variable
ncm_qhaemorrhageitem_appearance QGraphicsHaemorrhageItemPrivate::appearance_;


QGraphicsHaemorrhageItemPrivate::QGraphicsHaemorrhageItemPrivate()
: is_hovered_over_(false),
  is_drawing_outline_(false)
{
}

//
//: Haemorrhage that this item represents
void QGraphicsHaemorrhageItemPrivate::set_haemorrhage_data(
    ncm_haemorrhage* haemorrhage_data)
{
  haemorrhage_data_ = haemorrhage_data;
}

ncm_haemorrhage* QGraphicsHaemorrhageItemPrivate::haemorrhage_data()
{
  assert(haemorrhage_data_ != 0);
  return haemorrhage_data_;
}

ncm_haemorrhage const* QGraphicsHaemorrhageItemPrivate::haemorrhage_data() const
{
  assert(haemorrhage_data_ != 0);
  return haemorrhage_data_;
}

//
//: Cache mouse cursor position
void QGraphicsHaemorrhageItemPrivate::cache_cursor(const QPointF& cursor_pos)
{
  cursor_pos_ = cursor_pos;
}

//
//: Scene that holds this item
ncm_qscene* QGraphicsHaemorrhageItemPrivate::ncm_scene()
{
  Q_Q(const QGraphicsHaemorrhageItem);
  return dynamic_cast<ncm_qscene*>(q->scene());
}
ncm_qscene const* QGraphicsHaemorrhageItemPrivate::ncm_scene() const
{
  Q_Q(const QGraphicsHaemorrhageItem);
  return dynamic_cast<ncm_qscene const*>(q->scene());
}

//
//: Return opacity (dependent on active/inactive status)
int QGraphicsHaemorrhageItemPrivate::opacity() const
{
  // set opacity
  if (is_selected())
    return appearance_.opacity_selected();
  else
    return appearance_.opacity();
}

//
//: Return line width (dependent on vessel properties)
int QGraphicsHaemorrhageItemPrivate::line_width() const
{
  return appearance_.line_width();
}

//
//: Set colour and opacity for arterial limb
QBrush QGraphicsHaemorrhageItemPrivate::brush() const
{
  // set colour
  QColor rgba;
  if (is_selected())
    rgba.setRgb(appearance_.colour_selected());
  else
    rgba.setRgb(appearance_.colour());

  rgba.setAlpha(opacity());

  return QBrush(rgba);
}

//
//: Set additional pen properties (i.e. width) for arterial limb
QPen QGraphicsHaemorrhageItemPrivate::pen() const
{
  return QPen(brush(), line_width());
}

//
//: Path for the placeholder only
QPainterPath QGraphicsHaemorrhageItemPrivate::placeholder_path() const
{
  QPainterPath path;
  const double placeholder_radius = appearance_.placeholder_radius() * 
                                    ncm_scene()->imageRect().height();

  // represent vessel by a circle centred on the only marked point
  QPointF centre(haemorrhage_data()->anchor().x(),
                 haemorrhage_data()->anchor().y());
  path.addRect(centre.x() - placeholder_radius, centre.y() - placeholder_radius,
               placeholder_radius*2, placeholder_radius*2);

  return path;
}

//
//: Path for the vessel centreline
QPainterPath QGraphicsHaemorrhageItemPrivate::outline_path() const
{
  //unsigned n_points = haemorrhage_data()->n_points();

  QPainterPath path;
  //if (n_points > 0)
  //{
  //  path.moveTo(QPointF(vessel_data()->point(0)->x(),
  //                      vessel_data()->point(0)->y()));
  //  for (unsigned i = 1; i < n_points; ++i)
  //  {
  //    path.lineTo(QPointF(vessel_data()->point(i)->x(),
  //                        vessel_data()->point(i)->y()));
  //  }
  //}

  return path;
}

//
//: Redraw everything within the bounding rectangle
void QGraphicsHaemorrhageItemPrivate::updateBoundingRect() 
{
  Q_Q(QGraphicsHaemorrhageItem);

  // recompute bounding rectangle and redraw
  q->prepareGeometryChange();
  boundingRect = QRect();
  q->update(q->boundingRect());
}


//: Status indicators

//
//: Is vesselitem selected in list?
bool QGraphicsHaemorrhageItemPrivate::is_selected() const
{
  //return (vessel_data() == ncm_scene()->selected_vessel());
  return false;
}

bool QGraphicsHaemorrhageItemPrivate::can_be_dragged() const
{
  return false;
}

bool QGraphicsHaemorrhageItemPrivate::can_be_deleted() const
{
  return (ncm_scene()->edit_mode() == ncm_qscene::ModeAddHaemorrhages &&
          is_hovered_over_);
}

bool QGraphicsHaemorrhageItemPrivate::can_draw_outline() const
{
  //return (
  //        ncm_scene()->edit_mode() == ncm_qscene::ModeVesselPath &&
  //        !(d->vessel_data()->path_is_defined()) && 
  //        is_selected()
  //       );
  return false;
}

bool QGraphicsHaemorrhageItemPrivate::can_delete_outline() const
{
  //return (
  //        ncm_scene()->edit_mode() == ncm_qscene::ModeVesselPath &&
  //        d->vessel_data()->path_is_defined() && 
  //        (is_selected() || d->is_hovered_over_)
  //       );
  return false;
}
