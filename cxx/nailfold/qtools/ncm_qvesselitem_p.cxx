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

#include "QtGui/private/qgraphicsitem_p.h"

#include <nailfold/ncm_annotation.h>
#include <nailfold/qtools/ncm_qvesselitem_appearance.h>

//
//  Private data class
//

class QGraphicsVesselItemPrivate : public QGraphicsItemPrivate
{
    Q_DECLARE_PUBLIC(QGraphicsVesselItem)

public:
    QGraphicsVesselItemPrivate();
    ~QGraphicsVesselItemPrivate();

    //: Return pointer to ncm_scene (as opposed to default scene() that returns 
    //  a QGraphicsScene
    ncm_qscene* ncm_scene();
    ncm_qscene const* ncm_scene() const;

    //: Set/get vessel data via a pointer
    void set_vessel_data(ncm_vessel* const vessel_data);
    ncm_vessel* vessel_data();
    ncm_vessel const* vessel_data() const;

    //: Create apexitems
    QGraphicsApexItem* create_apexitem_from(ncm_apex* apex, 
                                            bool interactive = true);
    void create_apexitem_at(const QPointF& first_point);

    //: Get appearance properties for this specific vessel as shown onscreen
    QColor colour() const;
    int opacity() const;
    int line_width() const;
    Qt::PenStyle line_style() const;

    QBrush brush() const;
    QPen placeholder_pen() const;
    QPen path_pen() const;

    QPainterPath placeholder_path() const;
    QPainterPath vessel_path() const;

    //: Cache points of interest in the image (e.g. mouse position, nearest
    //  line pixel to the mouse position, etc.)
    void update_nearest_points_to(const QPointF& pos);

    //: Current state of the mouse interaction
    bool is_adding_arterial_points() const;
    bool is_adding_venous_points() const;

    //: Functions reflecting what you can and can't do with this vessel
    bool is_selected() const;
    bool can_set_size() const;
    bool can_set_shape() const;
    bool can_be_deleted() const;
    bool can_add_path() const;
    bool can_delete_path() const;
    bool can_trim_path() const;
    bool can_add_venous_points() const;
    bool can_add_arterial_points() const;
    bool can_add_apex() const;
    bool can_add_haemorrhage() const;


    // Variables

    //: Appearance options for *all* vessels (hence it is static)
    static ncm_qvesselitem_appearance appearance_;

    bool is_hovered_over_;
    bool is_adding_arterial_points_;
    bool is_adding_venous_points_;

    mutable QRectF boundingRect;

    // cached points of interest (cursor position and points closest to it that
    // lie on a line, edge, or the corresponding vessel path)
    QPointF cursor_pos_;
    QPoint nearest_line_point_;
    QPoint nearest_edge_point_;
    QPoint nearest_vessel_point_;

    //: List of apexitems attached to this vessel
    vcl_vector<QGraphicsApexItem*> ApexItems_;

private:

    // pointer to the class containing the geometry and properties of the vessel
    // This is private so that you are forced to use vessel_data() that checks
    // the validity of the pointer
    ncm_vessel* vessel_data_;
};

// Define the static member variable
ncm_qvesselitem_appearance QGraphicsVesselItemPrivate::appearance_;


QGraphicsVesselItemPrivate::QGraphicsVesselItemPrivate()
: is_hovered_over_(false),
  is_adding_arterial_points_(false),
  is_adding_venous_points_(false),
  ApexItems_(0),
  cursor_pos_(QPointF(0.0,0.0)),
  nearest_line_point_(QPoint(0,0)),
  nearest_edge_point_(QPoint(0,0)),
  nearest_vessel_point_(QPoint(0,0))
{
}

QGraphicsVesselItemPrivate::~QGraphicsVesselItemPrivate()
{
}

//
//: Return handle to ncm_qscene rather than default QGraphicsScene
ncm_qscene* QGraphicsVesselItemPrivate::ncm_scene()
{
  Q_Q(QGraphicsVesselItem);
  return dynamic_cast<ncm_qscene*>(q->scene());
}
const ncm_qscene* QGraphicsVesselItemPrivate::ncm_scene() const
{
  Q_Q(const QGraphicsVesselItem);
  return dynamic_cast<const ncm_qscene*>(q->scene());
}

//
//: Associated vessel
void QGraphicsVesselItemPrivate::set_vessel_data(ncm_vessel* const vessel_data)
{
  // FIXME: What if vessel_data_ is nonzero? Delete old data?

  vessel_data_ = vessel_data;
  
  if (vessel_data_ == NULL)
    return;

  // Create any existing apices (path is constructed when needed)
  const unsigned n_apices = vessel_data_->n_apices();
  for (unsigned i = 0; i < n_apices; ++i)
    create_apexitem_from(vessel_data_->apex(i), /* interactive = */ false);
}

ncm_vessel* QGraphicsVesselItemPrivate::vessel_data()
{
  assert(vessel_data_ != NULL);
  return vessel_data_;
}

const ncm_vessel* QGraphicsVesselItemPrivate::vessel_data() const
{
  assert(vessel_data_ != NULL);
  return vessel_data_;
}

//
//: Is vesselitem selected in list?
bool QGraphicsVesselItemPrivate::is_selected() const
{
  return (vessel_data() == ncm_scene()->selected_vessel());
}

//
//: Create a new apex associated with the vessel, and a corresponding apexitem
//  so that we can see it
QGraphicsApexItem* QGraphicsVesselItemPrivate::create_apexitem_from(ncm_apex* apex,
                                                                    bool interactive /* = true */)
{
  Q_Q(QGraphicsVesselItem);

  // create an ncm_qapexitem attached to the new ncm_apex
  QGraphicsApexItem* new_apexitem = 
      new QGraphicsApexItem(apex, interactive, /* parent = */ q);
  ApexItems_.push_back(new_apexitem);

  // Connect the apex's signals to the scene
  QObject::connect(new_apexitem, SIGNAL(lengthChanged(double)),
                   ncm_scene(), SIGNAL(apexLengthChanged(double)) );

  return new_apexitem;
}

void QGraphicsVesselItemPrivate::create_apexitem_at(const QPointF& first_point)
{
  Q_Q(QGraphicsVesselItem);

  // sanity check
  assert(can_add_apex());

  // create an ncm_apex attached to the ncm_vessel associated with this vessel
  ncm_apex* new_apex = 
      vessel_data()->add_apex_at(first_point.x(), first_point.y());

  QGraphicsApexItem* new_apexitem = create_apexitem_from(new_apex,
                                                         /* interactive = */ true);

  // transfer mouse focus to apexitem
  q->ungrabMouse();
  new_apexitem->grabMouse();
}

//
//: Set colour
QColor QGraphicsVesselItemPrivate::colour() const
{
  // Set colour, depending on edit mode and selected status
  QColor rgba;

  switch (ncm_scene()->edit_mode())
  {
    case ncm_qscene::ModeClassifyImage:
      // Probably best not to show vessels here but we'll define colours anyway
      // (fall through)
    case ncm_qscene::ModeAddVessels:
      if (vessel_data()->properties().is_distal())
        rgba.setRgb(appearance_.colour_normal());
      else
        rgba.setRgb(appearance_.colour_undefined());
      break;

    case ncm_qscene::ModeSetVesselSize:
      if (vessel_data()->properties().is_size_undefined())
        rgba.setRgb(appearance_.colour_undefined());
      else if (vessel_data()->properties().is_size_enlarged())
        rgba.setRgb(appearance_.colour_enlarged());
      else if (vessel_data()->properties().is_size_giant() ||
               vessel_data()->properties().is_size_irregular())
        rgba.setRgb(appearance_.colour_giant());
      else
        rgba.setRgb(appearance_.colour_normal());
      break;

    case ncm_qscene::ModeSetVesselShape:
      if (vessel_data()->properties().is_shape_undefined())
        rgba.setRgb(appearance_.colour_undefined());
      else if (vessel_data()->properties().is_shape_tortuous())
        rgba.setRgb(appearance_.colour_tortuous());
      else if (vessel_data()->properties().is_shape_ramified())
        rgba.setRgb(appearance_.colour_ramified());
      else
        rgba.setRgb(appearance_.colour_normal());
      break;

    case ncm_qscene::ModeLabelApices:
      // fall through
    case ncm_qscene::ModeDrawVesselPath:
      if (is_selected())
        rgba.setRgb(appearance_.colour_selected());
      else
        rgba.setRgb(appearance_.colour_normal());
      break;

    case ncm_qscene::ModeAddHaemorrhages:
      rgba.setRgb(appearance_.colour_undefined());
      break;

		case ncm_qscene::ModeDisplayAll:
      if (vessel_data()->properties().is_shape_auto())
			{
				if (vessel_data()->properties().is_distal())
					rgba.setRgb(appearance_.colour_auto_distal());
				else
					rgba.setRgb(appearance_.colour_auto_nondistal());
			}
      else
        rgba.setRgb(appearance_.colour_undefined());
      break;

    default:
      // All options accounted for - we shouldn't be here.
      assert(false);
  }

  return rgba;
}

//
//: Return opacity (dependent on active/inactive status)
int QGraphicsVesselItemPrivate::opacity() const
{  
  // Set opacity, depending on edit mode and selected status

  switch (ncm_scene()->edit_mode())
  {
    case ncm_qscene::ModeLabelApices:
      // fall through
    case ncm_qscene::ModeDrawVesselPath:
      if (is_selected())
        return appearance_.opacity_selected();
      else
        return appearance_.opacity();
      break;

    default:
      return appearance_.opacity();
  }
}

//
//: Return line width (dependent on vessel properties)
int QGraphicsVesselItemPrivate::line_width() const
{
  const double normal_line_width = appearance_.line_width() * 
                                   ncm_scene()->height();

  if (vessel_data()->properties().is_size_enlarged())
    return normal_line_width * appearance_.scale_enlarged();

  else if (vessel_data()->properties().is_size_giant() || 
           vessel_data()->properties().is_size_irregular())
    return normal_line_width * appearance_.scale_enlarged() 
                             * appearance_.scale_giant();

  else
    return normal_line_width;
}

//
//: Return line width (dependent on vessel properties)
Qt::PenStyle QGraphicsVesselItemPrivate::line_style() const
{
  //if (vessel_data()->properties().is_tortuous())
  //  return Qt::DashLine;
  //else if (vessel_data()->properties().is_ramified())
  //  return Qt::DotLine;
  //else if (vessel_data()->properties().is_something_else())
  //  return Qt::DashDotLine	
  //else
    return Qt::SolidLine;
}

//
//: Set colour and opacity for arterial limb
QBrush QGraphicsVesselItemPrivate::brush() const
{
  QColor rgba(colour());
  rgba.setAlpha(opacity());

  return QBrush(rgba);
}

//
//: Set additional pen properties (i.e. width) for placeholder
QPen QGraphicsVesselItemPrivate::placeholder_pen() const
{
  if (is_selected() && 
      appearance_.outline_when_selected())
    return QPen(brush(), 0);
  else
    return path_pen();
}
//
//: Set additional pen properties (i.e. width) for path
QPen QGraphicsVesselItemPrivate::path_pen() const
{
  // Always draw line at full thickness
  return QPen(brush(), line_width(), line_style());
}

//
//: Update nearest points (i.e. nearest points on a line/edge/vessel)
void QGraphicsVesselItemPrivate::update_nearest_points_to(const QPointF& pos)
{
  // current mouse location
  cursor_pos_ = pos;

  // nearest point on a detected centreline
  nearest_line_point_ = ncm_scene()->line_pixel_nearest_to(pos);

  // nearest point on a detected edge
  nearest_edge_point_ = ncm_scene()->edge_pixel_nearest_to(pos);

  // nearest point on the vessel path (if it exists)
  if (vessel_data()->path_is_defined())
  {
    const vgl_point_2d<double>* nearest_vessel_point = 
        vessel_data()->point_nearest_to(pos.x(),pos.y());
    nearest_vessel_point_ = QPoint(nearest_vessel_point->x(),
                                   nearest_vessel_point->y());
  }
  else
  {
    nearest_vessel_point_ = QPoint();
  }
}

//
//: Path for the placeholder only
QPainterPath QGraphicsVesselItemPrivate::placeholder_path() const
{
  QPainterPath path;
  const double placeholder_radius = appearance_.placeholder_radius() * 
                                    ncm_scene()->imageRect().height();

  const bool is_giant = vessel_data()->properties().is_size_giant();

  QPointF centre(vessel_data()->anchor().x(),
                 vessel_data()->anchor().y());
  if (is_giant)
  {
    // represent vessel by a 'G' centred on the only marked point
    const QRectF rect(centre.x()-placeholder_radius,
                      centre.y()-placeholder_radius,
                      2*placeholder_radius, 2*placeholder_radius);
    path.moveTo(centre);
    path.arcTo(rect,
               /* startAngle = */ 0, 
               /* sweepLength = */ -300);
  }
  else
  {
    // represent vessel by a circle centred on the only marked point
    path.addEllipse(centre, placeholder_radius, placeholder_radius);
  }

  return path;
}

//
//: Path for the vessel centreline
QPainterPath QGraphicsVesselItemPrivate::vessel_path() const
{
  unsigned n_points = vessel_data()->n_points();

  QPainterPath path;
  if (n_points > 0)
  {
    path.moveTo(QPointF(vessel_data()->point(0)->x(),
                        vessel_data()->point(0)->y()));
    for (unsigned i = 1; i < n_points; ++i)
    {
      path.lineTo(QPointF(vessel_data()->point(i)->x(),
                          vessel_data()->point(i)->y()));
    }
  }

  return path;
}

//
//: True if we are adding arterial points
bool QGraphicsVesselItemPrivate::is_adding_arterial_points() const
{
  return is_adding_arterial_points_;
}

//
//: True if we are adding points
bool QGraphicsVesselItemPrivate::is_adding_venous_points() const
{
  return is_adding_venous_points_;
}

//
// Private methods
//

bool QGraphicsVesselItemPrivate::can_set_size() const
{
  const bool correct_mode = 
    ncm_scene()->edit_mode() == ncm_qscene::ModeSetVesselSize;
  
  return (correct_mode &&
          is_hovered_over_
         );
}

bool QGraphicsVesselItemPrivate::can_set_shape() const
{
  const bool correct_mode = 
    ncm_scene()->edit_mode() == ncm_qscene::ModeSetVesselShape;
  
  // Giant vessels are Normal shape by definition.
  const bool is_giant =
    vessel_data()->properties().is_size_giant();

  return (correct_mode &&
          is_hovered_over_ &&
          !is_giant
         );
}

bool QGraphicsVesselItemPrivate::can_be_deleted() const
{
  return (
          (ncm_scene()->edit_mode() == ncm_qscene::ModeAddVessels) &&
          //!(vessel_data()->path_is_defined()) && 
          //(vessel_data()->n_apices() == 0) && 
          //is_selected() &&
          is_hovered_over_
         );
}

bool QGraphicsVesselItemPrivate::can_add_path() const
{
  return (
          (ncm_scene()->edit_mode() == ncm_qscene::ModeDrawVesselPath) &&
          !(vessel_data()->path_is_defined()) && 
          is_selected()
         );

}

bool QGraphicsVesselItemPrivate::can_delete_path() const
{
  return (
          (ncm_scene()->edit_mode() == ncm_qscene::ModeDrawVesselPath) &&
          (vessel_data()->path_is_defined()) && 
          is_selected() && 
          is_hovered_over_
         );
}

bool QGraphicsVesselItemPrivate::can_trim_path() const
{
  return (
          (ncm_scene()->edit_mode() == ncm_qscene::ModeDrawVesselPath) &&
          vessel_data()->path_is_defined() && 
          is_selected() &&
          is_hovered_over_
         );
}

bool QGraphicsVesselItemPrivate::can_add_venous_points() const
{
  // must be in right mode and path must exist already
  if (ncm_scene()->edit_mode() != ncm_qscene::ModeDrawVesselPath ||
      !(vessel_data()->path_is_defined()))
    return false;

  // get current mouse location as a vgl_point_2d
  vgl_point_2d<double> cursor_pos_vgl(cursor_pos_.x(), cursor_pos_.y());

  // if very close to the venous endpoint then add venous points
  // due to the earlier test, venous_endpoint() is guaranteed to be non-NULL
  const vgl_vector_2d<double> pos_to_venous = 
      *(vessel_data()->venous_endpoint()) - cursor_pos_vgl;

  const double dist_to_venous = pos_to_venous.length();

  // TODO: Create a user-defined variable to set this sensitivity
  const int threshold = line_width();

  return (dist_to_venous < threshold);
}

bool QGraphicsVesselItemPrivate::can_add_arterial_points() const
{
  // TODO: consider comparing distance from mouse pointer to both endpoints and
  // allowing addition to whichever is nearer.

  // Must be in right mode and path must exist already
  // Also defines a preference for adding venous over arterial points. This is 
  // an arbitrary decision but ensures that the order of if..elseif..else 
  // statements has less effect.
  if (ncm_scene()->edit_mode() != ncm_qscene::ModeDrawVesselPath ||
      !(vessel_data()->path_is_defined()) ||
      can_add_venous_points())
    return false;

  // get current mouse location as a vgl_point_2d
  vgl_point_2d<double> cursor_pos_vgl(cursor_pos_.x(), cursor_pos_.y());

  // if very close to the arterial endpoint then add arterial points
  // due to the earlier test, arterial_endpoint() is guaranteed to be non-NULL
  const vgl_vector_2d<double> pos_to_arterial = 
      *(vessel_data()->arterial_endpoint()) - cursor_pos_vgl;

  const double dist_to_arterial = pos_to_arterial.length();

  // TODO: Create a user-defined variable to set this sensitivity
  const int threshold = line_width();

  return (dist_to_arterial < threshold);
}

bool QGraphicsVesselItemPrivate::can_add_apex() const
{
  // get current mouse location as a vgl_point_2d
  vgl_point_2d<double> cursor_pos_vgl(cursor_pos_.x(), cursor_pos_.y());

  // Check whether this vessel is the one nearest the current mouse position
  ncm_vessel const* nearest_vessel = 
      ncm_scene()->annotation()->vessel_nearest_to(cursor_pos_vgl);
  const bool this_is_nearest = (vessel_data() == nearest_vessel);

  return (
          (ncm_scene()->edit_mode() == ncm_qscene::ModeLabelApices) &&
          this_is_nearest &&
          is_selected()
         );
}

bool QGraphicsVesselItemPrivate::can_add_haemorrhage() const
{
  return (ncm_scene()->edit_mode() == ncm_qscene::ModeAddHaemorrhages) &&
          false;
}
