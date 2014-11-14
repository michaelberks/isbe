//:
// \file
// \brief View including ability to zoom in and out
//        Based on Tim's qvcr_zoom_view
// \author Phil Tresadern

#include "ncm_qscene.h"

#include "ncm_qvesselitem.h"
#include "ncm_qhaemorrhageitem.h"
#include "ncm_qglobal.h"
#include "ncm_qimagehandler.h"

#include <vcl_cmath.h>
#include <vcl_algorithm.h>
#include <vcl_iostream.h>

#include <mbl/mbl_index_sort.h>

#include <nailfold/ncm_vessel.h>
#include <nailfold/ncm_annotation.h>

#include <nailfold/qtools/ncm_qscenegrid.h>

#include <QFileDialog>
#include <QMessageBox>
#include <QtGui>

//
//  Public methods
//

ncm_qscene::ncm_qscene(QWidget * parent, 
                       ncm_annotation* markup /* = 0*/)
: QGraphicsScene(parent), // Base class construction
  snap_to_lines_(false),
  annotation_(markup),
  edit_mode_(ModeUndefined),
  lines_visible_(false),
  edges_visible_(false),
  selected_vessel_(NULL),
  border_size_(0.0),
  image_processor_(NULL)
{
  // add pixmaps now
  raw_pixmap_item_ = addPixmap(QPixmap());
  line_pixmap_item_ = addPixmap(QPixmap());
  edge_pixmap_item_ = addPixmap(QPixmap());

  // set visibility
  line_pixmap_item_->setVisible(lines_visible_);
  edge_pixmap_item_->setVisible(edges_visible_);

  setSceneRect(QRectF(0,0,100,100));

  sceneGrid_ = new QGraphicsSceneGrid(this);
  sceneGrid_->setGridSpacing(500.0);
}

//
//: Destructor
ncm_qscene::~ncm_qscene()
{
  // manually delete all graphicsitems before the scene is destroyed
  delete_all_vesselitems();
  delete_all_haemorrhageitems();

  delete sceneGrid_;
}

//
//: Read access functions
bool ncm_qscene::is_snapping() const
{
  return lines_visible_ && snap_to_lines_;
}

ncm_annotation* ncm_qscene::annotation()
{
  return annotation_;
}
ncm_annotation const* ncm_qscene::annotation() const
{
  return annotation_;
}
ncm_qscene::editMode ncm_qscene::edit_mode() const
{
  return edit_mode_;
}

//
//: Write access functions
void ncm_qscene::set_image_processor(ncm_qimagehandler* image_processor)
{
  image_processor_ = image_processor;
  QObject::connect(image_processor_, SIGNAL(rawImageChanged()),
                   this, SLOT(update_raw_pixmapitem()));
  QObject::connect(image_processor_, SIGNAL(lineImageChanged()),
                   this, SLOT(update_line_pixmapitem()));
  QObject::connect(image_processor_, SIGNAL(edgeImageChanged()),
                   this, SLOT(update_edge_pixmapitem()));
}
ncm_qimagehandler* ncm_qscene::image_processor()
{
  return image_processor_;
}

void ncm_qscene::set_annotation(ncm_annotation* const markup)
{
  annotation_ = markup;
}

//: Return the bounding rectangle of the image (without a border)
QRectF ncm_qscene::imageRect() const
{
  if (image_processor_ != NULL)
  {
    const int width = image_processor_->width();
    const int height = image_processor_->height();
    return QRectF(0, 0, width, height);
  }
  else
    return QRectF(0, 0, 0, 0);
}

//: Set the size of a border to surround the image
void ncm_qscene::set_border_size(double border_size /* = 0.0 */)
{
  if (border_size >= 0.0)
    border_size_ = border_size;
}

void ncm_qscene::fit_to_image()
{
  const int width = image_processor_->width() + 2*border_size_;
  const int height = image_processor_->height() + 2*border_size_;
  setSceneRect(QRectF(-border_size_, -border_size_, width, height));
}


//: Gets nearest non-zero pixel position from supplied coordinates
QPoint ncm_qscene::line_pixel_nearest_to(const QPointF& pos) const
{
  return image_processor_->line_pixel_nearest_to(pos);
}
QPoint ncm_qscene::edge_pixel_nearest_to(const QPointF& pos) const
{
  return image_processor_->edge_pixel_nearest_to(pos);
}

void ncm_qscene::clear()
{
  // Remove all vessels and haemorrhages
  clear_markup();

  // Hide bitmaps in background
  hide_bitmaps();
}

void ncm_qscene::clear_markup()
{
  // Remove all vessels and haemorrhages
  delete_all_vesselitems();
  delete_all_haemorrhageitems();

	if (annotation() != NULL)
		annotation_->clear();
}

void ncm_qscene::hide_bitmaps()
{  
	// Hide bitmaps in background
  raw_pixmap_item_->hide();
  line_pixmap_item_->hide();
  edge_pixmap_item_->hide();
	
}

void ncm_qscene::set_grid_spacing(double grid_spacing)
{
  if (sceneGrid_ != NULL)
    sceneGrid_->setGridSpacing(grid_spacing);
}
void ncm_qscene::refresh_grid()
{
  if (sceneGrid_ != NULL)
    sceneGrid_->refresh();
}

void ncm_qscene::update_from_annotation()
{
  // Create new vesselitems from the associated markup data
  delete_all_vesselitems();
  const unsigned n_vessels = annotation_->n_vessels();
  for (unsigned i = 0; i < n_vessels; ++i)
  {
    // Tell the vesselitem what its scene is during (rather than after)
    // construction.
    QGraphicsVesselItem* dummy = 
        new QGraphicsVesselItem(annotation_->vessel(i), 
                                /* parent = */ NULL, 
                                /* scene = */ this);
  }

  // Create new haemorrhageitems from the associated markup data
  delete_all_haemorrhageitems();
  const unsigned n_haemorrhages = annotation_->n_haemorrhages();
  for (unsigned i = 0; i < n_haemorrhages; ++i)
    addItem(new QGraphicsHaemorrhageItem(annotation_->haemorrhage(i)));
}

//
//: Return (const) pointer to (const) selected vessel
ncm_vessel const* ncm_qscene::selected_vessel() const
{
  return selected_vessel_;
}

//
//: Find vesselitem associated with selected vessel
QGraphicsVesselItem* ncm_qscene::selected_vesselitem()
{
  if (selected_vessel() == NULL)
    return NULL;

  for (int i = 0; i < items().size(); ++i)
  {
    QGraphicsVesselItem* vesselitem = 
        qgraphicsitem_cast<QGraphicsVesselItem*>(items()[i]);

    const bool item_is_vesselitem = (vesselitem != NULL);
    if (item_is_vesselitem && vesselitem->is_selected())
      return vesselitem;
  }

  // We shouldn't be here - every vessel should have an associated vesselitem
  // so throw an error
  assert(false);

  return NULL;
}

//
//: New vessel properties
ncm_qvessel_properties& ncm_qscene::new_vessel_properties()
{
  return new_vessel_properties_;
}
const ncm_qvessel_properties& ncm_qscene::new_vessel_properties() const
{
  return new_vessel_properties_;
}

//
//: Hide all graphics items in the scene
void ncm_qscene::hideItems()
{
  raw_pixmap_item_->hide();
  line_pixmap_item_->hide();
  edge_pixmap_item_->hide();
}

//
//: Show all graphics items in the scene
void ncm_qscene::showItems()
{
  raw_pixmap_item_->show();
  line_pixmap_item_->setVisible(lines_visible_);
  edge_pixmap_item_->setVisible(edges_visible_);
}

//
//: Delete all vessel items
void ncm_qscene::delete_all_vesselitems()
{
  // delete vessels in reverse order, since items() is updated with every call
  // to removeItem() such that vessels change position in the vector and do not
  // get deleted.
  for (int i = items().length()-1; i >= 0; --i)
  {
    QGraphicsVesselItem* vessel_item = 
        qgraphicsitem_cast<QGraphicsVesselItem*>(items()[i]);
    const bool item_is_vesselitem = (vessel_item != NULL);

    // if item is a vesselitem then detach it from the scene and delete it
    // (including any associated entry in the markup)
    if (item_is_vesselitem)
      delete vessel_item;
  }
}

//
//: Delete all vessel items
void ncm_qscene::delete_all_haemorrhageitems()
{
  // delete haemorrhages in reverse order, since items() is updated with every 
  // call to removeItem() such that vessels change position in the vector and 
  // do not get deleted.
  for (int i = items().length()-1; i >= 0; --i)
  {
    QGraphicsHaemorrhageItem* haemorrageitem = 
        qgraphicsitem_cast<QGraphicsHaemorrhageItem*>(items()[i]);
    const bool item_is_haemorrageitem = (haemorrageitem != NULL);

    // if item is a haemorrageitem then detach it from the scene and delete it
    // (including any associated entry in the markup)
    if (item_is_haemorrageitem)
      delete haemorrageitem;
  }
}

//
//: VesselItem at given position
QGraphicsVesselItem* ncm_qscene::vesselitem_at(const QPointF& pos) const
{
  return qgraphicsitem_cast<QGraphicsVesselItem*>(itemAt(pos));
}

//
//: Helper functions that allow a QGraphicsItem to tell the scene to emit a 
//  signal
void ncm_qscene::emit_vessel_deleted()
{
  emit vesselDeleted();
}
void ncm_qscene::emit_vessel_changed()
{
  emit vesselChanged();
}

//
// Public slots
//

//: Set the edit mode by its index
void ncm_qscene::set_edit_mode(int edit_mode_index)
{
  // Store current mode for comparison
  const int old_mode = edit_mode_;

  if ((ModeFirst <= edit_mode_index) && (edit_mode_index <= ModeLast))
    edit_mode_ = static_cast<editMode>(edit_mode_index);
  else
  {
    vcl_cerr << "ncm_qscene::set_edit_mode: " << vcl_endl;
    vcl_cerr << "  Unknown edit mode index: " << edit_mode_index << vcl_endl;
  }

  if (edit_mode_ != old_mode)
  {
    emit editModeChanged(edit_mode_index);

    // Because items look different in different modes, their bounding boxes
    // need to be updated
    for (int i = 0; i < items().size(); ++i)
    {
      QGraphicsVesselItem* vesselitem =
          qgraphicsitem_cast<QGraphicsVesselItem*>(items()[i]);

      if (vesselitem != NULL)
        vesselitem->updateBoundingRect();
    }
  }
}

void ncm_qscene::set_lines_visible(bool lines_visible)
{
  lines_visible_ = lines_visible;

  // update graphics item
  line_pixmap_item_->setVisible(lines_visible_);
}

void ncm_qscene::set_edges_visible(bool edges_visible)
{
  edges_visible_ = edges_visible;

  // update graphics item
  edge_pixmap_item_->setVisible(edges_visible_);
}

void ncm_qscene::set_vessels_visible(bool vessels_visible)
{
  vessels_visible_ = vessels_visible;

  // update graphics items
  for (int i = 0; i < items().size(); ++i)
  {
    QGraphicsVesselItem* vessel = 
          qgraphicsitem_cast<QGraphicsVesselItem*>(items()[i]);

    if (vessel != 0)
      vessel->setVisible(vessels_visible_);
  }

  //update();
}


void ncm_qscene::set_snapping(bool snap_to_lines)
{
  snap_to_lines_ = snap_to_lines;
}

void ncm_qscene::set_grid_visible(bool grid_visible)
{
  if (sceneGrid_ != NULL)
    sceneGrid_->setVisible(grid_visible);
}

void ncm_qscene::select_vessel(int vessel_index)
{
  // Store to check if changed
  ncm_vessel* const old_selection = selected_vessel_;

  if (annotation() != NULL)
    selected_vessel_ = annotation()->vessel(vessel_index);
  else
    selected_vessel_ = NULL;

  // Emit signal only if the selection has changed
  if (selected_vessel_ != old_selection)
    emit selectionChanged(vessel_index);
}

void ncm_qscene::update_raw_pixmapitem()
{
  if (image_processor_ != NULL)
  {
    raw_pixmap_item_->setPixmap(image_processor_->raw_pixmap());
    raw_pixmap_item_->show();
  }
}
void ncm_qscene::update_line_pixmapitem()
{
  if (image_processor_ != NULL)
  {
    line_pixmap_item_->setPixmap(image_processor_->line_pixmap());
    line_pixmap_item_->setVisible(lines_visible_);
  }
}
void ncm_qscene::update_edge_pixmapitem()
{
  if (image_processor_ != NULL)
  {
    edge_pixmap_item_->setPixmap(image_processor_->edge_pixmap());
    edge_pixmap_item_->setVisible(edges_visible_);
  }
}

//
//  Protected events
//

void ncm_qscene::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  // Look for underlying items that will accept this event. If something 
  // accepted the event then we leave it to them. Used to:
  //   Set vessel properties (e.g. enlarged, bushy)
  //   Delete vessel
  //   Append/delete/trim vessel path
  //   Delete/modify apex
  //   Drag and drop haemorrhage?
  QGraphicsScene::mousePressEvent(event);
  if (event->isAccepted())
    return;

  // Otherwise, if a vessel has been selected from the list then consider
  // passing the event to it manually. Used to:
  //   Add path points
  //   Add apices
  if (selected_vesselitem() != NULL)
  {
    // Send subsequent mouse events to selected vesselitem
    selected_vesselitem()->grabMouse(); 
    sendEvent(selected_vesselitem(), event);

    if (event->isAccepted())
      return;
  }

  // Otherwise, we've clicked on empty space with nothing selected. Used to:
  //   Add vessel
  //   Add haemorrhage
  switch (edit_mode_)
  {
    case ModeClassifyImage:
      // do nothing
      break;

    case ModeAddVessels:
      if (event->button() == Qt::LeftButton)
        create_vesselitem_at(event->scenePos());
      break;

    case ModeSetVesselSize:
      // fall through
    case ModeSetVesselShape:
      // fall through
    case ModeDrawVesselPath:
      // fall through
    case ModeLabelApices:
      // All of these cases should have been handled by this point
      break;

    case ModeAddHaemorrhages:
      if (event->button() == Qt::LeftButton)
        create_haemorrhageitem_at(event->scenePos());
      break;

		case ModeDisplayAll:
			// Do nothing
			break;

    default:
      // why are we even here? Throw an error
      assert(false);
  }
}

void ncm_qscene::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  emit mouseMoved(event->scenePos().x(), event->scenePos().y());

  // Pass on mouse event to items underneath the cursor
  // If an item accepts the event then we're done
  QGraphicsScene::mouseMoveEvent(event);
  if (event->isAccepted())
    return;
}

void ncm_qscene::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
  // pass on mouse event to items underneath the cursor
  QGraphicsScene::mouseReleaseEvent(event);
}

//
//  Private methods
//

//
//: Create a new vesselitem at given coordinates
void ncm_qscene::create_vesselitem_at(const QPointF& anchor)
{
  // Check markup exists
  assert(annotation_ != NULL);

  //if ( is_snapping() )
  //  anchor = line_pixel_nearest_to(anchor);

  // create a vessel in the markup with specified properties
  ncm_vessel* new_vessel = 
      annotation_->create_vessel_at(anchor.x(), anchor.y(),
                                    &new_vessel_properties_);

  // create a corresponding vesselitem and add to scene
  QGraphicsVesselItem* new_vesselitem = new QGraphicsVesselItem(new_vessel);
  addItem(new_vesselitem);

  // Let other objects know that a vessel has been added
  emit vesselAdded();
}

void ncm_qscene::create_haemorrhageitem_at(const QPointF& anchor)
{
  // Check markup exists
  assert(annotation_ != NULL);

  //if ( is_snapping() )
  //  anchor = line_pixel_nearest_to(anchor);

  // create a vessel in the markup
  ncm_haemorrhage* new_haemorrhage = 
      annotation_->create_haemorrhage_at(anchor.x(), anchor.y());

  // create a corresponding vesselitem and add to scene
  QGraphicsHaemorrhageItem* new_haemorrhageitem = 
      new QGraphicsHaemorrhageItem(new_haemorrhage);
  addItem(new_haemorrhageitem);

  // Let other objects know that a vessel has been added
  //emit haemorrhage_added();
}