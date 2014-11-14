#ifndef __ncm_qscene_h__
#define __ncm_qscene_h__

//:
// \file
// \brief View including ability to zoom in and out
//        Adapted from Tim's qvcr_zoom_view
// \author Phil Tresadern

#include "ncm_qvessel_properties.h"

#include <QGraphicsView>
#include <QKeyEvent>
#include <QObject>

#include <vcl_vector.h>
#include <vcl_string.h>

#include <vil/vil_image_view.h>

#include <nailfold/vil_connected_components.h>

//class vcl_string;
class ncm_vessel;
class ncm_annotation;
//class vil_connected_components;
class ncm_qimagehandler;

class QGraphicsVesselItem;
class QGraphicsHaemorrhageItem;
class QGraphicsSceneGrid;

class ncm_qscene : public QGraphicsScene
{
  Q_OBJECT // needed if we want to use signals or slots

// INTERFACE

public:
  //: Data type to signify the current mode of editing
  //  (The indices also correspond to the page numbers of the Toolbox)
  //  Consider making this private and adding access functions such as
  //  bool isModeAddVessel()
  enum editMode { ModeUndefined = -1,
                  ModeClassifyImage = 0,
                  ModeAddVessels = 1,
                  ModeSetVesselSize = 2,
                  ModeSetVesselShape = 3,
                  ModeLabelApices = 4,
                  ModeDrawVesselPath = 5, 
                  ModeAddHaemorrhages = 6,
									ModeDisplayAll = 7,
  
                  ModeFirst = ModeClassifyImage,
                  ModeLast = ModeDisplayAll};

  //: Default constructor
  ncm_qscene(QWidget* parent = 0,
             ncm_annotation* markup = 0);

  //: Destruct
  ~ncm_qscene();

  void hideItems();
  void showItems();

  //: Delete all vessel items
  void delete_all_vesselitems();
  void delete_all_haemorrhageitems();

  //: Set the image processor being used
  void set_image_processor(ncm_qimagehandler* image_processor);
  ncm_qimagehandler* image_processor();

  //: Return the bounding rectangle of the image (without any border)
  QRectF imageRect() const;
  
  //: Set the size of a border to surround the image
  void set_border_size(double border_size = 0.0);

  //: Gets nearest non-zero pixel position from supplied coordinates
  QPoint line_pixel_nearest_to(const QPointF& pos) const;
  QPoint edge_pixel_nearest_to(const QPointF& pos) const;
  
  //  Access functions for ncm_qscene properties
  bool is_snapping() const;

  void set_annotation(ncm_annotation* const markup);
  ncm_annotation* annotation();
  ncm_annotation const* annotation() const;

  editMode edit_mode() const;

  //: Selected vessel
  ncm_vessel const* selected_vessel() const;
  QGraphicsVesselItem* selected_vesselitem();

  //: New vessel properties
  ncm_qvessel_properties& new_vessel_properties();
  const ncm_qvessel_properties& new_vessel_properties() const;

  //: VesselItem at given position
  QGraphicsVesselItem* vesselitem_at(const QPointF& pos) const;

  //: Update scene to match markup
  void update_from_annotation();

  //: Helper functions that allow a QGraphicsItem (e.g. a vessel) to emit a
  //  qscene's signal. This looks clumsy but is the best way I can think of.
  //  Signals are protected and therefore can't be called from other
  //  (non-child) classes. Normally, you would create a signal from the other
  //  class and connect the two signals; QGraphicsItem, however, does not 
  //  inherit from QObject and therefore cannot generate signals.
  void emit_vessel_deleted();
  void emit_vessel_changed();

signals:
  void vesselAdded();
  void vesselDeleted();
  void vesselChanged();
  void editModeChanged(int);
  void selectionChanged(int);
  void apexLengthChanged(double);
  void mouseMoved(double, double);

public slots:

  //: Resize boundary of scene to fit raw_pixmapitem_ perfectly
  void fit_to_image();

  //: Set edit mode (label vessel, label apices, label properties, etc)
  void set_edit_mode(int edit_mode_index);

  //: Delete vessels and haemorrhages, and clear background images
  void clear();
	void clear_markup();
	void hide_bitmaps();

  //: Recompute grid lines
  void set_grid_spacing(double grid_spacing);
  void refresh_grid();

  // on/off options
  void set_lines_visible(bool lines_visible);
  void set_edges_visible(bool edges_visible);
  void set_vessels_visible(bool vessels_visible);
  void set_snapping(bool snap_to_lines);
  void set_grid_visible(bool grid_visible);

  // select a specific vessel (necessary for marking path and labelling apices)
  void select_vessel(int vessel_index);

  void update_raw_pixmapitem();
  void update_line_pixmapitem();
  void update_edge_pixmapitem();


// IMPLEMENTATION

protected:

  virtual void mousePressEvent(QGraphicsSceneMouseEvent *event);
  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  virtual void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);

private:

  //: Make changes to scene and annotation
  void create_vesselitem_at(const QPointF& anchor);
  void create_haemorrhageitem_at(const QPointF& anchor);

  //: Edit mode types: should correspond to items in UI toolBox widget
  editMode edit_mode_;

  //: Properties to apply to any new vessel (or selected vessels)
  ncm_qvessel_properties new_vessel_properties_;

  //: Pointer to image processor
  ncm_qimagehandler* image_processor_;

  //: Class that draws a grid over the scene
  QGraphicsSceneGrid* sceneGrid_;

  //: Pixmaps that make up the background of the scene
  QGraphicsPixmapItem* raw_pixmap_item_;
  QGraphicsPixmapItem* line_pixmap_item_;
  QGraphicsPixmapItem* edge_pixmap_item_;

  //: Annotation corresponding to this scene
  ncm_annotation* annotation_;

  //: Currently selected vessel (possibly NULL)
  ncm_vessel* selected_vessel_;

  //: True if detected centrelines should be visible
  bool lines_visible_;

  //: True if detected edges should be visible
  bool edges_visible_;

  //: True if vessels should be visible
  bool vessels_visible_;

  //: True if we want to snap to lines
  bool snap_to_lines_;

  unsigned n_vessels_;

  //: Size of the border to add to the image in the scene
  double border_size_;
};

#endif // __ncm_qscene_h__
