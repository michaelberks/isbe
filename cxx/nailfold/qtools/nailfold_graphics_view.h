#ifndef nailfold_graphics_view_h_
#define nailfold_graphics_view_h_

//:
// \file
// \brief View including ability to zoom in and out
//        Adapted from Tim's qvcr_zoom_view
// \author Phil Tresadern

#include <QObject>
#include <QGraphicsView>
#include <QKeyEvent>


#include <vcl_vector.h>

//: View including ability to zoom in and out
// Use +/- buttons or mouse wheel to zoom
class nailfold_graphics_view : public QGraphicsView
{
  Q_OBJECT

public:
  //: Default constructor
  nailfold_graphics_view(QWidget* parent = 0);

  //: Construct and set scene
  nailfold_graphics_view(QGraphicsScene* scene, QWidget* parent = 0);

  //: Destruct
  ~nailfold_graphics_view();

  //: Set view to display all of scene rectangle
  void show_all_scene();

  //: Reverse zoom direction with respect to current value
  void reverse_zoom() { zoom_dir_ = -zoom_dir_; }

  //: Compiles a list of all nonzero pixels in background image of current scene
  void get_vessel_pixels(QGraphicsPixmapItem* pixmap_item);

  //: Set all vessels to red
  void clear_nearby_lines();

  //: Whether data has been modified
  bool isModified() { return is_modified_; }

  //: Reset modified status (e.g. after saving)
  void resetModified() { is_modified_ = false; }

public slots:
  //: Invoke file dialog and save screenshot
  void save_screenshot();

  //: Save screenshot to named file
  void save_screenshot(const QString& image_path);

  //: Set snapping on or off
  void set_snapping(bool use_snapping) { use_snapping_ = use_snapping; }

protected:
  //: Update view transform to retain visible area
  void resizeEvent (QResizeEvent* event);

  //: React to key presses
  void keyPressEvent(QKeyEvent *event);

  //: React to mouse button press
  void mousePressEvent(QMouseEvent *event);

  //: React to mouse movement
  void mouseMoveEvent(QMouseEvent *event);

  //: React to mouse button release
  void mouseReleaseEvent(QMouseEvent *event);

  //: Use mouse-wheel movement to scale view
  void wheelEvent(QWheelEvent *event);

  //: Scale the view by given factor
  void scaleView(qreal scaleFactor);


private:
  //: Sets various properties
  void set_defaults();

  //: Gets nearest non-zero pixel position from supplied coordinates
  QPoint pixel_nearest_to(QPointF pos);

  //: Find nearest vessel item to pos
  QGraphicsPathItem* vessel_nearest_to(QPoint pos);

  //: Add a point to a vessel
  void addVesselPoint(QPoint pos);



  //: Determines zoom direction with respect to wheel movement
  int zoom_dir_;

  //: Whether we're snapping to lines
  bool use_snapping_;

  //: Markup has been modified since last save
  bool is_modified_;

  //: Pen style for drawing vessels
  QPen vessel_pen_;

  //: Pen style for drawing vessels
  QPen selected_vessel_pen_;

  //: Workspace variable for drawing lines
  QPolygonF polyline_;

  //: Pointer to line currently being drawn
  QGraphicsPathItem* current_polyline_;

  //: Pointers to nearby lines
  vcl_vector<QGraphicsPathItem*> nearby_lines_;

  //: Nonzero pixels
  vcl_vector<QPoint> nonzero_pixels_;

  //: Pointers to objects that need deleting later
  //  (I'm having a problem deleting paths as I go)
  vcl_vector<QGraphicsItem*> trash_;
};


#endif
