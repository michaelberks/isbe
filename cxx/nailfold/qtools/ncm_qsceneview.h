#ifndef ncm_qsceneview_h_
#define ncm_qsceneview_h_

//:
// \file
// \brief View including ability to zoom in and out
//        Adapted from Tim's qvcr_zoom_view
// \author Phil Tresadern

#include <QGraphicsView>
#include <QKeyEvent>
#include <QObject>

#include <vcl_vector.h>

class ncm_qscene;

//: View including ability to zoom in and out
// Use +/- buttons or mouse wheel to zoom
class ncm_qsceneview : public QGraphicsView
{
  Q_OBJECT

public:

  //: Allow both Pan and Markup at same time (toggled using <Ctrl>)
  enum ViewMode { ModePan, ModePanMarkup, ModeMarkup };

  //: Default constructor
  ncm_qsceneview(QWidget* parent = 0);

  //: Construct and set scene
  ncm_qsceneview(QGraphicsScene* scene, QWidget* parent = 0);

  //: Destruct
  ~ncm_qsceneview();

  //: Reverse zoom direction with respect to current value
  void reverse_zoom();

  //: Get current mode
  ViewMode current_mode();

  //: Set to Pan mode
  void usePanMode();

  //: Set to Pan/Markup mode
  void usePanMarkupMode();

  //: Set to Markup mode
  void useMarkupMode();

  //: True if using Pan or PanMarkup mode
  bool isPanning() const;

  double fitWidthScale() const;
  double fitHeightScale() const;
  double fitScale() const;
  double maxScale() const;
  void setScale(double new_scale);
  void setMaxScale(double max_scale);

  ////: Zoom factor in real terms (1.0 = fit in view)
  //void setZoomFactor(double zoomFactor);
  //double zoomFactor() const;

  //: Return pointer to ncm_scene (as opposed to default scene() that returns a
  //  QGraphicsScene
  ncm_qscene* ncm_scene();

signals:
  //: Signal emitted if the view is moved or zoomed
  void changed();

  //: Contrast has changed automatically
  void contrast_changed();
  
  //: Signal emitted if the mode changes
  void viewModeChanged(int);

  //: Signal emitted if zoom level changes

  //: Signals corresponding to progressing through a sequence
  //  e.g. pressing the right arrow to go to the next vessel
  void next();
  void previous();

public slots:
  //: Automatically set contrast parameters
  void setAutoContrast(bool auto_contrast = true);

  //: Set view to display all of scene rectangle
  void fit_both();
  void fit_width();
  void fit_height();

  //: Invoke file dialog and save screenshot
  void save_screenshot();

  //: Save screenshot to named file
  void save_screenshot(const QString& image_path);

  //: Centre on a particular vessel
  void centre_on_vessel(int vessel_index);

  //: Focus on a particular vessel
  void focus_on_vessel(int vessel_index, 
                       double zoom = 40.0, double translate_y = 0.0);

protected:

  //: Update view transform to retain visible area
  void resizeEvent(QResizeEvent* event);

  //: What to do every time the view gets repainted
  void paintEvent(QPaintEvent* event);

  //: React to key presses
  void keyPressEvent(QKeyEvent *event);
  void keyReleaseEvent(QKeyEvent *event);

  //: Use mouse-wheel movement to scale view
  void wheelEvent(QWheelEvent *event);

  //: Scale the view by given factor
  void scaleView(qreal scaleFactor);

  //: Scale the view about a given point by given factor
  void scaleViewAbout(const QPoint& fixed_point, qreal scaleFactor);


private:
  //: Sets various properties
  void set_defaults();

  //: Cursor graphic to use when marking up
  QCursor markup_cursor_;

  //: Determines zoom direction with respect to wheel movement
  int zoom_dir_;

  //: Whether automatic contrast setting is on
  //  If it is, emit changed() at every repaint
  bool auto_contrast_;

  //: Maximum scale of transformation matrix
  //  Anything above this will be truncated to to max_scale_
  double max_scale_;
    
  ViewMode current_mode_;
};


#endif
