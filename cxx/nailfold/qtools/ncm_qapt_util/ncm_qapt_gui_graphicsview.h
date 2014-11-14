#ifndef ncm_qapt_gui_graphicsview_h_
#define ncm_qapt_gui_graphicsview_h_

//:
// \file
// \brief View including ability to zoom in and out
//        Adapted from Tim's qvcr_zoom_view
// \author Phil Tresadern

#include <QGraphicsView>
#include <QObject>

//: View including ability to zoom in and out
// Use +/- buttons or mouse wheel to zoom
class ncm_qapt_gui_graphicsview : public QGraphicsView
{
  Q_OBJECT

public:

  //: Default constructor
  ncm_qapt_gui_graphicsview(QWidget* parent = 0);

  //: Destruct
  ~ncm_qapt_gui_graphicsview();

signals:

  void clicked(double x, double y);

public slots:

protected:

  //: React to mouse button press/move/release
  void mousePressEvent(QMouseEvent *ev);

  //: Use mouse-wheel movement to scale view
  void wheelEvent(QWheelEvent *ev);

private:

};


#endif
