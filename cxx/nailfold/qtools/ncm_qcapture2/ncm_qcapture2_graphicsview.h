#ifndef ncm_qcapture2_graphicsview_h_
#define ncm_qcapture2_graphicsview_h_

//:
// \file
// \brief View including ability to zoom in and out
//        Adapted from Tim's qvcr_zoom_view
// \author Phil Tresadern

#include <QGraphicsView>
#include <QObject>

//: View including ability to zoom in and out
// Use +/- buttons or mouse wheel to zoom
class ncm_qcapture2_graphicsview : public QGraphicsView
{
  Q_OBJECT

public:

  //: Default constructor
  ncm_qcapture2_graphicsview(QWidget* parent = 0);

  //: Destruct
  ~ncm_qcapture2_graphicsview();

signals:

  void clicked(double x, double y);
  void zoomed(int delta);

public slots:

protected:

  //: React to mouse button press/move/release
  void mousePressEvent(QMouseEvent *ev);

  //: Use mouse-wheel movement to scale view
  void wheelEvent(QWheelEvent *ev);

private:

};


#endif
