//:
// \file
// \brief View including ability to zoom in and out
//        Based on Tim's qvcr_zoom_view
// \author Phil Tresadern

#include "ncm_qcapture_graphicsview.h"

#include <vcl_iostream.h>

#include <QMouseEvent>

//
//: Constructor
ncm_qcapture_graphicsview::ncm_qcapture_graphicsview(QWidget * parent)
  : QGraphicsView(parent)
{
}

//
//: Destructor
ncm_qcapture_graphicsview::~ncm_qcapture_graphicsview()
{
}

//
//  Event handlers
//

////
////: Respond to mouse presses
//void ncm_qcapture_graphicsview::mousePressEvent(QMouseEvent *ev)
//{
//  // Emit a clicked event with (x,y) coordinates.
//  QPointF scenePos(mapToScene(ev->pos()));
//  emit clicked(scenePos.x(), scenePos.y());
//
//  // Call parent event handler.
//  QGraphicsView::mousePressEvent(ev);
//
//  ev->accept();
//}

//
//: Use mouse-wheel movement to scale view
void ncm_qcapture_graphicsview::wheelEvent(QWheelEvent *ev)
{
  emit zoomed(ev->delta());

  ev->accept();
}

//
// Private functions
//

