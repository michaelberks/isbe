//:
// \file
// \brief View including ability to zoom in and out
//        Based on Tim's qvcr_zoom_view
// \author Phil Tresadern

#include "ncm_qapt_gui_graphicsview.h"

#include <vcl_iostream.h>

#include <QMouseEvent>

//
//: Constructor
ncm_qapt_gui_graphicsview::ncm_qapt_gui_graphicsview(QWidget * parent)
  : QGraphicsView(parent)
{
}

//
//: Destructor
ncm_qapt_gui_graphicsview::~ncm_qapt_gui_graphicsview()
{
}

//
//  Event handlers
//

//
//: Respond to mouse presses
void ncm_qapt_gui_graphicsview::mousePressEvent(QMouseEvent *ev)
{
  // Emit a clicked event with (x,y) coordinates.
  QPointF scenePos(mapToScene(ev->pos()));
  emit clicked(scenePos.x(), scenePos.y());

  // Call parent event handler.
  QGraphicsView::mousePressEvent(ev);
}

//
//: Use mouse-wheel movement to scale view
void ncm_qapt_gui_graphicsview::wheelEvent(QWheelEvent *ev)
{
  // TODO
}

//
// Private functions
//

