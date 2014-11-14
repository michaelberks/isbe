#include "ncm_qstitch_graphicsview.h"

#include <vcl_iostream.h>

#include <QMouseEvent>

//
//  Public methods
//

//: Default constructor
ncm_qstitch_graphicsview::ncm_qstitch_graphicsview(
  QWidget* parent /* = NULL */)
: QGraphicsView(parent),
  rectItem_(NULL)
{
  setMouseTracking(true);
}

//
//  Protected events
//

void ncm_qstitch_graphicsview::mouseMoveEvent(
  QMouseEvent* event)
{
  QPointF scenePos = mapToScene(event->pos());
  emit mouseMoved(scenePos);
}

void ncm_qstitch_graphicsview::mousePressEvent(
  QMouseEvent* event)
{
  QPointF scenePos = mapToScene(event->pos());
  emit mousePressed(scenePos);
}

//
//  Private methods
//
