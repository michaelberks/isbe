//:
// \file
// \brief A grid overlay for the NCM image, giving an indication of scale in
//        intuitive units (e.g. mm)
// \author Phil Tresadern

#include "ncm_qscenegrid.h"
#include "ncm_qglobal.h"

#include <nailfold/qtools/ncm_qscene.h>

//#include <QtGui>
#include <QGraphicsScene>

QGraphicsSceneGrid::QGraphicsSceneGrid(ncm_qscene *scene)
: spacing_(-1.0),
  scene_(scene),
  is_visible_(false)
{
}

QGraphicsSceneGrid::~QGraphicsSceneGrid()
{
  clearGrid();
}

//: Set the spacing (in pixels) between grid lines.
//  If a negative spacing is given, clear the grid altogether.
void QGraphicsSceneGrid::setGridSpacing(double spacing)
{
  // Do nothing if there is no change.
  if (spacing == spacing_)
    return;

  spacing_ = spacing;

  // Get rid of the old grid lines
  clearGrid();

  // Create a new grid if needed
  if (spacing_ > 0.0)
    createGrid();
}

//: Get the spacing (in pixels) between grid lines.
double QGraphicsSceneGrid::gridSpacing() const
{ 
  return spacing_;
}

void QGraphicsSceneGrid::setVisible(bool visible)
{
  is_visible_ = visible;

  // Set every line to visible
  for (unsigned i = 0; i < gridLines_.size(); i++)
    gridLines_[i]->setVisible(is_visible_);
}

void QGraphicsSceneGrid::show()
{
  setVisible(true);
}

void QGraphicsSceneGrid::hide()
{
  setVisible(false);
}

//: Clear old grid and create a new one.
//  Used when scene changes size, for example.
void QGraphicsSceneGrid::refresh()
{
  clearGrid();
  createGrid();
}

//
//  Private methods
//

void QGraphicsSceneGrid::createGrid()
{
  const double sceneWidth = scene_->imageRect().width();
  const double sceneHeight = scene_->imageRect().height();

  const unsigned nLinesX = sceneWidth / spacing_;
  const unsigned nLinesY = sceneHeight / spacing_;

  QGraphicsLineItem* gridLine;

  QColor gridColor(NcmQt::green);
  gridColor.setAlpha(128);
  QPen gridPen(gridColor);

  // Create vertical grid lines.
  double x = spacing_;
  for (unsigned i = 0; i < nLinesX; i++)
  {
    gridLine = scene_->addLine(QLineF(x, 0.0, x, sceneHeight), gridPen);
    gridLine->setVisible(is_visible_);
    gridLines_.push_back(gridLine);
    x += spacing_;
  }

  // Create horizontal grid lines.
  double y = spacing_;
  for (unsigned i = 0; i < nLinesY; i++)
  {
    gridLine = scene_->addLine(QLineF(0.0, y, sceneWidth, y), gridPen);
    gridLine->setVisible(is_visible_);
    gridLines_.push_back(gridLine);
    y += spacing_;
  }
}

void QGraphicsSceneGrid::clearGrid()
{
  // For every line, detach from the scene and delete the item
  for (unsigned i = 0; i < gridLines_.size(); i++)
  {
    scene_->removeItem(gridLines_[i]);
    delete gridLines_[i];
  }

  gridLines_.resize(0);
}