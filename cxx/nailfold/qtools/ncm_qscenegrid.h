#ifndef __ncm_qscenegrid_h__
#define __ncm_qscenegrid_h__

//:
// \file
// \brief QGraphicsItem-derived class for creating a vessel visualization that
//        will go in a QGraphicsScene (via addItem()).
// \author Phil Tresadern

#include <QGraphicsItem>

//#include <nailfold/qtools/ncm_qglobal.h>

#include <vcl_vector.h>

// Most graphics item classes are derived from the QAbstractGraphicsShapeItem
// class, rather than directly from QGraphicsItem. Using the d-pointer pattern,
// however, the QAbstractGraphicsShapeItemPrivate structure is defined in the
// source file (not the header) so we can't derive the 
// QGraphicsSceneGridPrivate structure from it. I had to abandon the 
// QAbstractGraphicsShapeItem altogether for this reason.

// Forward declarations
class ncm_qscene;

class QGraphicsSceneGrid
{
public:
  QGraphicsSceneGrid(ncm_qscene *scene);
  ~QGraphicsSceneGrid();

  void setGridSpacing(double spacing);
  double gridSpacing() const;

  void setVisible(bool visible);
  void show();
  void hide();
  void refresh();

private:
  void createGrid();
  void clearGrid();

  //: Parent scene
  ncm_qscene* scene_;

  //: Spacing (in pixels) between grid lines
  double spacing_;

  //: Whether the items are visible or not
  bool is_visible_;

  //: Vector of line items
  vcl_vector<QGraphicsLineItem*> gridLines_;
};

#endif // __ncm_qscenegrid_h__
