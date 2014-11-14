//:
// \file
// \brief QLabel adapted to provide a clickable thumbnail of an image
// \author Phil Tresadern

#include "ncm_qthumbnail.h"

//
//: Constructor
ncm_qthumbnail::ncm_qthumbnail(const QPixmap& pixmap, 
                               const QString& filename,
                               QWidget* parent /* = 0 */)
: filename_(filename.toStdString()),
  pressed_(false)
{
  const int pixmap_width = pixmap.width();

  setFixedSize(pixmap_width, pixmap_width*0.75);
  setPixmap(pixmap);
  setToolTip(filename);
}

//
//: What to do when the mouse is pressed over this thumbnail
void ncm_qthumbnail::mousePressEvent(QMouseEvent *event)
{
  pressed_ = true;
}

//
//: What to do when the mouse is released over this thumbnail
void ncm_qthumbnail::mouseReleaseEvent(QMouseEvent *event)
{
  // Let the receiving widget handle what happens (e.g. load an image)
  if (pressed_)
  {
    pressed_ = false;
    emit clicked(filename_);
  }
}