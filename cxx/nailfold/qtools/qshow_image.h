#ifndef __QSHOW_IMAGE_H__
#define __QSHOW_IMAGE_H__

#include <vil/vil_image_view.h>

#include <QApplication>
#include <QMainWindow>

inline
//void qshow_image(const vil_image_view<vxl_byte>& image)
void qshow_image()
{
  QApplication a;
  QMainWindow w;
  w.setGeometry(100, 100, 100, 100);
  w.show();
  return a.exec();
}

#endif __QSHOW_IMAGE_H__