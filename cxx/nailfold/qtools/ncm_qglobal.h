#ifndef __ncm_qglobal_h__
#define __ncm_qglobal_h__

//:
// \file
// \brief  Common global identifiers for use throughout ncm_qt
// \author Phil Tresadern

#include <QRgb>

struct NcmQt {
  //: Standard colours
  static const QRgb transparent = 0x00000000; // = qRgba(0,0,0,0);

  static const QRgb black =   0xff000000; // = qRgba(0,0,0,255);
  static const QRgb white =   0xffffffff; // = qRgba(255,255,255,255);
  static const QRgb grey =    0xff808080; // = qRgba(128,128,128,255);

  static const QRgb red =     0xffff0000; // = qRgba(255,0,0,255);
  static const QRgb green =   0xff00ff00; // = qRgba(0,255,0,255);
  static const QRgb blue =    0xff0000ff; // = qRgba(0,0,255,255);

  static const QRgb cyan =    0xff00ffff; // = qRgba(255,0,0,255);
  static const QRgb magenta = 0xffff00ff; // = qRgba(0,255,0,255);
  static const QRgb yellow =  0xffffff00; // = qRgba(255,255,0,255);

  // Manchester purple :)
  static const QRgb purple =  0xff6d009d; // = qRgba(109,0,157,255);

  //: Identifiers for custom QGraphicsItem classes
  static const int vesselitem_typenum = 1;
  static const int apexitem_typenum = 2;
  static const int haemorrhageitem_typenum = 3;
};

#endif // __ncm_qglobal_h__
