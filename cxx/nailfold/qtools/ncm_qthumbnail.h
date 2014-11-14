#ifndef ncm_qthumbnail_h_
#define ncm_qthumbnail_h_

//:
// \file
// \brief QLabel adapted to provide a clickable thumbnail of an image
// \author Phil Tresadern

#include <vcl_string.h>

#include <QObject>
#include <QLabel>
#include <QKeyEvent>

class ncm_qthumbnail : public QLabel
{
  Q_OBJECT

public:

  //: Default constructor
  ncm_qthumbnail(const QPixmap& pixmap, 
                 const QString& filename,
                 QWidget* parent = 0);

signals:
  void clicked(vcl_string);

public slots:

protected:
  //: React to mouse activity
  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);

private:

  vcl_string filename_;
  bool pressed_;
};


#endif
