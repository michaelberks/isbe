#ifndef __NCM_QSTITCH_GRAPHICSVIEW_H__
#define __NCM_QSTITCH_GRAPHICSVIEW_H__

#include <QGraphicsView>

class ncm_qstitch_graphicsview : public QGraphicsView
{
  Q_OBJECT // needed if we want to use signals or slots

public: // methods

  //: Default constructor
  ncm_qstitch_graphicsview(
    QWidget* parent = NULL);

signals:

  void mouseMoved(QPointF);
  void mousePressed(QPointF);

public slots:

protected: // events

  void mouseMoveEvent(
    QMouseEvent* event);

  void mousePressEvent(
    QMouseEvent* event);

protected: // methods

protected: // variables

private slots:

private: // methods

private: // variables

  QGraphicsItem* rectItem_;

};

#endif __NCM_QSTITCH_GRAPHICSVIEW_H__
