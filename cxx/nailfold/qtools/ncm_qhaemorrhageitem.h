#ifndef __ncm_qhaemorrageitem_h__
#define __ncm_qhaemorrageitem_h__

//:
// \file
// \brief QGraphicsItem-derived class for creating a visualization of a
//        haemorrhage
// \author Phil Tresadern

#include <QKeyEvent>
#include <QObject>
#include <QGraphicsItem>

#include <vcl_vector.h>

#include <nailfold/ncm_vessel.h>

#include <nailfold/qtools/ncm_qglobal.h>
//#include <nailfold/qtools/ncm_qhaemorrhageitem_appearance.h>

class ncm_vessel;
class ncm_haemorrhage;
class ncm_qscene;
class QPainterPath;
//class QGraphicsVesselItem;
//class QScene;
class ncm_qhaemorrhageitem_appearance;

// Most graphics item classes are derived from the QAbstractGraphicsShapeItem
// class, rather than directly from QGraphicsItem. Using the d-pointer pattern,
// however, the QAbstractGraphicsShapeItemPrivate structure is defined in the
// source file (not the header) so we can't derive the 
// QGraphicsHaemorrhageItemPrivate structure from it. I had to abandon the 
// QAbstractGraphicsShapeItem altogether for this reason.

// Forward declaration of private data structure
class QGraphicsHaemorrhageItemPrivate;

class /*Q_GUI_EXPORT*/ QGraphicsHaemorrhageItem : public QGraphicsItem
{
public:
    QGraphicsHaemorrhageItem(ncm_haemorrhage* const haemorrhage_data,
                             QGraphicsItem *parent = 0
#ifndef Q_QDOC
                      // ### obsolete argument
                      , QGraphicsScene *scene = 0
#endif
    );

    ~QGraphicsHaemorrhageItem();

    QPainterPath path() const;

    //: True if adding points at the moment
    bool is_adding_points() const;

    virtual QRectF boundingRect() const;
    virtual QPainterPath shape() const;
    virtual void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = 0);

    //: Appearance options for *all* vessels (hence it is static)
    static ncm_qhaemorrhageitem_appearance& appearance();

    // Need to use UserType when defining new classes of QGraphicsItem
    // type numbers are defined in ncm_qglobal.h
    enum { Type = UserType + NcmQt::haemorrhageitem_typenum };
    int type() const;

protected:
    bool supportsExtension(Extension extension) const;
    void setExtension(Extension extension, const QVariant &variant);
    QVariant extension(const QVariant &variant) const;

    virtual void hoverEnterEvent(QGraphicsSceneHoverEvent *event);
    virtual void hoverMoveEvent(QGraphicsSceneHoverEvent *event);
    virtual void hoverLeaveEvent(QGraphicsSceneHoverEvent *event);

    virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
    virtual void mousePressEvent(QGraphicsSceneMouseEvent *event);
    virtual void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);


private:
    Q_DISABLE_COPY(QGraphicsHaemorrhageItem)
    Q_DECLARE_PRIVATE(QGraphicsHaemorrhageItem)
};

#endif // __ncm_qhaemorrageitem_h__
