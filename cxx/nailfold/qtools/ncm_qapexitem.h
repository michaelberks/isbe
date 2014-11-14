#ifndef __ncm_qapexitem_h__
#define __ncm_qapexitem_h__

//:
// \file
// \brief QGraphicsItem-derived class for creating a visualization of the two
//        points either side of the apex (during annotation).
// \author Phil Tresadern

#include <QGraphicsItem>

#include <nailfold/ncm_vessel.h>

#include <nailfold/qtools/ncm_qglobal.h>

class ncm_qapexitem_appearance;
class QGraphicsVesselItem;

// Most graphics item classes are derived from the QAbstractGraphicsShapeItem
// class, rather than directly from QGraphicsItem. Using the d-pointer pattern,
// however, the QAbstractGraphicsShapeItemPrivate structure is defined in the
// source file (not the header) so we can't derive the 
// QGraphicsApexItemPrivate structure from it. I had to abandon the 
// QAbstractGraphicsShapeItem altogether for this reason.

// Forward declaration of private data structure
class QGraphicsApexItemPrivate;

class /*Q_GUI_EXPORT*/ QGraphicsApexItem : public QObject, public QGraphicsItem
{
  Q_OBJECT
  Q_INTERFACES(QGraphicsItem)

public:
    QGraphicsApexItem(ncm_apex* const apex_data,
                      bool interactive = false,
                      QGraphicsItem *parent = 0
#ifndef Q_QDOC
                      // ### obsolete argument
                      , QGraphicsScene *scene = 0
#endif
    );

    ~QGraphicsApexItem();

    QPainterPath path() const;

    virtual QRectF boundingRect() const;
    virtual QPainterPath shape() const;
    virtual void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = 0);

    //: Appearance options for *all* apices (hence it is static)
    static ncm_qapexitem_appearance& appearance();

    // Need to use UserType when defining new classes of QGraphicsItem
    // type numbers are defined in ncm_qglobal.h
    enum { Type = UserType + NcmQt::apexitem_typenum };
    int type() const;

signals:
    void lengthChanged(double);

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
    Q_DISABLE_COPY(QGraphicsApexItem)
    //Q_DECLARE_PRIVATE(QGraphicsApexItem)

    // Because this class inherits from both QObject and QGraphicsItem, both of
    // which have a member called d_ptr, we need to be explicit about which
    // d_ptr we want. I've therefore done a copy-paste on the Q_DECLARE_PRIVATE
    // macro from QObject and made the necessary local changes.
    inline QGraphicsApexItemPrivate* d_func() { 
      return reinterpret_cast<QGraphicsApexItemPrivate *>(qGetPtrHelper(QGraphicsItem::d_ptr)); 
    }
    inline const QGraphicsApexItemPrivate* d_func() const { 
      return reinterpret_cast<const QGraphicsApexItemPrivate *>(qGetPtrHelper(QGraphicsItem::d_ptr)); 
    }
    friend class QGraphicsApexItemPrivate;
};

#endif // __ncm_qapexitem_h__
