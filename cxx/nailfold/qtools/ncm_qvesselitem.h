#ifndef __ncm_qvesselitem_h__
#define __ncm_qvesselitem_h__

//:
// \file
// \brief QGraphicsItem-derived class for creating a vessel visualization that
//        will go in a QGraphicsScene (via addItem()).
// \author Phil Tresadern

#include <QGraphicsItem>

#include <nailfold/qtools/ncm_qglobal.h>

class ncm_vessel;
class ncm_qvesselitem_appearance;

class QGraphicsApexItem;

// Most graphics item classes are derived from the QAbstractGraphicsShapeItem
// class, rather than directly from QGraphicsItem. Using the d-pointer pattern,
// however, the QAbstractGraphicsShapeItemPrivate structure is defined in the
// source file (not the header) so we can't derive the 
// QGraphicsVesselItemPrivate structure from it. I had to abandon the 
// QAbstractGraphicsShapeItem altogether for this reason.

// Forward declaration of private data structure
class QGraphicsVesselItemPrivate;

class /*Q_GUI_EXPORT*/ QGraphicsVesselItem : public QGraphicsItem
{
public:
    QGraphicsVesselItem(QGraphicsItem *parent = 0
#ifndef Q_QDOC
                      // ### obsolete argument
                      , QGraphicsScene *scene = 0
#endif
    );

    QGraphicsVesselItem(ncm_vessel* const vessel_data,
                        QGraphicsItem *parent = 0
#ifndef Q_QDOC
                      // ### obsolete argument
                      , QGraphicsScene *scene = 0
#endif
    );

    ~QGraphicsVesselItem();

    QPainterPath path() const;
    QPainterPath normal_path() const;

    virtual QPainterPath shape() const;
    virtual QRectF boundingRect() const;
    virtual void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = 0);

    //: Update the bounding rectangle (used when switching between edit modes
    //  where the appearance of the vesselitem changes)
    void updateBoundingRect();

    //: True if adding points at the moment
    bool is_representing(const ncm_vessel* const vessel) const;
    bool is_selected() const;

    bool is_adding_points() const;

    void remove_apexitem(const QGraphicsApexItem* apexitem);

    //: Appearance options for *all* vessels (hence it is static)
    static ncm_qvesselitem_appearance& appearance();

    // Need to use UserType when defining new classes of QGraphicsItem
    // type numbers are defined in ncm_qglobal.h
    enum { Type = UserType + NcmQt::vesselitem_typenum };
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
    Q_DECLARE_PRIVATE(QGraphicsVesselItem)
    Q_DISABLE_COPY(QGraphicsVesselItem)
};

#endif // __ncm_qvesselitem_h__
