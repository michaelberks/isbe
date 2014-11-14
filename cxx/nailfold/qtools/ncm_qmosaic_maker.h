#ifndef __ncm_qmosaic_maker_h__
#define __ncm_qmosaic_maker_h__


#include <nailfold/ncm_mosaic_maker.h>

#include <QObject>
#include <QPoint>

class ncm_qmosaic_maker : public QObject, public ncm_mosaic_maker
{
  Q_OBJECT // needed if we want to use signals or slots

public:

  //:
  ncm_qmosaic_maker(
    unsigned input_ni, 
    unsigned input_nj);

  QPoint position(unsigned i);

  //: Set the number of frames that should be stitched.
  void setNumToStitch(
    unsigned n_to_stitch);

  //: Return the number of frames stitched so far.
  unsigned numStitched();

signals:

  void mosaicUpdated();
  void mosaicFinished();

public slots:

  //:
  void makeMosaicFrom(
    const vcl_vector<vcl_string> filenames,
    const vcl_vector<int> di,
    const vcl_vector<int> dj);

  //:
  void makeMosaicFrom(
    const vcl_vector<vcl_string> filenames,
    const vcl_vector<QPoint> displacements);

  //:
  virtual
  void add_image(
    const vil_image_view<vxl_byte> img,
    unsigned di, unsigned dj);

  //: Return a list of the images that contain a box of size width x height,
  //  centred at x, y.
  vcl_vector<unsigned> images_at(
    QPointF pos,
    double width = 128,
    double height = 128);

private slots:

private: // methods

  //:
  void set_displacements(
    const vcl_vector<QPoint>);

  //:
  void set_displacements(
    const vcl_vector<int> di,
    const vcl_vector<int> dj);

  //:
  void cache_positions();

private: // variables

  vcl_vector<QPoint> displacements_;
  vcl_vector<QPoint> positions_;

  //: Number of frames to align.
  unsigned n_to_stitch_;

  //: Index of currently aligned image.
  unsigned n_stitched_;

};

#endif // __ncm_qmosaic_maker_h__