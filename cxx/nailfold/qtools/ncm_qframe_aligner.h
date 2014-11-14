#ifndef __ncm_qframe_aligner_h__
#define __ncm_qframe_aligner_h__

//:
// \file
// \brief Qt extension to ncm_frame_aligner
//
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_vector.h>
#include <vcl_string.h>

#include <nailfold/ncm_frame_aligner.h>

#include <QObject>
#include <QPoint>

//  To enable signals and slots, the object must be a derived class of QObject.
//  To do this I've used multiple inheritance, though this is known to be
//  not for the fainthearted. We'll see how it goes. Interestingly, the order
//  of the base classes (QObject then ncm_vessel_properties) makes a difference
//  in this instance. Also note the keyword 'public' before both parent classes
class ncm_qframe_aligner : public QObject, 
                           public ncm_frame_aligner 
{
  Q_OBJECT // needed if we want to use signals or slots

//  INTERFACE

public: // methods

  //: Constructor.
  ncm_qframe_aligner();

  //: Set the number of frames that should be aligned.
  virtual void setNumToAlign(
    unsigned n_to_align) = 0;

  //: Return the number of frames aligned so far.
  unsigned numAligned();
	unsigned numToAlign();

  //: Manually set estimated displacements.
  //  e.g. from motor positions
  void setDisplacements(
    vcl_vector<QPoint> displacements);

  //: Manually set estimated displacement at index i.
  //  e.g. from motor positions
  void setDisplacement(
    unsigned i,
    QPoint displacement_i);

  //: Set all displacements to (0,0).
  void resetDisplacements();
  
  //: Return vector of estimated displacements.
  vcl_vector<QPoint> displacements();

  //: Return displacement at index i.
  QPoint displacement(
    unsigned i);

  //: Return gain at index i.
  double gain(
    unsigned i);

  //: Return offset at index i.
  double offset(
    unsigned i);

  //: Return the current image.
  vil_image_view<vxl_byte> current_image();

  //: 
	virtual void reset() = 0;
  

signals:

  void framesAligned(int di, int dj,
                     double scale, double offset, 
                     double mse);
  void alignmentFinished();

public slots:

  //: Align the next pair of frames to be done.
  virtual void alignNextFrame() = 0;

private slots:

protected: 
	//Functions
	void reset(vil_image_view<vxl_byte> img);
	void alignNextFrame(vil_image_view<vxl_byte> img, bool is_last_frame);
	// variables

  //: Current estimate of displacements between images.
  //  First item is always (0,0).
  vcl_vector<QPoint> displacements_;

  //: Cumulative gain that matches current frame to first frame.
  vcl_vector<double> gains_;
  
  //: Cumulative offset that matches current frame to first frame.
  vcl_vector<double> offsets_;

  //: Number of frames to align.
  unsigned n_to_align_;

  //: Index of currently aligned image.
  unsigned n_aligned_;

  //: True if the sequence is valid.
  bool sequence_is_valid_;

private:

};
  
#endif __ncm_qframe_aligner_h__