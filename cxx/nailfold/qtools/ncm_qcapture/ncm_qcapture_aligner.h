#ifndef __ncm_qcapture_aligner_h__
#define __ncm_qcapture_aligner_h__

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

#include <nailfold/qtools/ncm_qframe_aligner.h>

#include <QObject>
#include <QPoint>

//  MB: To allow two versions of the frame aligner - one that aligns from files and one that aligns
// from images kept in memory on the fly, I've subclassed ncm_qframe_aligner
class ncm_qcapture_aligner : public ncm_qframe_aligner 
{
  Q_OBJECT // needed if we want to use signals or slots

//  INTERFACE

public: // methods

  //: Constructor.
  ncm_qcapture_aligner(int max_length = 1e4);

	//: Set the pixel resolution of the motors
	void setPixelsPerMm(double pmm);

	//: Set the number of frames that should be aligned.
  void setNumToAlign(unsigned n_to_align);

  //: 
  void reset();

public slots:

  //: Align the next pair of frames to be done.
  void alignNextFrame();

private slots:

private: // variables
	//: Maxium length
	const int MAX_LENGTH;

	double pixels_per_mm_;

	//: Use the motor positions?
	bool use_motor_positions_;
	double motor_x0_;
	double motor_y0_;

  //: List of filenames (with path) to images to be aligned.
  vcl_vector<vcl_string> filenames_;

};
  
#endif //__ncm_qcapture_aligner_h__