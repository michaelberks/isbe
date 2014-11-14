#include "ncm_qframe_aligner.h"
 
#include <vil/vil_load.h>

#include <vil/algo/vil_gauss_reduce.h>
 
//
// Public methods
//
 
ncm_qframe_aligner::ncm_qframe_aligner()
: sequence_is_valid_(false),
  n_aligned_(0),
	n_to_align_(0)
{
}
 
//: Return the number of frames aligned so far.
unsigned ncm_qframe_aligner::numAligned()
{
  return n_aligned_;
}

unsigned ncm_qframe_aligner::numToAlign()
{
  return n_to_align_;
}
 
//: Manually set estimated displacements.
//  e.g. from motor positions
void ncm_qframe_aligner::setDisplacements(
  vcl_vector<QPoint> displacements)
{
  displacements_ = displacements;
}
 
//: Manually set estimated displacement at index i.
//  e.g. from motor positions
void ncm_qframe_aligner::setDisplacement(
  unsigned i,
  QPoint displacement_i)
{
  if (i < displacements_.size())
    displacements_[i] = displacement_i;
}
 
//: Set all displacements to (0,0).
void ncm_qframe_aligner::resetDisplacements()
{
	displacements_.clear();
	displacements_.resize(n_to_align_, QPoint(0.0,0.0));
	n_aligned_ = 0;//MB: filenames_.size()
  //displacements_ = vcl_vector<QPoint>(n_to_align_);//MB:filenames_.size(), QPoint(0,0)
}
  
//: Return vector of displacements.
vcl_vector<QPoint> ncm_qframe_aligner::displacements()
{
  return displacements_;
}
  
//: Return displacement at index i.
QPoint ncm_qframe_aligner::displacement(
  unsigned i)
{
  if (i < displacements_.size())
    return displacements_[i];
  else
    return QPoint();
}
 
//: Return gain at index i.
double ncm_qframe_aligner::gain(
  unsigned i)
{
  if (i < gains_.size())
    return gains_[i];
  else
    return -1.0;
}
 
//: Return offset at index i.
double ncm_qframe_aligner::offset(
  unsigned i)
{
  if (i < offsets_.size())
    return offsets_[i];
  else
    return -1.0;
}
 
//: Return to initial state where destination is first image.
void ncm_qframe_aligner::reset(vil_image_view<vxl_byte> img)
{
  set_destination(img);
  set_size(img.ni(), img.nj());

  // First frame is 'aligned' with the origin.
  n_aligned_ = 1;

  // Reset photometric alignment parameters.
	gains_.clear();
	offsets_.clear();
	
  gains_.resize(n_to_align_, 1.0);//MB: filenames_.size()
  offsets_.resize(n_to_align_, 0.0); //MB: filenames_.size()

}
 
//
// Public slots
//
 
//: Align the current frame with the previous one.
void ncm_qframe_aligner::alignNextFrame(vil_image_view<vxl_byte> img, bool is_last_frame)
{
  if (n_aligned_ < n_to_align_)
  {
    set_source(img);

    // Use an existing estimated displacement if it exists.
    QPoint d0 = displacements_[n_aligned_];

		//If the existing displacement (e.g. from motor pos) is greater than the frame size
		//don't try and align the frame or scale the grey-levels - we'll just have to put
		//it where it is and hope it matches (could warn the user here or make this behaviour optional?)

		int di = d0.x(), dj = d0.y();
		double mse = 0;

		if (vcl_abs(di) < (img.ni()/2) && vcl_abs(dj) < (img.nj()/2))
		{
			align_src_to_dest(di, dj);
			get_displacements(di, dj);

			displacements_[n_aligned_] = QPoint(di, dj);

			mse = mean_squared_error();

			double gain, offset;
			align_greylevels(gain, offset);

			// Update cumulative scale and offset.
			// (Offset uses previous scale so must be updated first.)
			offsets_[n_aligned_] = offsets_[n_aligned_-1] + 
														 gains_[n_aligned_-1] * offset;
			gains_[n_aligned_] = gains_[n_aligned_-1] * gain;
		}

    swap_src_with_dest();

    // Emit a signal to indicate that the frames have been aligned, along with 
    // the corresponding translation.
    emit framesAligned(di, dj, 
                       gains_[n_aligned_], offsets_[n_aligned_], mse);

    ++n_aligned_;
  }

  if (is_last_frame)
	{
		//Discard any excess from the gains, offsets and displacements
		gains_.resize(n_aligned_);
		offsets_.resize(n_aligned_);
		displacements_.resize(n_aligned_);

    emit alignmentFinished();
	}
}
 
//
// Private methods
//
