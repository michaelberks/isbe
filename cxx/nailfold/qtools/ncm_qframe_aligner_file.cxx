#include "ncm_qframe_aligner_file.h"
 
#include <vil/vil_load.h>

#include <vil/algo/vil_gauss_reduce.h>
 
//
// Public methods
//
 
ncm_qframe_aligner_file::ncm_qframe_aligner_file()
{
}
 
void ncm_qframe_aligner_file::setFilenames(
  vcl_vector<vcl_string> filenames)
{
  sequence_is_valid_ = (filenames.size() > 1);

  if (sequence_is_valid_)
  {
    filenames_ = filenames;
		n_to_align_ = filenames_.size();
    reset();
		resetDisplacements();
  }
  else
  {
    filenames_.clear();
		n_to_align_ = 0;
  }
}

vcl_vector<vcl_string> ncm_qframe_aligner_file::filenames() const
{
  return filenames_;
}
 
//: Set the number of frames that should be aligned.
void ncm_qframe_aligner_file::setNumToAlign(
  unsigned n_to_align)
{
  if (n_to_align <= filenames_.size())
    n_to_align_ = n_to_align;
  else
    n_to_align_ = filenames_.size();
}
 
//: Return to initial state where destination is first image.
void ncm_qframe_aligner_file::reset()
{
	if (n_to_align_ < 2)
		return;

  vil_image_view<vxl_byte> img = vil_load(filenames_[0].c_str());
	ncm_qframe_aligner::reset(img);

}
 
//
// Public slots
//
 
//: Align the current frame with the previous one.
void ncm_qframe_aligner_file::alignNextFrame()
{
  if (n_aligned_ < n_to_align_)
  {
    vil_image_view<vxl_byte> img = vil_load(filenames_[n_aligned_].c_str());
		ncm_qframe_aligner::alignNextFrame(img, n_aligned_+1 == n_to_align_);
	}
}
 
//
// Private methods
//
