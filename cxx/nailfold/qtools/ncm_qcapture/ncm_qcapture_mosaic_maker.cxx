#include "ncm_qcapture_mosaic_maker.h"
#include "ncm_qcapture_data_manager.h"
#include <vil/vil_load.h>

//
// Public methods
//

ncm_qcapture_mosaic_maker::ncm_qcapture_mosaic_maker(
  unsigned input_ni,
  unsigned input_nj)
: ncm_mosaic_maker(input_ni, input_nj),
  displacements_(0),
  positions_(0),
	filenames_(0)
{
}

//
// Public slots
//

//: Create a mosaic from a sequence of images and their offsets.
void ncm_qcapture_mosaic_maker::makeMosaicFrom(
  const vcl_vector<QPoint> displacements)
{
	bool use_data_queue = filenames_.empty();

	if (use_data_queue)
		n_to_stitch_ = ncm_qcapture_data_manager::Instance()->mosaic_queue()->count();
	else
		n_to_stitch_ = filenames_.size();

  assert(n_to_stitch_ == displacements.size());

	vcl_cout << "Num frames to make mosiac from: " << n_to_stitch_ << vcl_endl;

  set_displacements(displacements);
  positions_.resize(0);

	n_stitched_ = 0;

  // Convert QPoints to vector<x> and vector<y>
  vcl_vector<int> di(n_to_stitch_), dj(n_to_stitch_);
  for (unsigned i = 0; i < n_to_stitch_; ++i)
  {
    di[i] = displacements[i].x();
    dj[i] = displacements[i].y();
		//vcl_cout << "Displacement for frame " << i << "x = " << di[i] << ", y = " << dj[i] << vcl_endl;
  }

	bool success = get_mosaic_limits(di, dj);

	
	for (unsigned i = 0; i < n_to_stitch_; ++i)
	{
		vil_image_view<vxl_byte> image;
		if (use_data_queue)
		{
			//Even if we can't mosaic, we need to pop all these frames from the queue
			QSharedPointer< ncm_video_frame > vid_frame = 
				ncm_qcapture_data_manager::Instance()->remove_from_mosaic_queue();
			image = *(vid_frame->frame());
		}
		else //Load frame from file
		{
			image = vil_load(filenames_[i].c_str());
		}

		if (success)
		{
			//Adds image to mosaic and increments n_stitched_
			add_image(image, di[i], dj[i]);
		}
	}

  emit mosaicFinished(success);

  cache_positions();
}
 
//:
void ncm_qcapture_mosaic_maker::add_image(
  const vil_image_view<vxl_byte> img, 
  unsigned int di, unsigned int dj)
{
  ncm_mosaic_maker::add_image(img, di, dj);
  ++n_stitched_;

  emit mosaicUpdated(n_stitched_);
}

//: Return a vector of indices
vcl_vector<unsigned> ncm_qcapture_mosaic_maker::images_at(
  QPointF pos,
  double width /* = 128 */,
  double height /* = 128 */
  )
{
  vcl_vector<unsigned> indices(0);

  if (!positions_.empty())
  {
    const double half_width = 0.5 * width;
    const double half_height = 0.5 * height;

    for (unsigned i = 0; i < positions_.size(); ++i)
    {
      const int left = positions_[i].x();
      const int right = positions_[i].x() + input_ni_;
      const int top = positions_[i].y();
      const int bottom = positions_[i].y() + input_nj_;

      const bool is_in_view =  left < (pos.x()-half_width) &&
                                      (pos.x()+half_width) < right &&
                                top < (pos.y()-half_height) && 
                                      (pos.y()+half_height) < bottom;

      if (is_in_view)
        indices.push_back(i);
    }
  }

  return indices;
}
 
//: 
QPoint ncm_qcapture_mosaic_maker::position(unsigned i)
{
  if (i < positions_.size())
    return positions_[i];
  else
    return QPoint();
}
 
//: Set the number of frames that should be stitched.
void ncm_qcapture_mosaic_maker::setNumToStitch(
  unsigned n_to_stitch)
{
  //if (n_to_stitch <= filenames_.size())
  //  n_to_align_ = n_to_align;
  //else
  //  n_to_align_ = filenames_.size();
}

void ncm_qcapture_mosaic_maker::setFilenames(
    vcl_vector<vcl_string> filenames)
{
	filenames_ = filenames;
}
 
//: Return the number of frames stitched so far.
unsigned ncm_qcapture_mosaic_maker::numStitched()
{
  return n_stitched_;
}
 
//
// Private methods
//

//
// Private slots
//
 
//: 
void ncm_qcapture_mosaic_maker::set_displacements(
  const vcl_vector<QPoint> displacements)
{
  displacements_ = displacements;
}
 
//:
void ncm_qcapture_mosaic_maker::set_displacements(
  const vcl_vector<int> di,
  const vcl_vector<int> dj)
{
  assert(di.size() == dj.size());

  displacements_.resize(di.size());
  for (unsigned i = 0; i < di.size(); ++i)
  {
    displacements_[i] = QPoint(di[i], dj[i]);
  }
}
 
//:
void ncm_qcapture_mosaic_maker::cache_positions()
{
  positions_.resize(displacements_.size());
  
  if (positions_.size() > 0)
    positions_[0] = QPoint(origin_i_, origin_j_);

  for (unsigned i = 1; i < displacements_.size(); ++i)
  {
    positions_[i] = positions_[i-1] + displacements_[i];
  }
}