#include "vil_suppress_non_max_dir.h"

#include <vcl_cmath.h>

#include <vil/vil_image_view.h>
#include <vil/vil_bilin_interp.h>

void vil_suppress_non_max_dir(const vil_image_view<double>& intensity,
                              const vil_image_view<double>& orientation,
                              vil_image_view<double>& peaks,
                              bool normals /* = false */)
{
  peaks.set_size(intensity.ni(), intensity.nj());

  // If we don't do this, the 1 pixel border around the edge is filled with a
  // random value that can be neither 0 nor 1 (I got 205). This messes up
  // everything when you subsequently want to process the image (e.g. stretch
  // values to the range 0..255)
  peaks.fill(0);

  // cache image properties
  const unsigned ni = intensity.ni();
  const unsigned nj = intensity.nj();
  const int intensity_istep = intensity.istep();
  const int intensity_jstep = intensity.jstep();
  const int orientation_istep = orientation.istep();
  const int orientation_jstep = orientation.jstep();
  const int peaks_istep = peaks.istep();
  const int peaks_jstep = peaks.jstep();

  // skip the first row
  double const* intensity_row_ptr = intensity.top_left_ptr() + intensity_jstep;
  double const* orientation_row_ptr = orientation.top_left_ptr() + orientation_jstep;
  double* peaks_row_ptr = peaks.top_left_ptr() + peaks_jstep;

  for (unsigned j = 1; j < nj-1; ++j)
  {
    // skip the first column
    double const* intensity_ptr = intensity_row_ptr + intensity_istep;
    double const* orientation_ptr = orientation_row_ptr + orientation_istep;
    double* peaks_ptr = peaks_row_ptr + peaks_istep;

    for (unsigned i = 1; i < ni-1; ++i)
    {
      double di, dj;
      if (normals)
      {
        di = vcl_cos(*orientation_ptr); 
        dj = vcl_sin(*orientation_ptr);
      }
      else
      {
        di = -vcl_sin(*orientation_ptr);
        dj = vcl_cos(*orientation_ptr);
      }

      // note that dj is negated in order to use the image coordinate frame
      if ( (*intensity_ptr > vil_bilin_interp(intensity, i+di, j-dj)) &&
           (*intensity_ptr > vil_bilin_interp(intensity, i-di, j+dj)) )
        *peaks_ptr = *intensity_ptr;

      intensity_ptr += intensity_istep;
      orientation_ptr += orientation_istep;
      peaks_ptr += peaks_istep;
    }

    intensity_row_ptr += intensity_jstep;
    orientation_row_ptr += orientation_jstep;
    peaks_row_ptr += peaks_jstep;
  }
}

// Do same but straight to a binary image
void vil_suppress_non_max_dir(const vil_image_view<double>& intensity,
                              const vil_image_view<double>& orientation,
                              vil_image_view<bool>& peaks,
                              bool normals /* = false */)
{
  peaks.set_size(intensity.ni(), intensity.nj());

  // If we don't do this, the 1 pixel border around the edge is filled with a
  // random value that can be neither 0 nor 1 (I got 205). This messes up
  // everything when you subsequently want to process the image (e.g. stretch
  // values to the range 0..255)
  peaks.fill(false);

  // cache image properties
  const unsigned ni = intensity.ni();
  const unsigned nj = intensity.nj();
  const int intensity_istep = intensity.istep();
  const int intensity_jstep = intensity.jstep();
  const int orientation_istep = orientation.istep();
  const int orientation_jstep = orientation.jstep();
  const int peaks_istep = peaks.istep();
  const int peaks_jstep = peaks.jstep();

  // skip the first row
  double const* intensity_row_ptr = intensity.top_left_ptr() + intensity_jstep;
  double const* orientation_row_ptr = orientation.top_left_ptr() + orientation_jstep;
  bool* peaks_row_ptr = peaks.top_left_ptr() + peaks_jstep;

  for (unsigned j = 1; j < nj-1; ++j)
  {
    // skip the first column
    double const* intensity_ptr = intensity_row_ptr + intensity_istep;
    double const* orientation_ptr = orientation_row_ptr + orientation_istep;
    bool* peaks_ptr = peaks_row_ptr + peaks_istep;

    for (unsigned i = 1; i < ni-1; ++i)
    {
      double di, dj;
      if (normals)
      {
        di = vcl_cos(*orientation_ptr); 
        dj = vcl_sin(*orientation_ptr);
      }
      else
      {
        di = -vcl_sin(*orientation_ptr);
        dj = vcl_cos(*orientation_ptr);
      }

      // note that dj is negated in order to use the image coordinate frame
      *peaks_ptr = 
        ((*intensity_ptr > vil_bilin_interp(intensity, i+di, j-dj)) &&
         (*intensity_ptr > vil_bilin_interp(intensity, i-di, j+dj)));

      intensity_ptr += intensity_istep;
      orientation_ptr += orientation_istep;
      peaks_ptr += peaks_istep;
    }

    intensity_row_ptr += intensity_jstep;
    orientation_row_ptr += orientation_jstep;
    peaks_row_ptr += peaks_jstep;
  }
}
