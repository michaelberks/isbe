#ifndef vil_suppress_non_max_dir_h_
#define vil_suppress_non_max_dir_h_

#include <vil/vil_image_view.h>

void vil_suppress_non_max_dir(const vil_image_view<double>& intensity,
                              const vil_image_view<double>& orientation,
                              vil_image_view<double>& peaks,
                              bool normals = false);

void vil_suppress_non_max_dir(const vil_image_view<double>& intensity,
                              const vil_image_view<double>& orientation,
                              vil_image_view<bool>& peaks,
                              bool normals = false);

#endif // vil_suppress_non_max_dir_h_