#ifndef vil_find_connected_components_h_
#define vil_find_connected_components_h_

#include <vcl_vector.h>

#include <vil/vil_image_view.h>

void vil_find_connected_components(
        const vil_image_view<bool>& binary_image,
        vil_image_view<int>& labels,
        vcl_vector<unsigned>& component_sizes);

#endif // vil_find_connected_components_h_