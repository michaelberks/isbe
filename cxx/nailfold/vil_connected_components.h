#ifndef vil_find_connected_components_h_
#define vil_find_connected_components_h_

#include <vcl_vector.h>

#include <vil/vil_image_view.h>

class vil_connected_components
{
public:
  // Constructors
  vil_connected_components();
  vil_connected_components(const vil_image_view<bool>& binary_image);

  //: Clear the image data
  void clear();
    
  //: Use a different image
  void set_image(const vil_image_view<bool>& binary_image);

  //: Number of components
  unsigned n_components();

  //: Vector of component sizes
  const vcl_vector<unsigned>& component_sizes();

  //: Access image of component labels
  const vil_image_view<unsigned>& label_image();

  //: Access image of component sizes
  const vil_image_view<vxl_byte>& size_image();

  //: Set binning of size image
  enum binning_strategy_type { BinningOff = 0,
                               BinningEqualWidths,
                               BinningEqualCounts };
  void set_binning_strategy(binning_strategy_type binning_strategy);

private:

  //:
  void get_neighbour_labels_of(unsigned i, unsigned j);

  void get_label_range(int& min_label, int& max_label);

  void fill_size_image();

  //: Neighbourhood size (4- or 8-connected)
  enum neighbourhood_type { FourConnected, EightConnected };
  neighbourhood_type neighbourhood_size_;

  //: True if size image is up to date
  bool size_image_is_valid_;

  //: Number of bins (default 255)
  unsigned n_bins_;

  //: Binning strategy
  binning_strategy_type binning_strategy_;


  //  Workspace variables

  //: Component labels
  vil_image_view<unsigned> label_image_;

  //: Component sizes
  vil_image_view<vxl_byte> size_image_;

  unsigned n_components_;
  vcl_vector<unsigned> component_sizes_;

  unsigned n_neighbours_;
  vcl_vector<int> neighbour_labels_;
};

#endif // vil_find_connected_components_h_