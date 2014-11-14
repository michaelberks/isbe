#include "vil_connected_components.h"

//#include <vcl_cmath.h>

#include <vil/vil_bilin_interp.h>
#include <vil/vil_convert.h>

#include <mbl/mbl_index_sort.h>

//
// PUBLIC MEMBER FUNCTIONS
//

//
//: Constructors
vil_connected_components::vil_connected_components()
: neighbourhood_size_(EightConnected),
  n_components_(0),
  binning_strategy_(BinningOff),
  n_bins_(255),
  size_image_is_valid_(false)
{
  label_image_.set_size(0,0);
  size_image_.set_size(0,0);

  neighbour_labels_.resize(4);
}

vil_connected_components::vil_connected_components(
    const vil_image_view<bool>& binary_image)
: neighbourhood_size_(EightConnected),
  n_components_(0),
  binning_strategy_(BinningOff),
  n_bins_(255),
  size_image_is_valid_(false)
{
  neighbour_labels_.resize(4);

  set_image(binary_image);
}

//
//: Clear any existing data
void vil_connected_components::clear()
{
  label_image_.clear();
  size_image_.clear();
}

//
//: Find connected components of binary image
void vil_connected_components::set_image(const vil_image_view<bool>& binary_image)
{
  // clear size image
  size_image_.set_size(0,0);

  // take a copy of the binary image
  vil_convert_cast(binary_image, label_image_);

  // parent node of each node
  // set parent of background component to dummy value (-1)
  vcl_vector<int> parent_of_node(1);
  parent_of_node[0] = -1;

  // root node of each of the (up to four) neighbouring pixels
  vcl_vector<int> root_of_neighbour(4);

  unsigned n_labels = 0;
  unsigned n_background = 0;

  const unsigned ni = label_image_.ni();
  const unsigned nj = label_image_.nj();

  for (unsigned j = 0; j < nj; ++j)
  {
    for (unsigned i = 0; i < ni; ++i)
    {
      // ignore background pixels
      if (label_image_(i,j) == 0)
        continue;

      // get labels of neighbouring foreground pixels
      get_neighbour_labels_of(i, j);

      if (n_neighbours_ == 0)
      {
        // neighbours are all background: create a new component
        n_labels++;
        label_image_(i,j) = n_labels;
        parent_of_node.push_back(0);
      }
      else
      {
        // get minimum and maximum label values of neighbours
        int min_label = n_labels+1;
        int max_label = -1;
        for (unsigned n = 0; n < n_neighbours_; n++)
        {
          if (neighbour_labels_[n] > max_label)
            max_label = neighbour_labels_[n];
          if (neighbour_labels_[n] < min_label)
            min_label = neighbour_labels_[n];
        }

        if (min_label == max_label)
        {
          // connected to one known component: assign to this component
          label_image_(i,j) = min_label;
        }
        else
        {
          // pixel connects to two different components: assign it to the 
          // component with lowest label
          label_image_(i,j) = min_label;

          // find root component for all neighbours
          root_of_neighbour = neighbour_labels_;
          int min_root = n_labels;
          for (unsigned n = 0; n < n_neighbours_; n++)
          {
            // replace label with that of its root
            while (parent_of_node[ root_of_neighbour[n] ] > 0)
              root_of_neighbour[n] = parent_of_node[ root_of_neighbour[n] ];

            // keep track of lowest root index
            if (root_of_neighbour[n] < min_root)
              min_root = root_of_neighbour[n];
          }

          // make the lowest root the parent of all neighbouring components
          for (unsigned n = 0; n < n_neighbours_; n++)
          {
            if (root_of_neighbour[n] != min_root)
              parent_of_node[ root_of_neighbour[n] ] = min_root;
          }
        } // if (min_label == max_label)
      } // if (n_neighbours_ == 0)
    } // for i
  } // for j

  // create a vector that assigns labels computed in the first pass
  // with the true connected components
  vcl_vector<int> root_of_node(parent_of_node.size());
  vcl_vector<int> true_label_of_node(parent_of_node.size());
  n_components_ = 0;

  // component 0 is always background so start from node = 1
  for (unsigned n = 1; n < parent_of_node.size(); n++)
  {
    if (parent_of_node[n] == 0)
    {
      // root node: assign to a new component
      n_components_++;
      true_label_of_node[n] = n_components_;
    }
    else
    {
      // branch/leaf node: trace root node of this branch/leaf and assign this
      // node to that component

      // trace up tree to find root of this node
      root_of_node[n] = n;
      while (parent_of_node[ root_of_node[n] ] > 0)
        root_of_node[n] = parent_of_node[ root_of_node[n] ];

      true_label_of_node[n] = true_label_of_node[ root_of_node[n] ];
    }
  }

  // update label image and record component sizes
  component_sizes_.assign(n_components_+1,0);
  component_sizes_[0] = ni*nj;
  for (unsigned j = 0; j < nj; ++j)
    for (unsigned i = 0; i < ni; ++i)
      if (label_image_(i,j) > 0)
      {
        label_image_(i,j) = true_label_of_node[ label_image_(i,j) ];
        component_sizes_[ label_image_(i,j) ]++;
        component_sizes_[0]--;
      }

  // schedule size image for a recompute (if needed)
  size_image_.set_size(label_image_.ni(),label_image_.nj());
  size_image_is_valid_ = false;
}

//
//: Number of components
unsigned vil_connected_components::n_components()
{
  return n_components_;
}

//
//: Vector of component sizes
const vcl_vector<unsigned>& vil_connected_components::component_sizes()
{
  return component_sizes_;
}

//
//: Access image of component labels
const vil_image_view<unsigned>& vil_connected_components::label_image()
{
  return label_image_;
}

//
//: Access image of component sizes
const vil_image_view<vxl_byte>& vil_connected_components::size_image()
{
  // compute size image at point of first use
  if (!size_image_is_valid_)
    fill_size_image();

  return size_image_;
}

//
//: Set whether size image is binned
void vil_connected_components::set_binning_strategy(
    binning_strategy_type binning_strategy)
{
  // if we're changing the value then schedule size_image_ for a recompute
  if (binning_strategy != binning_strategy_)
  {
    size_image_is_valid_ = false;
    binning_strategy_ = binning_strategy;
  }
}

//
// PRIVATE MEMBER FUNCTIONS
//

//
//: Get the component labels of neighbouring foreground pixels
void vil_connected_components::get_neighbour_labels_of(unsigned i, unsigned j)
{
  // clear entries in neighbours vector
  n_neighbours_ = 0;
  neighbour_labels_.assign(4,-1);

  // if not on left edge than add pixel to left
  if ((i > 0) && (label_image_(i-1,j) > 0))
  {
    neighbour_labels_[n_neighbours_] = label_image_(i-1,j);
    n_neighbours_++;
  }

  // if not on top edge than add pixel above
  if ((j > 0) && (label_image_(i,j-1) > 0))
  {
    neighbour_labels_[n_neighbours_] = label_image_(i,j-1);
    n_neighbours_++;
  }

  // if using 8-neighbourhood, add pixel above-left
  if ((neighbourhood_size_ == EightConnected) && (j > 0))
  {
    // if not on left edge than add pixel above-left
    if ((i > 0) && (label_image_(i-1,j-1) > 0))
    {
      neighbour_labels_[n_neighbours_] = label_image_(i-1,j-1);
      n_neighbours_++;
    }

    // if not on right edge than add pixel above-right
    if ((i < label_image_.ni()-1) && (label_image_(i+1,j-1) > 0))
    {
      neighbour_labels_[n_neighbours_] = label_image_(i+1,j-1);
      n_neighbours_++;
    }
  }
}

//
//: Fill size_image with sizes of components
void vil_connected_components::fill_size_image()
{
  size_image_.set_size(label_image_.ni(), label_image_.nj());
  size_image_.fill(0);

  // do nothing more (return) if no components have been detected
  if (n_components_ == 0)
    return;

  const unsigned ni = label_image_.ni();
  const unsigned nj = label_image_.nj();

  if (binning_strategy_ == BinningOff)
  {
    // Replace label with size of component
    for (unsigned j = 0; j < nj; j++)
      for (unsigned i = 0; i < ni; i++)
        if (label_image_(i,j) > 0)
          size_image_(i,j) = static_cast<vxl_byte>(component_sizes_[label_image_(i,j)]);
  }
  else
  {
    // take a copy of line size vector and sort it, storing destinations
    vcl_vector<unsigned> sorted_sizes = component_sizes_;
    vcl_vector<int> index_vector;

    // skip first element (corresponding to background)
    mbl_index_sort(&component_sizes_[1], n_components_, index_vector);

    // create a vector of bin indices (1..n_bins)
    vcl_vector<unsigned> bin_of_component(n_components_+1);
    switch (binning_strategy_)
    {
      case BinningEqualWidths:
        {
          // equally sized bin widths (varying number of lines in each bin)

          // size of biggest component
          unsigned max_size = component_sizes_[index_vector.back()+1];

          bin_of_component[0] = 0;
          for (unsigned i = 0; i < n_components_; ++i)
          {
            bin_of_component[ index_vector[i]+1 ] = 
                static_cast<unsigned>(1 + ((n_bins_-1) * component_sizes_[index_vector[i]+1]) / max_size);
          }
        }
        break;

      case BinningEqualCounts:
        // equal number of lines in each bin (varying bin widths)

        bin_of_component[0] = 0;
        for (unsigned i = 1; i <= n_components_; ++i)
          bin_of_component[ index_vector[i-1]+1 ] = 
              static_cast<unsigned>(1 + ((n_bins_-1) * i) / n_components_);
        break;

      default:
        // should never reach here
        vcl_cerr << "vil_connected_components::fill_size_image error:"
                 << "Unknown binning strategy" << vcl_endl;
    }

    // Replace label with bin index of component
    for (unsigned j = 0; j < nj; j++)
      for (unsigned i = 0; i < ni; i++)
        if (label_image_(i,j) > 0)
          size_image_(i,j) = static_cast<vxl_byte>(bin_of_component[label_image_(i,j)]);
  }

  size_image_is_valid_ = true;
}