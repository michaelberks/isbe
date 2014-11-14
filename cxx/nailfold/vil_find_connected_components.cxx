#include "vil_find_connected_components.h"

#include <vcl_cmath.h>

#include <vil/vil_image_view.h>
#include <vil/vil_bilin_interp.h>
#include <vil/vil_convert.h>

//
// Get the component labels of neighbouring pixels
void get_neighbours(const vil_image_view<int>& labels,
                    int i, int j,
                    vcl_vector<int>& neighbours,
                    unsigned& n_neighbours)
{
  n_neighbours = 0;

  // clear entries in neighbours vector
  neighbours.assign(4,-1);

  unsigned neighbourhood_size = 8;

  // if not on left edge than add pixel to left
  if ((i > 0) && (labels(i-1,j) > 0))
  {
    neighbours[n_neighbours] = labels(i-1,j);
    n_neighbours++;
  }

  // if not on top edge than add pixel above
  if ((j > 0) && (labels(i,j-1) > 0))
  {
    neighbours[n_neighbours] = labels(i,j-1);
    n_neighbours++;
  }

  // if using 8-neighbourhood, add pixel above-left
  if ((neighbourhood_size == 8) && (j > 0))
  {
    // if not on left edge than add pixel above-left
    if ((i > 0) && (labels(i-1,j-1) > 0))
    {
      neighbours[n_neighbours] = labels(i-1,j-1);
      n_neighbours++;
    }

    // if not on right edge than add pixel above-right
    if ((i < labels.ni()-1) && (labels(i+1,j-1) > 0))
    {
      neighbours[n_neighbours] = labels(i+1,j-1);
      n_neighbours++;
    }
  }
}

//
// Find root nodes of neighbouring components
void get_roots(vcl_vector<int>& roots,
               unsigned n_neighbours,
               const vcl_vector<int>& parent_of)
{
  for (unsigned i = 0; i < n_neighbours; i++)
  {
    // replace component index with index of its root
    while (parent_of[ roots[i]-1 ] > 0)
      roots[i] = parent_of[ roots[i]-1 ];
  }
}

//
//: Find minimum and maximum label values for neighbouring pixels
void get_label_range(const vcl_vector<int> neighbours,
                     unsigned n_neighbours,
                     int& min_label,
                     int& max_label)
{
  // get minimum and maximum label values of neighbours
  min_label = -1;
  max_label = -1;
  for (unsigned n = 0; n < n_neighbours; n++)
  {
    // ignore neighbours that are part of background
    if (neighbours[n] != 0)
    {
      if (min_label == -1)
      {
        // if this is the first nonzero neighbour
        // set both min and max to this value
        min_label = neighbours[n];
        max_label = neighbours[n];
      }
      else
      {
        // otherwise update min and max values accordingly
        if (neighbours[n] > max_label)
          max_label = neighbours[n];
        if (neighbours[n] < min_label)
          min_label = neighbours[n];
      }
    }
  }
}

//
//: Find connected components of binary image
void vil_find_connected_components(
        const vil_image_view<bool>& binary_image,
        vil_image_view<int>& labels,
        vcl_vector<unsigned>& component_sizes)
{
  // take a copy of the binary image
  vil_convert_cast(binary_image, labels);

  vcl_vector<int> parent_of;

  unsigned n_labels = 0;
  unsigned n_neighbours = 0;

  // storage for neighbouring components and their roots in the tree
  vcl_vector<int> neighbours(4);
  vcl_vector<int> roots(4);

  for (unsigned j = 0; j < labels.nj(); ++j)
  {
    for (unsigned i = 0; i < labels.ni(); ++i)
    {
      if (labels(i,j) > 0)
      {
        get_neighbours(labels, i, j, 
                       neighbours, n_neighbours);

        int min_label = -1;
        int max_label = -1;
        get_label_range(neighbours, n_neighbours, 
                        min_label, max_label);

        if (min_label == -1)
        {
          // no non-background neighbours found
          // create a new component
          n_labels++;
          labels(i,j) = n_labels;
          parent_of.push_back(0);
        }
        else if (min_label == max_label)
        {
          // connected to one known component
          // assign to this component
          labels(i,j) = min_label;
        }
        else
        {
          // pixel connects to two different components
          // assign it to the component with lowest label
          labels(i,j) = min_label;

          roots = neighbours;
          get_roots(roots, n_neighbours, parent_of);

          // get lowest root index
          int min_root = -1;
          for (unsigned n = 0; n < n_neighbours; n++)
          {
            if ((n == 0) || (roots[n] < min_root))
              min_root = roots[n];
          }

          // make this root the parent of all the others
          for (unsigned n = 0; n < n_neighbours; n++)
          {
            if (roots[n] != min_root)
              parent_of[ roots[n]-1 ] = min_root;
          }
        }
      }
    }
  }

  // create a vector that assigns labels computed in the first pass
  // with the true connected components
  vcl_vector<int> true_labels(parent_of.size());
  int n_components = 0;
  for (unsigned i = 0; i < parent_of.size(); i++)
  {
    if (parent_of[i] == 0)
    {
      // root node
      n_components++;
      true_labels[i] = n_components;
    }
    else
    {
      // trace back up the tree to find this node's root

      // check that neither this node nor its parent is a root node
      while ((parent_of[i] > 0) && (parent_of[ parent_of[i]-1 ] > 0))
        parent_of[i] = parent_of[ parent_of[i]-1 ];

      true_labels[i] = true_labels[ parent_of[i]-1 ];
    }
  }

  // update label image and record component sizes
  component_sizes.assign(n_components,0);
  for (unsigned j = 0; j < labels.nj(); ++j)
    for (unsigned i = 0; i < labels.ni(); ++i)
      if (labels(i,j) > 0)
      {
        labels(i,j) = true_labels[ labels(i,j)-1 ];
        component_sizes[labels(i,j)-1]++;
      }
}
