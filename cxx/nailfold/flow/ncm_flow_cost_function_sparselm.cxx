#include "ncm_flow_cost_function_sparselm.h"

#include <vil/vil_math.h>
#include <vil/algo/vil_sobel_3x3.h>
#include <vil/algo/vil_convolve_1d.h>

#include <nailfold/flow/ncm_flow_null_robustifier.h>
#include <nailfold/flow/ncm_flow_field.h>

#include <splm.h>
 
//
//  Public methods
//
 
//: Constructor.
ncm_flow_cost_function_sparselm::ncm_flow_cost_function_sparselm(
  vcl_vector< vil_image_view<float> > const* image_stack,
  vcl_vector< vil_image_view<float> > const* warped_image_stack /* = NULL */,
  ncm_flow_robustifier const* robustifier /* = NULL */)
: ncm_flow_cost_function(image_stack, warped_image_stack, robustifier)
{
  // All of the construction is handled by the base class.
}
 
//: Destructor.
ncm_flow_cost_function_sparselm::~ncm_flow_cost_function_sparselm()
{
}
 
//: Compute the vector of residuals.
//  Because the optimizer will not accept a member function for f(), I've had
//  to provide a pointer to the function as additional data (*adata) and
//  dereference it within this function.

// static
void ncm_flow_cost_function_sparselm::f(
  double* parameters, 
  double* x_hat, 
  int nvars, 
  int nobs, 
  void* adata)
{  
  ncm_flow_cost_function_sparselm* func = 
    static_cast<ncm_flow_cost_function_sparselm*>(adata);

  ncm_flow_field flow_field_to_test(func->ni_, func->nj_);

  flow_field_to_test.set_from_vector(parameters);

  unsigned n_added = 0;
  long f_offset = 0;

  n_added = func->brightness_constancy(flow_field_to_test, x_hat, f_offset);
  f_offset += n_added;
  
  n_added = func->spatial_smoothness(flow_field_to_test, x_hat, f_offset);
  f_offset += n_added;

  //f = robustifier_->f(f_in);
}

//: Specify the Jacobian nonzeros
//  Because the optimizer will not accept a member function for f(), I've had
//  to provide a pointer to the function as additional data (*adata) and
//  dereference it within this function.

// static
void ncm_flow_cost_function_sparselm::crs_jacobian_sparsity_ST(
  double* parameters, 
  struct splm_crsm* jacobian_map, 
  int nvars, 
  int nobs, 
  void* adata)
{
  ncm_flow_cost_function_sparselm* func = 
    static_cast<ncm_flow_cost_function_sparselm*>(adata);

  // Allocate and fill-in a ST Jacobian
  struct splm_stm jac_st;
  splm_stm_alloc(&jac_st, nobs, nvars, jacobian_map->nnz);

  unsigned row_index = 0;

  // Fill the sparsity map with values for the brightness constancy.
  unsigned u_index = 0;
  unsigned v_index = 1;
  for (unsigned j = 0; j < func->nj_; ++j)
  {
    for (unsigned i = 0; i < func->ni_; ++i)
    {
      for (unsigned k = 0; k < func->n_images_; ++k)
      {
        splm_stm_nonzero(&jac_st, row_index, u_index);
        splm_stm_nonzero(&jac_st, row_index, v_index);
        ++row_index;
      }

      u_index += 2;
      v_index += 2;
    }
  }

  // Fill the sparsity map with values for the smoothness constraints.
  ncm_flow_field ff(func->ni_, func->nj_);
  for (unsigned j = 0; j < func->nj_; ++j)
  {
    for (unsigned i = 0; i < func->ni_; ++i)
    {
      if (i > 0)
      {
        splm_stm_nonzero(&jac_st, row_index, ff.index_of_u(i,j));
        splm_stm_nonzero(&jac_st, row_index, ff.index_of_u(i-1,j));
        ++row_index;

        splm_stm_nonzero(&jac_st, row_index, ff.index_of_v(i,j));
        splm_stm_nonzero(&jac_st, row_index, ff.index_of_v(i-1,j));
        ++row_index;
      }

      if (j > 0)
      {
        splm_stm_nonzero(&jac_st, row_index, ff.index_of_u(i,j));
        splm_stm_nonzero(&jac_st, row_index, ff.index_of_u(i,j-1));
        ++row_index;

        splm_stm_nonzero(&jac_st, row_index, ff.index_of_v(i,j));
        splm_stm_nonzero(&jac_st, row_index, ff.index_of_v(i,j-1));
        ++row_index;
      }
    }
  }

  // Convert to CCS and free memory
  splm_stm2crsm(&jac_st, jacobian_map);
  splm_stm_free(&jac_st);
}
 
// static
void ncm_flow_cost_function_sparselm::crs_jacobian(
  double* parameters, 
  struct splm_crsm* jacobian_map, 
  int nvars, 
  int nobs, 
  void* adata)
{
  ncm_flow_cost_function_sparselm* func = 
    static_cast<ncm_flow_cost_function_sparselm*>(adata);

  // Fill the jacobian with values for the brightness constancy.
  unsigned u_index = 0;
  unsigned v_index = 1;
  unsigned row_index = 0;
  unsigned nz = 0;

  const unsigned ni = func->ni_;
  const unsigned nj = func->nj_;
  const unsigned n_images = func->n_images_;

  for (unsigned j = 0; j < nj; ++j)
  {
    for (unsigned i = 0; i < ni; ++i)
    {
      for (unsigned k = 0; k < n_images; ++k)
      {
        jacobian_map->rowptr[row_index] = nz;

        jacobian_map->colidx[nz] = u_index;
        jacobian_map->val[nz] = func->Ix_[k](i,j);
        ++nz;

        jacobian_map->colidx[nz] = v_index;
        jacobian_map->val[nz] = func->Iy_[k](i,j);
        ++nz;

        ++row_index;
      }

      u_index += 2;
      v_index += 2;
    }
  }

  // Fill the jacobian with values for the smoothness constraints.
  ncm_flow_field ff(ni, nj);
  const double weight = (1.0 - func->alpha_);
  const double weight_hv = weight * (n_images / 6);

  for (unsigned j = 0; j < nj; ++j)
  {
    for (unsigned i = 0; i < ni; ++i)
    {
      if (i < ni-1)
      {
        jacobian_map->rowptr[row_index] = nz;
        jacobian_map->colidx[nz] = ff.index_of_u(i,j);
        jacobian_map->val[nz] = weight_hv;
        ++nz;
        jacobian_map->colidx[nz] = ff.index_of_u(i+1,j);
        jacobian_map->val[nz] = -weight_hv;
        ++nz;
        ++row_index;

        jacobian_map->rowptr[row_index] = nz;
        jacobian_map->colidx[nz] = ff.index_of_v(i,j);
        jacobian_map->val[nz] = weight_hv;
        ++nz;
        jacobian_map->colidx[nz] = ff.index_of_v(i+1,j);
        jacobian_map->val[nz] = -weight_hv;
        ++nz;
        ++row_index;
      }

      if (j < nj-1)
      {
        jacobian_map->rowptr[row_index] = nz;
        jacobian_map->colidx[nz] = ff.index_of_u(i,j);
        jacobian_map->val[nz] = weight_hv;
        ++nz;
        jacobian_map->colidx[nz] = ff.index_of_u(i,j+1);
        jacobian_map->val[nz] = -weight_hv;
        ++nz;
        ++row_index;

        jacobian_map->rowptr[row_index] = nz;
        jacobian_map->colidx[nz] = ff.index_of_v(i,j);
        jacobian_map->val[nz] = weight_hv;
        ++nz;
        jacobian_map->colidx[nz] = ff.index_of_v(i,j+1);
        jacobian_map->val[nz] = -weight_hv;
        ++nz;
        ++row_index;      
      }
    }
  }

  // Add pointer to row N+1
  jacobian_map->rowptr[row_index] = nz;
}

//// static
//void ncm_flow_cost_function_sparselm::ccs_jacobian(
//  double* parameters, 
//  struct splm_ccsm* jacobian_map, 
//  int nvars, 
//  int nobs, 
//  void* adata)
//{
//  ncm_flow_cost_function_sparselm* func = static_cast<ncm_flow_cost_function_sparselm*>(adata);
//
//  // Fill the jacobian with values for the brightness constancy.
//  unsigned u_index = 0;
//  unsigned v_index = 1;
//  unsigned col_index = 0;
//  unsigned nz = 0;
//
//  const unsigned ni = func->ni_;
//  const unsigned nj = func->nj_;
//  const unsigned n_images = func->n_images_;
//
//  ncm_flow_field ff(ni, nj);
//  const double weight = (1.0 - func->alpha_);
//  const double weight_hv = weight * (n_images / 6);
//
//  const unsigned n_obs_brightness = func_->n_obs_brightness();
//
//  for (unsigned j = 0; j < nj; ++j)
//  {
//    for (unsigned i = 0; i < ni; ++i)
//    {
//      // k = 0 => dr/du, k = 1 => dr/dv
//      for (unsigned k = 0; k < 2; ++k)
//      {
//        jacobian_map->colptr[col_index] = nz;
//
//        unsigned row_index;
//
//        // Brightness penalties
//        row_index = j*ni + i;
//        for (unsigned k = 0; k < n_images; ++k)
//        {
//          jacobian_map->rowidx[nz] = row_index;
//          jacobian_map->val[nz] = func->Ix_[k](i,j);
//          ++nz;
//          ++row_index;
//        }
//
//        // Smoothness penalties
//        if (i > 0)
//        {
//          row_index = n_obs_brightness + 
//          jacobian_map->rowidx[nz] = row_index + 2*ff.index_of_u(i,j);
//          jacobian_map->val[nz] = weight_hv;
//          ++nz;
//          jacobian_map->rowidx[nz] = ff.index_of_u(i-1,j);
//          jacobian_map->val[nz] = -weight_hv;
//          ++nz;
//          ++row_index;
//        }
//
//        if (j < nj-1)
//        {
//          jacobian_map->rowidx[nz] = ff.index_of_u(i,j);
//          jacobian_map->val[nz] = weight_hv;
//          ++nz;
//          jacobian_map->rowidx[nz] = ff.index_of_u(i,j+1);
//          jacobian_map->val[nz] = -weight_hv;
//          ++nz;
//          ++row_index;
//        }
//      }
//    }
//  }
//}
 