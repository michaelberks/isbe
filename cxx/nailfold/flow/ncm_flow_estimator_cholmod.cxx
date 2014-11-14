#include "ncm_flow_estimator_cholmod.h"

#include <nailfold/flow/ncm_flow_field.h>

#include <cholmod.h>
#include <cholmod_internal.h>
 
//
//  Public methods
//
 
//: Estimate the flow field from the input image stack.
void ncm_flow_estimator_cholmod::estimate(
  ncm_flow_field& flow_field,
  bool update)
{
  ncm_flow_cost_function_cholmod func(
    image_stack_,
    warped_image_stack_,
    robustifier_);

  const unsigned n_parameters = flow_field.n_parameters();
  vnl_vector<double> parameter_vector(n_parameters, 0.0);
  //flow_field.get_to_vector(parameter_vector.data_block());

  // Determine the optimal parameters
  minimize(func, parameter_vector);

  if (update)
    flow_field.update_from_vector(parameter_vector.data_block());
  else
    flow_field.set_from_vector(parameter_vector.data_block());
}

//
//  Private methods
//
 
//: 
void ncm_flow_estimator_cholmod::minimize(
  ncm_flow_cost_function_cholmod func,
  vnl_vector<double>& parameters)
{
  double const Jte_mag_threshold = 1e-3;
  double const delta_p_threshold = 1e-3;
  double const d_error_threshold = 1e-3;

  cholmod_common c;
  set_cholmod_options(c);

  cholmod_start(&c);

  vnl_vector<double> observed(func.n_observations(), 0.0);
  vnl_vector<double> lowest_error = observed - func.f(parameters);
  double lowest_sse = lowest_error.squared_magnitude();

  double lambda = 1.0;
  double const lambda_multiplier = 10.0;
  assert(lambda_multiplier > 1.0);
  
  //vnl_sparse_matrix<double> J;
  //vnl_vector<double> Jte;
  //vnl_sparse_matrix<double> JtJ;

  const unsigned n_its_max = 5;
  for (unsigned it = 0; it < n_its_max; ++it)
  {
    // Compute Jacobian
    vnl_sparse_matrix<double> J;
    J = func.jacobian(parameters);

    // Compute Jte (and its cholmod counterpart)
    // (Jte is the projection of the vector from f(p) to obs, onto the
    // column space of J. This will be zero when the vector is normal to the 
    // space spanned by J, which occurs when we're as close as we can get.)

    vnl_vector<double> Jte;
    J.transpose().mult(observed - func.f(parameters), Jte);
    if (Jte.squared_magnitude() < Jte_mag_threshold)
      break;
    cholmod_dense* Jte_cm = NULL;
    convert_vector_to_dense(Jte, Jte_cm, &c);

    // Compute JtJ and its diagonal.
    vnl_sparse_matrix<double> JtJ;
    JtJ = J.transpose() * J;
    vnl_vector<double> JtJ_diag;
    JtJ.diag_AtA(JtJ_diag);

    // Search for lambda (within reason) that gives best results.
    // The higher the lambda, the more the solution for delta_p will be biased
    // toward zero.
    bool lambda_is_decreasing = false;
    vnl_vector<double> delta_p;
    vnl_vector<double> better_parameters;
    double d_error = 0.0;
    while ((1e-6 <= lambda) && (lambda <= 1e6))
    {
      // Multiply diagonal of JtJ by (1+lambda) for preconditioning
      vnl_vector<double> lambda_JtJ_diag = lambda * JtJ_diag;

      for (unsigned i = 0; i < JtJ.rows(); ++i)
        JtJ(i,i) += lambda_JtJ_diag[i];

      solve_system(JtJ, Jte_cm, c, delta_p);

      // Check the SSE corresponding to f(p+dp).
      vnl_vector<double> err = observed - func.f(parameters + delta_p);
      double eps = err.squared_magnitude();

      if (eps < lowest_sse)
      {
        lambda /= lambda_multiplier;

        lambda_is_decreasing = true;
        better_parameters = parameters + delta_p;
        d_error += lowest_sse - eps;
        lowest_sse = eps;
      }
      else
      {
        lambda *= lambda_multiplier;

        if (lambda_is_decreasing)
          break;
      }

      // Restore JtJ to its original value.
      for (unsigned i = 0; i < JtJ.rows(); ++i)
        JtJ(i,i) -= lambda_JtJ_diag[i];
    }

    cholmod_free_dense(&Jte_cm, &c);

    if (lambda < 1e-6)
      lambda = 1e-6;
    else if (lambda > 1e6)
      lambda = 1e6;

    if (!better_parameters.empty())
    {
      delta_p = better_parameters - parameters;
      parameters = better_parameters;
    }

    // Finish if the change in parameters is below a threshold.
    if (delta_p.squared_magnitude() < delta_p_threshold)
      break;

    // Finish if the reduction in squared error is below a threshold.
    if (vcl_abs(d_error) < d_error_threshold)
      break;
  }

  cholmod_finish(&c);
}
 
//:
bool ncm_flow_estimator_cholmod::converged()
{
  return false;
}
 
//:
void ncm_flow_estimator_cholmod::convert_matrix_to_sparse(
  const vnl_sparse_matrix<double> in,
  cholmod_sparse*& out,
  cholmod_common* common)
{
  // Count the number of nonzeros.
  unsigned nnz = 0;
  in.reset();
  while (in.next())
    ++nnz;
  
  if (out == NULL)
  {
    out = cholmod_allocate_sparse(
      in.rows(), in.cols(), nnz,
      /* sorted = */ true,
      /* packed = */ true,
      /* stype = */ 1,
      /* xtype = */ CHOLMOD_REAL,
      common);
  }

  // Copy data from vnl_sparse_matrix.
  double* xptr = static_cast<double*>(out->x);
  Int* iptr = static_cast<Int*>(out->i);
  Int* pptr = static_cast<Int*>(out->p);

  unsigned ind = 0;
  unsigned col = 0;
  in.reset();
  while (in.next())
  {
    // This looks wrong but isn't: we're reading values from a CRS matrix and 
    // writing to a CCS matrix. Since the matrix is symmetric, it's ok to
    // copy directly from the rows of one to the columns of the other.

    for (/* no init */; col < in.getrow(); ++col)
    {
      pptr[col+1] = ind;
    }
    iptr[ind] = in.getcolumn(); // row index of in
    xptr[ind] = in.value();

    ++ind;
  }

  // Fill in any remaining entries in pptr
  for (/* no init */; col < in.cols(); ++col)
  {
    pptr[col+1] = ind;
  }
}
 
//:
void ncm_flow_estimator_cholmod::convert_vector_to_dense(
  const vnl_vector<double> in,
  cholmod_dense*& out,
  cholmod_common* common)
{
  out = cholmod_allocate_dense(
    in.size(), 1,
    /* d = */ in.size(),
    /* xtype = */ CHOLMOD_REAL,
    common);

  // Copy data from vnl_sparse_matrix.
  double* xptr = static_cast<double*>(out->x);

  for (unsigned i = 0; i < in.size(); ++i)
    xptr[i] = in[i];
}
 
//:
void ncm_flow_estimator_cholmod::set_cholmod_options(
  cholmod_common& c)
{
  //c.nmethods = 1 ;
  //c.method[0].ordering = CHOLMOD_AMD;
}
 
//:
void ncm_flow_estimator_cholmod::solve_system(
  // Input
  vnl_sparse_matrix<double> const& JtJ,
  cholmod_dense* Jte_cm,
  // Input/Output
  cholmod_common& c,
  // Output
  vnl_vector<double>& delta_p)
{
  // Factorize JtJ into LLt
  cholmod_sparse* JtJ_cm = NULL;
  convert_matrix_to_sparse(JtJ, JtJ_cm, &c);
  cholmod_factor* L = cholmod_analyze(JtJ_cm, &c);
  cholmod_factorize(JtJ_cm, L, &c);
  cholmod_free_sparse(&JtJ_cm, &c);
  
  // Solve for delta_p
  cholmod_dense* delta_p_cm = cholmod_solve(CHOLMOD_A, L, Jte_cm, &c);
  cholmod_free_factor(&L, &c);

  // Copy delta_p_cm into a vnl_vector
  delta_p.set_size(delta_p_cm->nrow);
  double* xptr = static_cast<double*>(delta_p_cm->x);
  for (unsigned i = 0; i < delta_p.size(); ++i)
    delta_p[i] = xptr[i];
  cholmod_free_dense(&delta_p_cm, &c);

  // Jte_cm is free()d outside of this function.
}
