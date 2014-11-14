#include "ncm_flow_minimizer.h"

#include <cholmod.h>

void ncm_flow_minimizer::minimize()
{
  cholmod_common c;
  cholmod_start(&c);

  // repeat until converged
  //   compute J
  //   compute JtJ

  //cholmod_sparse* A;
  //A = cholmod_read_sparse(stdin, &c);
  //cholmod_print_sparse(A, "A", &c);

  ////   compute Jte
  //cholmod_dense* Jte;
  //Jte = cholmod_ones(JtJ->nrow, 1, JtJ->xtype, &c);

  //   multiply diagonal of JtJ by (1+lambda)
  //   solve (JtJ + ...)d = Jte for d
  //   check for convergence

  cholmod_finish(&c);
}

//void ncm_flow_minimizer::evaluate(cost_function)
//{
//  if (JtJ == NULL || JtJ->stype == 0)		    /* A must be symmetric */
//  {
//    return;
//  }
//
//  // Factorize JtJ into LLt
//  cholmod_factor* L;
//  L = cholmod_analyze(JtJ, &c);
//  cholmod_factorize(JtJ, L, &c);
//
//  // Solve for delta_p
//  cholmod_dense* delta_p;
//  delta_p = cholmod_solve(CHOLMOD_A, L, Jte, &c);
//
//  // Compute the new e(psilon)
//  cholmod_dense* r;
//  r = cholmod_copy_dense(Jte, &c);
//
//  /* r = r-Ax */
//  double one [2] = {1,0}, m1 [2] = {-1,0} ;	    /* basic scalars */
//  cholmod_sdmult(JtJ, 0, m1, one, delta_p, r, &c);
//
//  /* print norm(r) */
//  printf ("norm(Jt.e - JtJ.delta_p) %8.1e\n",
//          cholmod_norm_dense(r, 0, &c)); 
//
//  /* free matrices */
//  cholmod_free_factor(&L, &c);
//  cholmod_free_sparse(&JtJ, &c);
//  cholmod_free_dense(&r, &c);
//  cholmod_free_dense(&x, &c);
//  cholmod_free_dense(&b, &c);
//}
