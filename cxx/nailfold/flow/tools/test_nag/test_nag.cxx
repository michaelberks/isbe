#include <vcl_iosfwd.h>
#include <vcl_cstddef.h>
#include <vcl_iostream.h>
#include <vcl_vector.h>
#include <vcl_sstream.h>
#include <vcl_cmath.h>

#include <mbl/mbl_parse_keyword_list.h>

// NAG libraries
#include <nag.h>
#include <nag_stdlib.h>
#include <nage04.h>
#include <nagx02.h>
 
#define DIV(n, m) ((n) / (m))
//#define MAX(n, m) (((n) >= (m)) ? (n) : (m))

/* Chained Rosenbrock function */
// initial point: p[i] =  -1.2 for i odd, 1.0 for i even, minimum at (1, 1, ..., 1)
void chainedRosenbrock(
  const double *p,
  double *hx, 
  int nParameters, int nResiduals,
  void *adata)
{
register int k, k1, i;

  for(k=0; k<nResiduals; ++k){
    k1=k+1; // k is zero-based, convert to one-based
    i=DIV(k1+1, 2) - 1; // convert i to zero-based
    if(k1%2==1) // k1 odd
      hx[k]=10.0*(p[i]*p[i]-p[i+1]);
    else // k1 even
      hx[k]=p[i]-1.0;
  }
}
//
///* analytic CRS Jacobian for chainedRosenbrock() */
//void chainedRosenbrock_anjacCRS(
//  double *p, 
//  struct splm_crsm *jac, 
//  int m, int n, 
//  void *adata)
//{
//  register int k, k1, i;
//  int L;
//
//  for(k=L=0; k<n; ++k)
//  {
//    jac->rowptr[k] = L;
//
//    k1 = k+1; // k is zero-based, convert to one-based
//    i = DIV(k1+1, 2) - 1; // convert i to zero-based
//
//    if(k1%2 == 1)
//    { // k1 odd, hx[k]=10*(p[i]*p[i]-p[i+1])
//      jac->val[L] = 20.0*p[i]; 
//      jac->colidx[L] = i;
//      ++L;
//
//      jac->val[L] = -10.0; 
//      jac->colidx[L] = i+1;
//      ++L;
//    }
//    else 
//    { // k1 even, hx[k]=p[i]-1.0
//      jac->val[L] = 1.0; 
//      jac->colidx[L] = i;
//      ++L;
//    }
//  }
//
//  jac->rowptr[n]=L;
//}
//
///* analytic CCS Jacobian for chainedRosenbrock() */
//void chainedRosenbrock_anjacCCS(double *p, struct splm_ccsm *jac, int m, int n, void *adata)
//{
//register int k, k1, i;
//int l;
//struct splm_crsm jac_crs;
//
//  /* allocate and fill-in a CRS Jacobian... */
//  splm_crsm_alloc(&jac_crs, n, m, jac->nnz);
//
//  for(k=l=0; k<n; ++k){
//    jac_crs.rowptr[k]=l;
//    k1=k+1; // k is zero-based, convert to one-based
//    i=DIV(k1+1, 2) - 1; // convert i to zero-based
//    if(k1%2==1){ // k1 odd, hx[k]=10*(p[i]*p[i]-p[i+1])
//      jac_crs.val[l]=20.0*p[i]; jac_crs.colidx[l++]=i;
//      jac_crs.val[l]=-10.0; jac_crs.colidx[l++]=i+1;
//    }
//    else { // k1 even, hx[k]=p[i]-1.0
//      jac_crs.val[l]=1.0; jac_crs.colidx[l++]=i;
//    }
//  }
//  jac_crs.rowptr[n]=l;
//
//  /* ...convert to CCS */
//  splm_crsm2ccsm(&jac_crs, jac);
//  splm_crsm_free(&jac_crs);
//}
//
///* analytic CCS Jacobian for chainedRosenbrock(). Jacobian is first constructed using sparse triplet format */
//void chainedRosenbrock_anjacCCS_ST(double *p, struct splm_ccsm *jac, int m, int n, void *adata)
//{
//register int k, k1, i;
//struct splm_stm jac_st;
//
//  /* allocate and fill-in a ST Jacobian... */
//  splm_stm_allocval(&jac_st, n, m, jac->nnz);
//
//  for(k=0; k<n; ++k){
//    k1=k+1; // k is zero-based, convert to one-based
//    i=DIV(k1+1, 2) - 1; // convert i to zero-based
//    if(k1%2==1){ // k1 odd, hx[k]=10*(p[i]*p[i]-p[i+1])
//      splm_stm_nonzeroval(&jac_st, k, i, 20.0*p[i]);
//      splm_stm_nonzeroval(&jac_st, k, i+1, -10.0);
//    }
//    else { // k1 even, hx[k]=p[i]-1.0
//      splm_stm_nonzeroval(&jac_st, k, i, 1.0);
//    }
//  }
//  /* ...convert to CCS */
//  splm_stm2ccsm(&jac_st, jac);
//  splm_stm_free(&jac_st);
//}
//
///* zero pattern of CCS Jacobian for chainedRosenbrock(). Jacobian is first constructed using sparse triplet format */
//void chainedRosenbrock_zpjacCCS_ST(double *p, struct splm_ccsm *jac, int m, int n, void *adata)
//{
//  register int k, k1, i;
//  struct splm_stm jac_st;
//
//  /* allocate and fill-in a ST Jacobian... */
//  splm_stm_alloc(&jac_st, n, m, jac->nnz);
//
//  for(k=0; k<n; ++k)
//  {
//    k1=k+1; // k is zero-based, convert to one-based
//    i=DIV(k1+1, 2) - 1; // convert i to zero-based
//
//    if(k1%2==1)
//    { 
//      // k1 odd, hx[k]=10*(p[i]*p[i]-p[i+1])
//      splm_stm_nonzero(&jac_st, k, i);
//      splm_stm_nonzero(&jac_st, k, i+1);
//    }
//    else 
//    { 
//      // k1 even, hx[k]=p[i]-1.0
//      splm_stm_nonzero(&jac_st, k, i);
//    }
//  }
//
//  /* ...convert to CCS */
//  splm_stm2ccsm(&jac_st, jac);
//  splm_stm_free(&jac_st);
//}



// Wrappers
void NAG_CALL lsqfun(
  Integer nResiduals, Integer nParameters, 
  const double x[], 
  double fvec[],
  Nag_Comm *comm)
{
  chainedRosenbrock(x, fvec, nParameters, nResiduals, NULL);
}

//: Main function
int main( int argc, char* argv[] )
{
  Nag_E04State *state = NULL;
  NagError fail;

  INIT_FAIL(fail);
  
  const int nParameters = 160;
  const int nResiduals = 2*(nParameters-1);
  const int nnz = 3*nParameters/2;
  const int tdfjac = nParameters;

  double* p = NULL;
  double* fvec = NULL;
  double* fjac = NULL;
  double fsumsqr = 0.0;

  if ((nResiduals >= 1) && 
      (nParameters <= nResiduals))
  {
    if (!(fjac = NAG_ALLOC(nResiduals*nParameters, double)) ||
        !(fvec = NAG_ALLOC(nResiduals, double)) ||
        !(p = NAG_ALLOC(nParameters, double)))
    {
      printf("Allocation failure\n");
      return 1;
    }
    else
    {
      vcl_cout << "Allocation success" << vcl_endl;
    }
  }

  for(unsigned i = 0; i < nParameters; i += 2)
  { 
    p[i] = -1.2;
    p[i+1] = 1.0;
  }

  Nag_E04_Opt options;
  nag_opt_init(&options); /* Initialise options structure */
  options.optim_tol = 10.0 * vcl_sqrt(nag_machine_precision);
  options.print_level = Nag_NoPrint;
  options.max_iter = 100;

  Nag_Comm comm;

  nag_opt_lsq_no_deriv(
    nResiduals, nParameters, lsqfun,
    p, &fsumsqr, fvec, 
    fjac, tdfjac, 
    &options, &comm, &fail);

  for (unsigned i = 0; i < nParameters; ++i)
    vcl_cout << p[i] << " ";
  vcl_cout << vcl_endl;

  return 0;
}
