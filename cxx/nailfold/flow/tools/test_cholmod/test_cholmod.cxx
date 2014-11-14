#include <cholmod.h>

int main(int argc, char* argv[])
{
  cholmod_common c ;
  cholmod_start(&c);			    /* start CHOLMOD */
  
  cholmod_sparse* A;
  A = cholmod_read_sparse (stdin, &c) ;	    /* read in a matrix */
  cholmod_print_sparse (A, "A", &c) ;		    /* print the matrix */

  if (A == NULL || A->stype == 0)		    /* A must be symmetric */
  {
    cholmod_free_sparse (&A, &c) ;
    cholmod_finish (&c) ;
    return (0) ;
  }

  cholmod_dense *b;
  b = cholmod_ones (A->nrow, 1, A->xtype, &c) ;   /* b = ones(n,1) */

  cholmod_factor *L ;
  L = cholmod_analyze (A, &c) ;		    /* analyze */
  cholmod_factorize (A, L, &c) ;		    /* factorize */

  /* solve Ax=b */
  cholmod_dense *x;
  x = cholmod_solve (CHOLMOD_A, L, b, &c) ;	    

  /* r = b */
  cholmod_dense *r ;
  r = cholmod_copy_dense (b, &c) ;

  /* r = r-Ax */
  double one [2] = {1,0}, m1 [2] = {-1,0} ;	    /* basic scalars */
  cholmod_sdmult (A, 0, m1, one, x, r, &c);

  /* print norm(r) */
  printf ("norm(b-Ax) %8.1e\n",
          cholmod_norm_dense (r, 0, &c)); 

  /* free matrices */
  cholmod_free_factor (&L, &c) ;		    
  cholmod_free_sparse (&A, &c) ;
  cholmod_free_dense (&r, &c) ;
  cholmod_free_dense (&x, &c) ;
  cholmod_free_dense (&b, &c) ;
  cholmod_finish (&c) ;			    /* finish CHOLMOD */

  return 0;
}