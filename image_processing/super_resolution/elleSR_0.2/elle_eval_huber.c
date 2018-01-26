#include "mex.h"

/* eout = elle_eval_huber(x,alp,bet,biv,bih); */

double alp;

double huber(double x) {
    if (x<0) {x=-x;}
    if (x<alp) { return (x*x); }
    else { return (2*x*alp - alp*alp); } 
    return -1;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   const double r2 = 1.414213562373095; /* ideally this would be sqrt(2), but the compiler is complaining.. */
   uint i=0, j=0;
   
   /* Check for proper number of arguments. */
   if(nrhs!=5) { 
       printf("Wrong number of inputs (Should have: x, alpha, beta, biv, bih)\n");
       printf("I see %i, and should see 5.\n",nrhs);
       return;
   }
   if(nlhs>1) { printf("Too many output arguments"); return; }
   
   if(mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) { printf("X (1st arg) must be of type double."); return;}
   if(mxGetClassID(prhs[1])!=mxDOUBLE_CLASS) { printf("alpha (2nd arg) must be of type double."); return;}
   if(mxGetClassID(prhs[2])!=mxDOUBLE_CLASS) { printf("alpha (3rd arg) must be of type double."); return;}
   if(mxGetClassID(prhs[3])!=mxDOUBLE_CLASS) { printf("biv (4th arg) must be of type double."); return;}
   if(mxGetClassID(prhs[4])!=mxDOUBLE_CLASS) { printf("bih (5th arg) must be of type double."); return;}
   if(mxGetM(prhs[0])>1) {printf("X must be a row vector."); return;}
   
   double *X = (double *) mxGetPr(prhs[0]);
   alp = (double) mxGetPr(prhs[1])[0];
   double bet = (double) mxGetPr(prhs[2])[0];
   int biv = (int) mxGetPr(prhs[3])[0];
   int bih = (int) mxGetPr(prhs[4])[0];

   plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
   double *eout = mxGetPr(plhs[0]);
   
   double etemp = 0;
   int c;
   for (i=0; i<biv; i++) {
       for (j=0; j<bih; j++) {
           c = j*biv+i;
           if (i<(biv-1)) { *eout += huber(X[c+1]-X[c]);
           if (j<(bih-1)) { *eout += huber((X[c+biv+1]-X[c])/r2); }
           }
           if (j<(bih-1)) { *eout += huber(X[c+biv]-X[c]);
           if (i>=1) { *eout += huber((X[c+biv-1]-X[c])/r2); }
           }
       }
   }
   *eout *= bet/2; /* now eout = nu*Huber(<x_neighbour_differences,alp>); */
}

