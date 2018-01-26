#include "mex.h"

/* g = elle_eval_huber_grad(x,alp,bet,biv,bih); */

double alp;

double hubergradx(const double x) {
    if (x>0) {
        if (x<alp) { return x; }
        else {return alp; }
    }
    else {
        if (x>-alp) { return x; }
        else {return -alp; }
    }
    return 0;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   const double r2 = 1.414213562373095; /* ideally this would be sqrt(2), but the compiler is complaining.. */
   unsigned i=0, j=0;
   
   /* Check for proper number of arguments. */
   if(nrhs!=5) {
       printf("Wrong number of inputs (Should have: X, alpha, beta, biv, bih)\n");
       printf("I see %i, and should see 5.\n",nrhs);
       return;
   }
   if(nlhs>1) { printf("Too many output arguments\n"); return; }
   
   if(mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) { printf("X (1st arg) must be of type double.\n"); return; }
   if(mxGetClassID(prhs[1])!=mxDOUBLE_CLASS) { printf("alpha (2nd arg) must be of type double."); return;}
   if(mxGetClassID(prhs[2])!=mxDOUBLE_CLASS) { printf("beta (3rd arg) must be of type double."); return;}
   if(mxGetClassID(prhs[3])!=mxDOUBLE_CLASS) { printf("biv (4th arg) must be of type double."); return;}
   if(mxGetClassID(prhs[4])!=mxDOUBLE_CLASS) { printf("bih (5th arg) must be of type double."); return;}
   
   double *X = (double *) mxGetPr(prhs[0]);
   alp = (double) mxGetPr(prhs[1])[0];
   double bet = (double) mxGetPr(prhs[2])[0];
   int biv = (int) mxGetPr(prhs[3])[0];
   int bih = (int) mxGetPr(prhs[4])[0];
   
   plhs[0] = mxCreateDoubleMatrix(1,biv*bih,mxREAL);
   double *G = mxGetPr(plhs[0]);
   
    /* huber_prior_gradx(G, X, biv, bih, nu) */
    /* adds the gradient of p(x) wrt x to the current values of the elements of G. */
   double *Y2 = (double *) mxCalloc(biv*bih*4,sizeof(double)); /* This is initialised to all zeros. */
   int oset1 = biv*bih;
   int oset2 = 2*oset1;
   int oset3 = oset1+oset2;
   int c;
   for (i=0; i<biv; i++) {
       for (j=0; j<bih; j++) {
           c = j*biv+i;
           if (i<(biv-1)) { Y2[c] = hubergradx(-X[c]+X[c+1]);
           if (j<(bih-1)) {  Y2[c+oset2] = hubergradx((-X[c]+X[c+biv+1])/r2); }
           }
           if (j<(bih-1)) { Y2[c+oset1] = hubergradx(-X[c]+X[c+biv]);
           if (i>=1) {  Y2[c+oset3] = hubergradx((-X[c]+X[c+biv-1])/r2); }
           }
       }
   }
   /* loop to go back from vectors of size (4n)*(n) to vectors of size n*n. */
   double tsum = 0;
   for (i=0; i<biv; i++) {
       for (j=0; j<bih; j++) {
           c = j*biv+i;
           tsum = 0;
           if (i>0) {tsum = tsum + Y2[c-1];}
           if (j>0) {tsum = tsum + Y2[c-biv+oset1];}
           if (i>0&&j>0) {tsum = tsum + Y2[c-biv-1+oset2]/r2;}
           if (i<(biv-1)&&j>0) {tsum = tsum + Y2[c-biv+1+oset3]/r2;}
           
           if (i<(biv-1)) {tsum = tsum - Y2[c];}
           if (j<(bih-1)) {tsum = tsum - Y2[c+oset1];}
           if (j<(bih-1)&&i<(biv-1)) {tsum = tsum - Y2[c+oset2]/r2;}
           if (j<(bih-1)&&i>0) {tsum = tsum - Y2[c+oset3]/r2;}
           
           G[c] += bet*tsum;
       }
   }
    mxFree(Y2);
}
