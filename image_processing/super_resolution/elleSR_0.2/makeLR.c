#include "mex.h"
#include "math.h"
#include "string.h"
#include "nrand.h"

/* [im,noise,orig] = makeLR(biv,bih,o); 
 *
 * INPUTS:
 *  biv, bih: vertical and horizontal sizes for the high-res image. 
 *  o: the super-res problem datastructure, which must be a MATLAB struct with the following entries:
 *    o.H is a 3x3 double representing the homography taking points in the LR image into the SR image.
 *    o.la is the multiplicative photometric parameter.
 *    o.lb is the additive photometric parameter.
 *    o.g is the sigma for the Gaussian PSF.
 *    o.n is the sigma for the additive Gaussian iid noise on each pixel.
 *    o.v vertical size for the image.
 *    o.h horizontal size for the image.
 *  NOTE: if o is larger than 1x1, only the first part is used (i.e. this function doens't 
 *        currently produce a stack of images, just one at a time).
 *
 * OUTPUTS:
 *  im: the generated low-resolution image, including image noise.
 *  noise: the noise image.
 *  orig: the low-resolution image without the added noise.
 */


void getLambdaGauss(const double *apars, const double sqgam, const double ismv, const double ismh, double *lambda, double *nuv, double *nuh, double *b11, double *b12, double *b22) {
    double kdenom = (ismh+1)*apars[6] + (ismv+1)*apars[7] + 1;
    *nuh = ((ismh+1)*apars[0] + (ismv+1)*apars[1]  + apars[2])/kdenom -1; /* the final -1 takes us back to indexed-from-1 coords. */
    *nuv = ((ismh+1)*apars[3] + (ismv+1)*apars[4]  + apars[5])/kdenom -1; /* the final -1 takes us back to indexed-from-1 coords. */
    kdenom = kdenom * kdenom;
    /* H = [h11,h12;h21,h22] is the Hessian of the transform given by apars, evaluated at
       the point (ismv,ismh) maps to under that transform. */
    double h11 = ((apars[0]*apars[7] - apars[1]*apars[6])*(ismv+1) + apars[0] - apars[2]*apars[6])/kdenom;
    double h12 = ((apars[1]*apars[6] - apars[0]*apars[7])*(ismh+1) + apars[1] - apars[2]*apars[7])/kdenom;
    double h21 = ((apars[3]*apars[7] - apars[4]*apars[6])*(ismv+1) + apars[3] - apars[5]*apars[6])/kdenom;
    double h22 = ((apars[4]*apars[6] - apars[3]*apars[7])*(ismh+1) + apars[4] - apars[5]*apars[7])/kdenom;
    /* I also know that the covariance of the psf will be H*Sig*H', where sig was the
       original covariance in the low-res image, which has a variance of sqgam. */
    double detH = h11*h22 - h12*h21;
    detH = 1/(sqgam*(detH*detH));
    *b11 = detH*(h21*h21+h22*h22); /* inv(H*Sig*H')(1,1) */
    *b12 = -2*detH*(h11*h21+h12*h22);/* inv(H*Sig*H')(1,2)+inv(H*Sig*H')(2,1) */
    *b22 = detH*(h11*h11+h12*h12); /* inv(H*Sig*H')(2,2) */
    /* The scaling along the semimajor axis of the new covariance ellipse will be given by the sqrt of the
       bigger eigenvalue of inv(Sig). I'm going to call that scaling 'lambda' in here. */
    *lambda = *b11+*b22;
    *lambda = 0.5*(*lambda - sqrt((*lambda)*(*lambda) - 4*((*b11)*(*b22)-0.25*(*b12)*(*b12))));
    *lambda = 3/sqrt(*lambda) + 1; /* Multiplying by 3 means everything with P>tiny lies within this distance of the transformed mu. */
    return;
}


/* getY gets the down-projected superim Yd. Useful when using the generative model in the forward direction. */
void getYgauss(double *Yd, const int biv, const int bih, const int smv, const int smh, const double sqgam, const double* Hph, const double* x) {
    double nuv, nuh, myscaler, wval, nv, nh, lambda, b11, b12, b22, eout=0;
    int lowh, highh, lowv, highv;
    uint ismh, ismv, ibih, ibiv;
    for (ismh=0; ismh < smh; ismh++) {
        for (ismv=0; ismv < smv; ismv++) {
            Yd[ismh*smv+ismv] = 0;
            getLambdaGauss(Hph, sqgam, ismv, ismh, &lambda, &nuv, &nuh, &b11, &b12, &b22);
            if ((nuh>=0)&&(nuh<bih)&&(nuv>=0)&&(nuv<biv)) { /* this means the mapped low-res pixel lands somewhere in the superimage. */
                lowh = (floor(nuh-lambda)>0) ? (int)floor(nuh-lambda) : 0;
                lowv = (floor(nuv-lambda)>0) ? (int)floor(nuv-lambda) : 0;
                highh = (ceil(nuh+lambda)<bih) ? (int)ceil(nuh+lambda) : bih;
                highv = (ceil(nuv+lambda)<biv) ? (int)ceil(nuv+lambda) : biv;
                myscaler = 0;
                /* Find all the values at pixel locations (and for a 1-pixel border around region). */
                for (ibih = lowh; ibih < highh; ibih++) {
                    for (ibiv = lowv; ibiv < highv; ibiv++) {
                        nv = ibiv-nuv;
                        nh = ibih-nuh;
                        wval = exp(-0.5*(b11*nh*nh + b12*nh*nv + b22*nv*nv));
                        myscaler += wval;
                        Yd[ismh*smv+ismv] += wval*x[ibih*biv+ibiv];
                    }
                }
                if (myscaler > 0.0000001) { Yd[ismh*smv+ismv] = (Hph[9]*(Yd[ismh*smv+ismv]/myscaler) + Hph[8]); }
                else{
                    printf("WARNING! Insufficient support for psf!\n"); 
                    printf("( %i , %i ) -> ( %g , %g ); lambda = %g \n",ismv,ismh,nuv,nuh,lambda);
                    printf("v %i to %i; h %i to %i \n",lowv,highv,lowh,highh);
                    printf("b = [ %g , %g , %g ] \n",b11,b12,b22);
                    return;
                }
            }
        }
    }
    return;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
    
   /* Check for proper number of arguments. */
    if(nrhs!=2) {
        printf("Wrong number of inputs (Should have:gtruth,o)\n");
        printf("I see %i, and should see 2.\n",nrhs);
        return;
    }
    if(nlhs!=3) { printf("Wrong number of outputs. I should see 3 (im, noise, orig), but I see %i.\n",nlhs); return; }
    if(mxGetClassID(prhs[0])!=mxDOUBLE_CLASS) { printf("X (1st arg) must be of type double."); return;}
    /* Now deal with input prhs[1], which is the K-element structure, o. */
    if (!mxIsStruct(prhs[1])) { printf("o (2nd arg) must be a structure."); return; }
    
    /* Grab inputs */
    double *X = (double *) mxGetPr(prhs[0]);
    int biv = (int) mxGetM(prhs[0]);
    int bih = (int) mxGetN(prhs[0]);;
    
    
    /* Check all the contents of the struct I need are actually there */
    mxArray *fa;
    const char *field_names[] = {"H", "la", "lb", "g", "n", "v", "h"};
    uint i;
    for (i=0; i<7; i++) {
        fa = mxGetField(prhs[1], 0, field_names[i]);
        if (fa==NULL) {
            mexPrintf("ERROR: there doesn't appear to be a \"%s\" element in the input struct.\n",field_names[i]);
            plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[0])[0]=0;
            plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[1])[0]=0;
            plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[2])[0]=0;
            return;
        }
    }
    
    double *H = mxGetPr(mxGetField(prhs[1], 0, "H"));
    int smv = (int)mxGetPr(mxGetField(prhs[1], 0, "v"))[0];
    int smh = (int)mxGetPr(mxGetField(prhs[1], 0, "h"))[0];
    double sqgam = pow(mxGetPr(mxGetField(prhs[1], 0, "g"))[0],2);
    double std = mxGetPr(mxGetField(prhs[1], 0, "n"))[0];
    
    /* Create output arrays */
    plhs[0] = mxCreateDoubleMatrix(smv,smh,mxREAL);
    double *im = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(smv,smh,mxREAL);
    double *noise = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(smv,smh,mxREAL);
    double *orig = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(1,10,mxREAL);
    double *M = (double *) mxCalloc(10,sizeof(double)); /* This is initialised to all zeros.  */
    M[0] = H[0]/H[8]; /* I want M to represent row-major rather than        */
    M[1] = H[3]/H[8]; /* column-major form, because that's what I           */
    M[2] = H[6]/H[8]; /* wrote the getLambdaGauss function to work with..   */
    M[3] = H[1]/H[8]; /*  .                                                 */
    M[4] = H[4]/H[8]; /*  .                                                 */
    M[5] = H[7]/H[8]; /*  .                                                 */
    M[6] = H[2]/H[8]; /*  .                                                 */
    M[7] = H[5]/H[8]; /* .. so these lines are doing a transpose.           */
    M[8] = mxGetPr(mxGetField(prhs[1], 0, "lb"))[0];
    M[9] = mxGetPr(mxGetField(prhs[1], 0, "la"))[0];
    
    long idum = -55;
    getYgauss(orig, biv, bih, smv, smh, sqgam, M, X);
    for (i=0; i<smv*smh; i++) {
        noise[i] = gasdev(&idum)*std;
        im[i] = noise[i] + orig[i];
    }
    mxFree(M);
    
    return;
}
