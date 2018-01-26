#include "mex.h"
#include "math.h"
#include "string.h"
const double stdlim = 5; /* controls how many standard deviations I truncate the PSF after. */

/*  [avim,msk,M] = getAvim(biv,bih,o);
 *
 * INPUTS:
 *  biv, bih: vertical and horizontal sizes for the average image. 
 *  o: the super-res problem datastructure, which must be a MATLAB struct with the following entries:
 *    o(i).H is a 3x3 double representing the homography taking points in the i^th LR image into the SR image.
 *    o(i).la is the multiplicative photometric parameter for image i.
 *    o(i).lb is the additive photometric parameter for image i.
 *    o(i).g is the sigma for the Gaussian PSF associated with image i.
 *    o(i).im is the i^th image. W doesn't depend on image pixels, but the image size is important.
 *
 * OUTPUTS:
 *  avim: the average image (which will be of size biv-by-bih).
 *  msk: a double matrix of size biv-by-bih indicating which pixels in avim were estimated.
 *    -- average image pixels with negligible support from any of the PSFs of the low-res
 *       image pixels are returned as zeros.
 *  M: a 10xK matrix of geometric registration values.
 *
 */



void getLambdaGauss(const double *apars, const double ismv, const double ismh, double *nuv, double *nuh, double *b11, double *b12, double *b22, double *deltav, double *deltah) {
    double sqgam = apars[10];
    double denom = (ismh+1)*apars[6] + (ismv+1)*apars[7] + 1;
    *nuh = ((ismh+1)*apars[0] + (ismv+1)*apars[1]  + apars[2])/denom -1; /* the final -1 takes us back to indexed-from-0 coords. */
    *nuv = ((ismh+1)*apars[3] + (ismv+1)*apars[4]  + apars[5])/denom -1; /* the final -1 takes us back to indexed-from-0 coords. */
    denom = denom * denom;
    /* H = [h11,h12;h21,h22] is the Hessian of the transform given by apars, evaluated at
       the point that (ismv,ismh) maps to under that transform. */
    double h11 = ((apars[0]*apars[7] - apars[1]*apars[6])*(ismv+1) + apars[0] - apars[2]*apars[6])/denom;
    double h12 = ((apars[1]*apars[6] - apars[0]*apars[7])*(ismh+1) + apars[1] - apars[2]*apars[7])/denom;
    double h21 = ((apars[3]*apars[7] - apars[4]*apars[6])*(ismv+1) + apars[3] - apars[5]*apars[6])/denom;
    double h22 = ((apars[4]*apars[6] - apars[3]*apars[7])*(ismh+1) + apars[4] - apars[5]*apars[7])/denom;
    /* I also know that the covariance of the psf will be H*Sig*H', where sig was the
       original covariance in the low-res image, which has a variance of sqgam. */
    double detH = h11*h22 - h12*h21;
    detH = 1/(apars[10]*(detH*detH));
    *b11 = detH*(h21*h21+h22*h22); /* inv(H*Sig*H')(1,1) */
    *b12 = -2*detH*(h11*h21+h12*h22);/* inv(H*Sig*H')(1,2)+inv(H*Sig*H')(2,1) */
    *b22 = detH*(h11*h11+h12*h12); /* inv(H*Sig*H')(2,2) */
    /* Now find the greatest h and v extend of this PSF that's within 4*gam in the original (low-res) PSF.
     * Derive this by expressing x in terms of y, using the quadratic formula, and finding the value of y^2 necessary to set the sqrt(b^2-4ac) part to zero. */
    denom = sqrt(4*(*b11)*(*b22)-(*b12)*(*b12));
    double gam = sqrt(apars[10]);
    *deltav = 2*stdlim*gam*sqrt(*b11)/denom; /* vertical extent of the new psf kernel */
    *deltah = 2*stdlim*gam*sqrt(*b22)/denom; /* vertical extent of the new psf kernel */
    return;
}

void avimFromoN10Gauss(mxArray **avim, mxArray **mask, const double *M, const int biv, const int bih, const mxArray **o) {
    mxArray *fa;
    double wval, *im, lambdavals[]={0,0,0,0,0,0};
    double myscaler = 0, lambda, sc, detH, b11, b12, b22, scgam, nh, nv, nuh, nuv, deltav, deltah;
    int lowh, highh, lowv, highv;
    uint smv, smh, ismv, ismh, ibiv, ibih, kj;
    int K = mxGetNumberOfElements(o[0]);
    avim[0] = mxCreateDoubleMatrix(biv,bih,mxREAL);
    double *av = mxGetPr(avim[0]);
    mask[0] = mxCreateDoubleMatrix(biv,bih,mxREAL);
    double *ms = mxGetPr(mask[0]);
    double *h = (double *) mxCalloc((int)(biv*bih),sizeof(double));
    
    for (kj=0; kj<K; kj++) { /* iterate over the K low-res images */
        fa = mxGetField(o[0], kj, "im");
        im = (double *) mxGetPr(fa);
        smv = mxGetM(fa);
        smh = mxGetN(fa);
        for (ismh=0; ismh < smh; ismh++) { /* iterate over low-res image's rows. */
            for (ismv=0; ismv < smv; ismv++) { /* scan down each column in turn. */
                getLambdaGauss(&M[kj*11], ismv, ismh, &nuv, &nuh, &b11, &b12, &b22, &deltav, &deltah); /* This call maps the location to the HR frame and find the params of the affine approximation to the PSF under projection. */
                if ((nuh>=0)&&(nuh<bih)&&(nuv>=0)&&(nuv<biv)) { /* this means the mapped low-res pixel lands somewhere in the superimage. */
                    lowh = (floor(nuh-deltah)+1>0) ? (int)floor(nuh-deltah+1) : 0;
                    lowv = (floor(nuv-deltav)+1>0) ? (int)floor(nuv-deltav+1) : 0;
                    highh = (ceil(nuh+deltah)<bih) ? (int)ceil(nuh+deltah) : bih;
                    highv = (ceil(nuv+deltav)<biv) ? (int)ceil(nuv+deltav) : biv; /* Find the right box in the HR image to scan over. */
                    myscaler = 0;
                    /* Find all the values at pixel locations (and for a 1-pixel border around region). */
                    for (ibih = lowh; ibih < highh; ibih++) {
                        for (ibiv = lowv; ibiv < highv; ibiv++) {
                            nv = ibiv-nuv;
                            nh = ibih-nuh;
                            h[ibih*biv+ibiv] = exp(-0.5*(b11*nh*nh + b12*nh*nv + b22*nv*nv));
                            myscaler += h[ibih*biv+ibiv];
                        }
                    }
                    if (myscaler>0.0000001) {
                        for (ibih = lowh; ibih < highh; ibih++) {
                            for (ibiv = lowv; ibiv < highv; ibiv++) {
                                h[ibih*biv+ibiv] = h[ibih*biv+ibiv]/myscaler;
                                av[ibih*biv+ibiv] += im[ismh*smv+ismv]*h[ibih*biv+ibiv];
                                ms[ibih*biv+ibiv] += h[ibih*biv+ibiv];
                            }
                        }
                    }
                }
            }
        }
    }
    for (ibih = 0; ibih < bih; ibih++) {
        for (ibiv = 0; ibiv < biv; ibiv++) {
            if (ms[ibih*biv+ibiv]>0.0000001) {
                av[ibih*biv+ibiv] = av[ibih*biv+ibiv]/ms[ibih*biv+ibiv];
                ms[ibih*biv+ibiv] = 0;
            }
            else {
                av[ibih*biv+ibiv] = 0;
                ms[ibih*biv+ibiv] = 1;
            }
        }
    }
    mxFree(h);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int smv, smh, nzmax, K, mpt = 0;
    double *im, *H;
    const char *field_name;
    const mxArray *fa;
    uint kj, i;
    
    
   /* Check for proper number of arguments. */
    if(nrhs<3) { printf("Error: Too few inputs provided. (I should see biv, bih, o)\n"); return; }
    else if(nlhs>3) { printf("Error: Too many output arguments.\n"); return; }
    int biv = (int) mxGetPr(prhs[0])[0];
    int bih = (int) mxGetPr(prhs[1])[0];
    
   /* Now deal with input prhs[2], which is the K-element structure, o. */
    if (!mxIsStruct(prhs[2])) { printf("Input 4 must be a structure."); return; }
    K = mxGetNumberOfElements(prhs[2]);
    
    
    /* Check all the contents of the struct I need are actually there */
    const char *field_names[] = {"H", "la", "lb", "g", "im"}; /* Make sure "im" is the last one to be called! */
    for (kj=0;kj<K;kj++) {
        for (i=0; i<5; i++) {
            fa = mxGetField(prhs[2], kj, field_names[i]);
            if (fa==NULL) {
                mexPrintf("ERROR: there doesn't appear to be a \"%s\" for element %i of the input struct array.\n",field_names[i],kj+1);
                plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[0])[0]=0;
                plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[1])[0]=0;
                plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[2])[0]=0;
                return;
            }
        }
    }
    plhs[2] = mxCreateDoubleMatrix(11,K,mxREAL);
    double *M = mxGetPr(plhs[2]);
    
    /* One more loop through the inputs to start reading stuff into outputs. */
    for (kj=0;kj<K;kj++) {
        fa = mxGetField(prhs[2], kj, "H");
        H = (double *) mxGetPr(fa);
        M[mpt++] = H[0]/H[8]; /* I want M to represent row-major rather than         */
        M[mpt++] = H[3]/H[8]; /* column-major form, because that's what I            */
        M[mpt++] = H[6]/H[8]; /* wrote the getLambdaGauss function to work with..    */
        M[mpt++] = H[1]/H[8]; /*  .                                                  */
        M[mpt++] = H[4]/H[8]; /*  .                                                  */
        M[mpt++] = H[7]/H[8]; /*  .                                                  */
        M[mpt++] = H[2]/H[8]; /*  .                                                  */
        M[mpt++] = H[5]/H[8]; /* .. so these lines are doing a transpose.            */
        M[mpt++] = *mxGetPr(mxGetField(prhs[2], kj, "la"));
        M[mpt++] = *mxGetPr(mxGetField(prhs[2], kj, "lb"));
        M[mpt++] = pow(*mxGetPr(mxGetField(prhs[2], kj, "g")),2);
    }
    
    avimFromoN10Gauss(&plhs[0], &plhs[1], M, biv, bih, &prhs[2]);
    return;
}
