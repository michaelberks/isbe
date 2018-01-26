#include "mex.h"
#include "math.h"
#include "string.h"
const double stdlim = 5; /* controls how many standard deviations I truncate the PSF after. */

/* [W,Y,La,Lb,M] = makeW(biv,bih,o); 
 *
 * INPUTS:
 *  biv, bih: vertical and horizontal sizes for the high-res image. 
 *  o: the super-res problem datastructure, which must be a MATLAB struct with the following entries:
 *    o(i).H is a 3x3 double representing the homography taking points in the i^th LR image into the SR image.
 *    o(i).la is the multiplicative photometric parameter for image i.
 *    o(i).lb is the additive photometric parameter for image i.
 *    o(i).g is the sigma for the Gaussian PSF associated with image i.
 *    o(i).im is the i^th image. W doesn't depend on image pixels, but the image size is important.
 *
 * OUTPUTS:
 *  W: the super-res system matrix. Note that this is the TRANSPOSE of most of my papers' "W" matrices, since it's
 *     rather easier to compute this way.
 *  Y: A KMx1 vector containing all the low-res pixel values, where "KM" is the total number of low-res pixels in the data.
 *  La: a KMx1 vector of multiplicative photometric params (yes, lots of repeated values).
 *  Lb: a KMx1 vector of additive photometric params.
 *  M: a 10xK matrix of geometric registration values.
 */



void getLambdaGauss(const double *apars, const double ismv, const double ismh, double *nuv, double *nuh, double *b11, double *b12, double *b22, double *deltav, double *deltah) {
    double sqgam = apars[10];
    double denom = (ismh+1)*apars[6] + (ismv+1)*apars[7] + 1;
    *nuh = ((ismh+1)*apars[0] + (ismv+1)*apars[1]  + apars[2])/denom -1; /* the final -1 takes us back to indexed-from-0 coords. */
    *nuv = ((ismh+1)*apars[3] + (ismv+1)*apars[4]  + apars[5])/denom -1; /* the final -1 takes us back to indexed-from-0 coords. */
    denom = denom * denom;
    /* H = [h11,h12;h21,h22] is the Hessian of the transform given by apars, evaluated at
       the point (ismv,ismh) maps to under that transform. */
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

void Wgauss(mxArray **W, double *La, double *Lb, const double *M, const int biv, const int bih, const int nzmax_in, const int K, const int *smvs, const int* smhs, const int KM) {
    int ypt = 0; /* this will be a pointer int the array Y holding the relevant low-res pixel values. */
    int hipt = 0, lowh, highh, lowv, highv, nzmax = nzmax_in;
    double *h = (double *) mxCalloc((int)(biv*bih),sizeof(double)); /* This is enough space to store a nonsparse row of W as I compute it. */
    double myscaler, b11, b12, b22, nh, nv, nuh, nuv, deltav, deltah, tval, thresh;
    const double *Hph;
    uint kj, ismv, ismh, ibiv, ibih;

    W[0] = mxCreateSparse(biv*bih,KM,nzmax,mxREAL);
    double *Wdata = (double*) mxGetPr(W[0]);
    int *Wi = (int*) mxGetIr(W[0]);
    int *Wj = (int*) mxGetJc(W[0]);
    
    for (kj=0; kj<K; kj++) { /* iterate over the K low-res images */
        thresh = stdlim*stdlim*M[kj*11+10];
        for (ismh=0; ismh < smhs[kj]; ismh++) { /* iterate over low-res image's rows. */
            for (ismv=0; ismv < smvs[kj]; ismv++) { /* scan down each column in turn. */
                La[ypt] = M[kj*11+8]; /* multiplicative photometric parameter. */
                Lb[ypt] = M[kj*11+9]; /* additive photometric parameter. */
                Wj[ypt] = hipt; /* data for this column starts at index hipt in Wdata and Wi. */
                ypt++;
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
                            tval = b11*nh*nh + b12*nh*nv + b22*nv*nv;
                            if (tval<thresh) {
                                h[ibih*biv+ibiv] = exp(-0.5*tval); /* This evaluates the Gaussian kernel value. */
                                myscaler = myscaler + h[ibih*biv+ibiv]; /* I'll need this value later for normalization. */
                            } else { h[ibih*biv+ibiv] = 0; }
                        }
                    }
                    if (myscaler > 0.0000001) { /* if all the values are low, I'll flag the pixel as unsupported rather than assume good values. */
                        for (ibih = lowh; ibih < highh; ibih++) { /* iterate back across the patch, storing the normalized values in W. */
                            for (ibiv = lowv; ibiv < highv; ibiv++) {
                                if (h[ibih*biv+ibiv]>0) {
                                    if (hipt>=nzmax) { /* If there's not enough space in W, increase its size. */
                                        int nzmaxnew = (int) ceil(1+1.1*nzmax);
                                        printf("Warning! Having to increase size of output array from %i to %i!\n",nzmax,nzmaxnew);
                                        nzmax = nzmaxnew;
                                        mxSetNzmax(W[0], nzmax);
                                        mxSetPr(W[0], (double*)mxRealloc(Wdata, nzmax*sizeof(double)));
                                        mxSetIr(W[0], (int*)mxRealloc(Wi, nzmax*sizeof(int)));
                                        Wdata  = mxGetPr(W[0]);
                                        Wi = mxGetIr(W[0]);
                                    }
                                    Wi[hipt] = ibih*biv+ibiv; /* vertical coordinate in W matrix. */
                                    Wdata[hipt] = h[ibih*biv+ibiv]/myscaler; /* data element for W matrix (normalized). */
                                    hipt++;
                                }
                            }
                        }
                    }
                    else{printf("WARNING! Insufficient support for psf!\n");} /* This happens if "myscale" is too low. */
                }
            }
        }
    }
    Wj[KM] = hipt; /* This puts nnz into the last element of the Wj array. */
    mxFree(h); /* free my working space. */
    nzmax = hipt; /* Now free up any extra memory that was allocated for W but unused. */
    mxSetNzmax(W[0], nzmax);
    mxSetPr(W[0], (double*)mxRealloc(Wdata, nzmax*sizeof(double)));
    mxSetIr(W[0], (int*)mxRealloc(Wi, nzmax*sizeof(int)));
    return; /* All done! */
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    int nzmax = 0, K = 0, ypt = 0, mpt = 0, KM = 0;
    double *im, *H;
    const char *field_name;
    const mxArray *fa;
    uint i, kj;
    
   /* Check for proper number of arguments. */
    if(nrhs<3) { printf("Error: Too few inputs provided. I should see biv,bih,o\n"); return; }
    else if(nlhs>5) { printf("Error: Too many output arguments.\n"); return; }
    int biv = (int) mxGetPr(prhs[0])[0];
    int bih = (int) mxGetPr(prhs[1])[0];
    
   /* Now deal with input prhs[2], which is the K-element structure, o. */
    if (!mxIsStruct(prhs[2])) { printf("Input 3 must be a structure."); return; }
    K = mxGetNumberOfElements(prhs[2]);
    
    int *smvs = (int *) mxCalloc((int)(K),sizeof(int));
    int *smhs = (int *) mxCalloc((int)(K),sizeof(int));
    
    /* Check all the contents of the struct I need are actually there */
    const char *field_names[] = {"H", "la", "lb", "g", "im"}; /* Make sure "im" is last in list!  */
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
        smvs[kj] = (int) mxGetM(fa); /* "fa" still points to the im. */
        smhs[kj] = (int) mxGetN(fa);
        KM += smvs[kj] * smhs[kj];
    }
    
    /* Now create the output matrices. */
    plhs[1] = mxCreateDoubleMatrix(KM,1,mxREAL);
    double *Y = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(KM,1,mxREAL);
    double *La = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(KM,1,mxREAL);
    double *Lb = mxGetPr(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix(11,K,mxREAL);
    double *M = mxGetPr(plhs[4]);

    /* One more loop through the inputs to start reading stuff into outputs. */
    for (kj=0;kj<K;kj++) {
        fa = mxGetField(prhs[2], kj, "im");
        im = (double *) mxGetPr(fa);
        for (i=0; i<smvs[kj]*smhs[kj]; i++) {
            Y[ypt++] = im[i];
        }
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
        nzmax += ((int)ceil(3.14*16*M[mpt-1]*(M[mpt-7]*M[mpt-11]-M[mpt-8]*M[mpt-10])))*(smvs[kj]*smhs[kj]); /* this is pi*r^2*zoom*[number of pixels], where r=3*gamma. */
    }
    
    Wgauss(&plhs[0], La, Lb, M, biv, bih, nzmax, K, smvs, smhs, KM);
    
    return;
}
