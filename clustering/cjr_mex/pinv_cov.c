/////////////////////////////////////////////////////////////////////////
//
// pinv_cov.c
//
// Description:
// ------------
//	
//	Quickly compute the pseudo-inverse of a **COVARIANCE** matrix.
//
//	Matlab function interface looks like:
//		X = pinv_symm(A)
//	where A is the matrix to invert.
//
// This function will not work on matrices that are not covariance
// matrices, although the function *will work* for covariance matrices
// which are not positive definite.  Almost no error checking is
// performed, so if you call this function on a matrix that is not a
// covariance matrix, expect crashes or (worse!) incorrect answers.
//
// Performance and Robustness:
// ---------------------------------
//
// The code performance is compared to Matlab's PINV function.
//
// Figures are given in the following format: x (y) [x]. The x figures
// refer to this code compiled on Windows 2000 with Microsoft's Visual
// C/C++ 6.0 compiler. The (y) figures refer to this code compiled on
// Windows 2000 with the GCC 3.2 compiler (MinGQ distribution). The
// GCC compiler yields slightly better results than the Microsoft and
// Lcc compilers.  The [z] figures are for GCC version 3.3 on Linux.
//
// The first speed comparison (Speed Comparison 1) is based on the
// average call time for this function vs. PINV over matrices ranging
// from (2 x 2) to (300 x 300) in size. This gives an idea of how well
// the code performs and how it scales with problem size. A better
// comparison (Speed Comparison 2) is to actually run the code on a
// real problem.
//
// No comparison is made to the INV function, as experience suggests
// that INV is unsuitable for inverting covariance matrices.
//
// The code operates in one of two modes, depending whether the input
// matrix is positive definite or not. In the case where the matrix is
// not positive definite, the input matrix is projected into a
// Principal Components space, where the generalised inverse is
// computed.  The resulting matrix is then projected back into the
// natural space. There is some overhead is doing this, but as a
// result, a smaller matrix needs to be inverted. Therefore, it is not
// possible to state that computing the inverse of a non positive
// definite matris is the worse case. However, experience shows that
// positive definite matrices are inverted faster. If the matrix you
// are inverting is not positive definite, you can expect speed
// increases over Matlab's PINV for matrices up to XXX in size; above
// that, this code is less efficient than the PINV function.  (But you
// should not be working with non positive definite matrices anyway!)
//
// Speed Comparison 1:
//
// In the "best case", where the matrix is positive definite, you can
// expect an average speed increase of around ??? (???) [1.34] times
// over Matlab's PINV.
//
// In the "worse case", where the matrix is not positive definite, you
// can expect an average speed increase of around ??? (???) [1.05]
// times over Matlab's PINV.
//
// Since one generally works in Principal Component spaces, it is
// likely that the "best case" will be far more common than the "worst
// case".
//
// Speed Comparison 2:
//
// As an example, using pinv_cov to compute the conditional
// distribution of a Gaussian Mixture Model (a simple 2-D toy
// problem), using this code rather than PINV yields an average speed
// increase of ??? (???) [1.936] times. So, although the speedup
// suggested above seems modest, on a real problem using the pinv_cov
// code nearly doubles the speed of the conditioning experiment.
// 
// You can use the TIME_PINV.M function (in CJR's CVS repository)
// function to test the speed of this code compared to Matlab's
// PINV. However, the real test is in your research code!
//
// Robustness:
//
// As far as robustness is concerned, the code appears to be robust
// (compared to Matlab's PINV) provided your covariances are of an
// order greater than 1x10^-5. However you should, as a matter of
// course, ensure that your data are appropriately scaled
// (i.e. standardisation to Z-scores).
//
// You can use the TEST_PINV_COV_ROBUSTNESS.M function (in CJR's CVS
// repository) function to test the robustness of this code compared
// to Matlab's PINV. However, the real test is in your research code!
//
// Compilation Notes:
// ---------------------
//
// Be sure to use a decent C compiler. Matlab's Lcc compiler is pretty
// poor: using Microsoft's Visual C/C++ compiler resulted in code that
// runs in about 2/3 of the time that the Lcc compiler's code. Using
// GCC 3.2 on Windows (the MinGW distribution) yields slightly faster
// code (but it is slightly harder to get this compiler going on
// Windows).
//
// + Under Microsoft Visual C/C++ on Windows 2000:
//
// Do not use the -inline compiler flag: this resulted in slightly
// slower code.
//
// Just compile using "mex pinv_symm.c".
//
// I found that, on the WIndows platform using VC/C++, running on a
// Pentium 4, specifying the -G5 compiler flag (optimise for the
// pentium, rather than for the P4) resulted in better performance.
//
// + Under GCC 3.2 (MinGW distribution) on Windows 2000:
//
// See Matthew Bertt's page on compiling MEX files under Windows using
// GCC at http://www.mrc-cbu.cam.ac.uk/Imaging/gnumex20.html
//
//
// Known Bugs:
// --------------
//
// There are no known bugs for the two classes of matrix that the code
// is designed to handle.
//
// DO NOT try to invert matrices that are not covariance matrices, as
// you WILL get unwanted results!
//
// I use code from Numerical Recipes to compute the
// eigendecomposition. This does not appear to be as fast as Matlab's
// interface to LAPACK. I tried getting this code to use LAPACK, but
// it was such a pain in the ass that I gave up. We still get good
// speed improvements; the robustness does not appear to be affected,
// and we don't have to mess around with that hairy beast
// LAPACK. LAPACK may be faster and better than NR, but at least I can
// *use* NR! (If you are a MEX/LAPACK whizz, please consider the
// mortals among us before changing this code to use LAPACK -- it's no
// use if one can't build the code!)
//
// Revision Hsitory:
// -----------------
//
//	Written by: Chris Rose
//	
// Version 1.0	-	26 Oct 2003.
// Version 1.1	-	29 Oct 2003.
//     Function now copes with non positive definite matrices.
// Version 1.2	- 30 Oct 2003
//     Now compiled and tested under GCC 3.2 on Windows
// Version 1.3	-	3 Nov 2003
//     Now compiled and tested on Linux.
//
/////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////
//				HEADERS
/////////////////////////////////////////////////////////////////////////

// Headers to include
#include "mex.h"
#include <math.h>


/////////////////////////////////////////////////////////////////////////
//				DEBUGGING
/////////////////////////////////////////////////////////////////////////

// Define this to turn on debugging
#undef DEBUGGING

// Function to display a matrix on screen for debugging purposes
#ifdef DEBUGGING

// Macro to help us narrow down crashes to line numbers ("debugging for cheapskates")
#define DBGMSG mexPrintf("Executing: %i\n", __LINE__)

// Function to display a matrix in the Matlab window
void display_matrix(mxArray* aMatrix)
{
  int num_rows = mxGetM(aMatrix);
  int num_cols = mxGetN(aMatrix);
  double* data = mxGetPr(aMatrix);
  int row, col;
  
  for(row = 1; row <= num_rows; row++)
    {
      for(col = 1; col <= num_cols; col++)
	{
	  mexPrintf("%f\t", data[MATRIX_INDEX(num_rows, row, col)]);
	}
      mexPrintf("\n");
    }
  mexPrintf("\n\n");
}
#endif


/////////////////////////////////////////////////////////////////////////
//				MACROS
/////////////////////////////////////////////////////////////////////////

// Macro to compute an index into a 1D array for a 2D matrix m is
// number of row in matrix; i and j are the row and col to index.
// Macro returns the 1D index.
//
// A macro is more efficient than a function call (function call does
// not seem to get compiled out).
//
#define MATRIX_INDEX(m, i, j) ((m * (j -1)) + i -1)

static double sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double pythag(double a, double b)
{
  double absa, absb;
  absa = fabs(a);
  absb = fabs(b);
  if (absa > absb)
    return absa * sqrt(1.0 + SQR(absb / absa));
  else
    return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

// If eigenvalues are smaller than this, PCA will be performed on the
// input matrix to ensure it is positive-definite.
#define EIGENVALUES_TOLERANCE 0.0001

/////////////////////////////////////////////////////////////////////////
//				HOUSEHOLDER REDUCTION
/////////////////////////////////////////////////////////////////////////

void householder(double* aInputMatrix, int aN, double* aD, double* aE)
// Householder reduction of real symmetric matrix.
// Inputs:
//    aInputMatrix
//    -- the (aN x aN) matrix to transform.
// Outputs:
//    aInputMatrix
//    -- the (aN x aN) orthogonal matrix effecting the transformation
//    aD
//    -- an aN vector containing the diagonal elements of the
//	 tridiagonal matrix
//    aE
//    -- an aN vector of the of the off-diagonal elements with aE[0]
//    -- equal to zero.
{
  int l, k, j, i;
  double scale, hh, h, g, f;
  double temp;
  
  for(i = aN; i >= 2; i--)
    {
      l = i - 1;
      h = scale = 0.0;
      if (l > 1)
	{
	  for(k = 1; k <= l; k++)
	    { scale += fabs(aInputMatrix[MATRIX_INDEX(aN, i, k)]);}
	  if (scale == 0.0)
	    {aE[i-1] = aInputMatrix[MATRIX_INDEX(aN, i, l)];}
	  else
	    {
	      for(k = 1; k <= l; k++)
		{
		  aInputMatrix[MATRIX_INDEX(aN, i, k)] /= scale;
		  temp = aInputMatrix[MATRIX_INDEX(aN, i, k)];
		     // Use temp to avoid indexing the same element twice.
		  h += (temp * temp);
		}
	      f = aInputMatrix[MATRIX_INDEX(aN, i, l)];
	      g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
	      aE[i-1] = scale * g;
	      h -= f * g;
	      aInputMatrix[MATRIX_INDEX(aN, i, l)] = f - g;
	      f = 0.0;
	      for(j = 1; j <= l; j++)
		{
		  aInputMatrix[MATRIX_INDEX(aN, j, i)] = aInputMatrix[MATRIX_INDEX(aN, i, j)] / h;
		  g = 0.0;
		  for(k = 1; k <= j; k++)
		    {g += aInputMatrix[MATRIX_INDEX(aN, j, k)] * aInputMatrix[MATRIX_INDEX(aN, i, k)];}
		  for(k = j + 1; k <= l; k++)
		    {g += aInputMatrix[MATRIX_INDEX(aN, k, j)] * aInputMatrix[MATRIX_INDEX(aN, i, k)];}
		  aE[j-1] = g / h;
		  f += aE[j-1] * aInputMatrix[MATRIX_INDEX(aN, i, j)];
		}
	      hh = f / (h + h);
	      for(j = 1; j <= l; j++)
		{
		  f = aInputMatrix[MATRIX_INDEX(aN, i, j)];
		  aE[j-1] = g = aE[j-1] - hh * f;
		  for(k = 1; k <= j; k++)
		    {aInputMatrix[MATRIX_INDEX(aN, j, k)] -=
		       (f * aE[k-1] + g * aInputMatrix[MATRIX_INDEX(aN, i, k)]);}
		}
	    }
	}
      else
	{aE[i-1] = aInputMatrix[MATRIX_INDEX(aN, i, l)];}
      aD[i-1] = h;
    }
  aD[0] = 0.0;
  aE[0] = 0.0;
  for(i = 1; i <= aN; i++)
    {
      l = i - 1;
      if (aD[i-1])
	{
	  for(j = 1; j <= l; j++)
	    {
	      g = 0.0;
	      for(k = 1; k <= l; k++)
		{g += aInputMatrix[MATRIX_INDEX(aN, i, k)] * aInputMatrix[MATRIX_INDEX(aN, k, j)];}
	      for(k = 1; k <= l; k++)
		{aInputMatrix[MATRIX_INDEX(aN, k, j)] -= g * aInputMatrix[MATRIX_INDEX(aN, k, i)];}
	    }
	}
      aD[i-1] = aInputMatrix[MATRIX_INDEX(aN, i, i)];
      aInputMatrix[MATRIX_INDEX(aN, i, i)] = 1.0;
      for(j = 1; j <= l; j++)
	{aInputMatrix[MATRIX_INDEX(aN, j, i)] = aInputMatrix[MATRIX_INDEX(aN, i, j)] = 0.0;}
    }
}


/////////////////////////////////////////////////////////////////////////
//				EIGENVECTORS AND -VALUES
/////////////////////////////////////////////////////////////////////////

mxArray* do_eigs(double* aD, double* aE, int aN, double* aInputMatrix)
// QL algorithm with implicit shifts to determine eigenvectors and -values
// of a real symmetric matrix reduced to tridiagonal form using householder().
// Inputs:
//   aD
//   -- an aN vector of the diagonal elements of the tridiagonal matrix.
//   aE
//   --	an aN vector containing the subdiagonal elements of the tridiagonal
//	matrix, where aE[0] is arbitrary. aE is destroyed on output.
//   aN
//   -- the size of the input matrix.
//   aInputMatrix
//   -- the matrix returned by householder().
// Outputs:
//   aD
//   -- returns the aN eigenvalues.
//   aInputMatrix
//   -- the i-th column of aInputMatrix contains the normalised eigenvector
//      corresponding to aD[i].
{
  // A pointer to a matrix to hold the eigenvalues, and a pointer to its data
  mxArray* mxEigenValues = NULL;
  double* EigenValues = NULL;
  
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;
  
  for(i = 2; i <= aN; i++)
    {aE[i-2] = aE[i-1];}
  aE[aN-1] = 0.0;
  for(l = 1; l <= aN; l++)
    {
      iter = 0;
      do
	{
	  for(m = l; m <= aN - 1; m++)
	    {
	      dd = fabs(aD[m-1]) + fabs(aD[m]);
	      if ((double) (fabs(aE[m-1]) + dd) == dd) break;
	    }
	  if (m != 1)
	    {
	      if (iter++ == 40)
		{mexErrMsgTxt("Too many iterations in QL algorithm");}
	      g = (aD[l] - aD[l-1]) / (2.0 * aE[l-1]);
	      r = pythag(g, 1.0);
	      g = aD[m-1] - aD[l-1] + aE[l-1] / (g + SIGN(r, g));
	      s = c = 1.0;
	      p = 0.0;
	      for(i = m-1; i >= l; i--)
		{
		  f = s * aE[i-1];
		  b = c * aE[i-1];
		  aE[i] = (r = pythag(f, g));
		  if (r == 0.0)
		    {
		      aD[i] -= p;
		      aE[m-1] = 0.0;
		      break;
		    }
		  s = f / r;
		  c = g / r;
		  g = aD[i] - p;
		  r = (aD[i-1] - g) * s + 2.0 * c * b;
		  aD[i] = g + (p = s * r);
		  g = c * r - b;
		  for(k = 1; k <= aN; k++)
		    {
		      f = aInputMatrix[MATRIX_INDEX(aN, k, i+1)];
		      aInputMatrix[MATRIX_INDEX(aN, k, i+1)] =
			s * aInputMatrix[MATRIX_INDEX(aN, k, i)] + c * f;
		      aInputMatrix[MATRIX_INDEX(aN, k, i)] =
			c * aInputMatrix[MATRIX_INDEX(aN, k, i)] - s * f;
		    }
		}
	      if (r == 0.0 && i >= 1) continue;
	      aD[l-1] -= p;
	      aE[l-1] = g;
	      aE[m-1] = 0.0;
	    }
	} while (m != l);
    }
  
  // Create a matrix to hold the eigenvalues on the diagonal, and populate it.
  mxEigenValues = mxCreateDoubleMatrix(aN, aN, mxREAL);
  EigenValues = mxGetPr(mxEigenValues);
  for(i = 1; i <= aN; i++)
    {EigenValues[MATRIX_INDEX(aN, i, i)] = aD[i-1];}

  // Return the eigenvalues in this matrix, the eigenvectors are
  // written the the memory location of aInputMatrix
  return mxEigenValues;
}

mxArray* EigenVectorsAndValues(double* aCovMat, int aNumRows)
// Compute the eigenvectors and -values for a covariance matrix.
// aCovMat is the covariance matrix, aNumRows is the number of rows in
// the covariance matrix. The function returns a diagonal matrix
// containing the eigenvalues, and the covariance matrix is replaced
// with a matrix of eigenvectors, where each eigenvector is a column.
// The ordering of the eigenvectors and values is smallest first,
// which is opposite to what Matlab returns.
{
  mxArray* eigen_values;
  
  // create some working space
  double* D = mxMalloc(sizeof(double) * aNumRows);
  double* E = mxMalloc(sizeof(double) * aNumRows);
  
  // Compute the householder reduction of aCovMat
  householder(aCovMat, aNumRows, D, E);
  
  // Compute the eigendecomposition
  eigen_values = do_eigs(D, E, aNumRows, aCovMat);
  
  // free the working space
  mxFree(D);
  mxFree(E);
  
  // now return the eigenvalues (the -vectors are now stored in the
  // same matrix as the original covariance matrix was)
  return eigen_values;
}

/////////////////////////////////////////////////////////////////////////
//				MATRIX MULTIPLICATION
/////////////////////////////////////////////////////////////////////////

mxArray* MatrixMatrixMult(mxArray* aLeftMatrix, mxArray* aRightMatrix)
// Multiply two matrices and return a newmatrix containing the result
{
  double* result;
  int i, j, k, counter;
  
  // Find out the sizes of the two inputs, so we can create the matrix to return.
  int num_rows_left = mxGetM(aLeftMatrix);
  int num_cols_left = mxGetN(aLeftMatrix);
  int num_cols_right = mxGetN(aRightMatrix);
    // num_rows_right is equal to num_cols_left, by definition of matrix multiplication.
  double* Left = mxGetPr(aLeftMatrix);
  double* Right = mxGetPr(aRightMatrix);
  
  // Create the matrix to return
  mxArray* mxResult = mxCreateDoubleMatrix(num_rows_left, num_cols_left, mxREAL);
  result = mxGetPr(mxResult);
  
  // Do the multiplication.
  for(i = 1; i <= num_rows_left; i++)
    {
      for(j = 1; j <= num_cols_right; j++)
	{
	  counter = MATRIX_INDEX(num_rows_left, i, j);
	  for(k = 1; k <= num_cols_left; k++)
	    {
	      result[counter] +=
		Left[MATRIX_INDEX(num_rows_left, i, k)] * Right[MATRIX_INDEX(num_cols_left, k, j)];
	    }
	}
    }
  return mxResult;
}

mxArray* MatrixTransposeMult(mxArray* aLeftMatrix, mxArray* aRightMatrix)
// Compute the matrix multiplication of the left matrix with the transpose of the right matrix
{
  double* result;
  int i, j, k, counter;
  
  // Find out the sizes of the two inputs, so we can create the matrix to return.
  int num_rows_left = mxGetM(aLeftMatrix);
  int num_rows_right = mxGetM(aRightMatrix);
  int num_cols_right = mxGetN(aRightMatrix);
  double* Left = mxGetPr(aLeftMatrix);
  double* Right = mxGetPr(aRightMatrix);
  
  // Create the result matrix
  mxArray* mxResult = mxCreateDoubleMatrix(num_rows_left, num_rows_right, mxREAL);
  result = mxGetPr(mxResult);
  
  // Do the multiplication
  for(i = 1; i <= num_rows_left; i++)
    {
      for(j = 1; j <= num_rows_right; j++)
	{
	  counter = MATRIX_INDEX(num_rows_left, i, j);
	  for(k = 1; k <= num_cols_right; k++)
	    {
	      result[counter] +=
		Left[MATRIX_INDEX(num_rows_left, i, k)] * Right[MATRIX_INDEX(num_rows_right, j, k)];
	    }
	}
    }
  return mxResult;
}


/////////////////////////////////////////////////////////////////////////
//				ENTRY POINT
/////////////////////////////////////////////////////////////////////////

// Here's the MEX function
void mexFunction(int nlhs, mxArray *plhs[], // the left-hand side
		 int nrhs, const mxArray *prhs[]) // the right-hand side
{
  // integers
  int num_rows_A, orig_num_rows_A; // num rows and cols of A
  int r = 0;
  int col_counter, row_counter, k;
  int num_positive_eigenvalues;
  
  // This will hold the result of the inverse.
  double* ret_array = NULL;
  
  // these will point to a copy of the input matrix and its data
  mxArray* mxA = NULL;
  double* A = NULL;
  
  // create other array pointers and pointers to their data
  mxArray* mxS = NULL;
  double* S = NULL;
  mxArray* mxV = NULL;
  double* V = NULL;
  mxArray* mxP = NULL; // these are used if we need to do a PCA on the input matrix
  double* P = NULL;
  
  // doubles
  const double eps = mxGetEps();
  const double eigenvalues_tolerance = EIGENVALUES_TOLERANCE;
  double temporary_result;
  double max_S; // will hold the maximum value of the diagonal of S
  
		// Here's the algorithm:
  
		// Copy the input matrix
  mxA = mxDuplicateArray(prhs[0]);
  A = mxGetPr(mxA);
  
  // Get the number of rows in A; remember, A is square!
  num_rows_A = mxGetM(mxA);
  orig_num_rows_A = num_rows_A;
  num_positive_eigenvalues = num_rows_A; // assume all eigenvalues are positive
  
  // Compute the Eigen-decomposition of A.
  mxV = mxDuplicateArray(mxA);
    // V will hold the eigenvalues, but needs to hold the covariance
    // matrix to start with
  V = mxGetPr(mxV);
  mxS = EigenVectorsAndValues(V, num_rows_A);
  S = mxGetPr(mxS);
  // V now holds the eigenvectors
  // S now holds the eigenvalues
  
#ifdef DEBUGGING
  mexPrintf("A is\n");
  display_matrix(mxA);
  mexPrintf("S is\n");
  display_matrix(mxS);
  mexPrintf("V is\n");
  display_matrix(mxV);
#endif
  
  if (num_rows_A > 1)
    {
      //  get the maximum value of the diagonal of S, and note which
      // (if any) of the elements of the diagonal are zero
      max_S = S[num_rows_A - 1]; // we know this because this is how the eigendecomposition works
      for(k = 1; k <= num_rows_A; k++)
	{
	  if (S[MATRIX_INDEX(num_rows_A, k, k)] < eigenvalues_tolerance)
	    {
	      num_positive_eigenvalues--;
	    }
	}
    }
  else if (1 == num_rows_A)
    {
      if (A[0] < 0.0)
	{
	  // this algorithm does not work for negative scalars
	  mexErrMsgTxt("Error in pinv_symm (MEX Function): The scalar passed was negative.");
	}
      max_S = S[0];
    }
  else
    {
      // we shouldn't get here: m should be 1 or greater!
      mexErrMsgTxt("Error in pinv_symm (MEX Function): The number of rows in the input matrix seems to be less than 1.");
    }
  
  // Test to see if the matrix is positive definite; if not, project A
  // into the Principal Component space, then compute the
  // eigendecomposition of this and set any other variables that need
  // changing to make this reformation of the problem work with the
  // general code that follows.
  if (num_positive_eigenvalues < num_rows_A)
    {
      // 1. Compute the new A and compute the PCA transformation
      // 1. matrix, P, which is used at the end of the MEX function.
      
      // Delete and create a new A
      mxDestroyArray(mxA);
      mxA = NULL;
      A = NULL;
      mxA = mxCreateDoubleMatrix(num_positive_eigenvalues, num_positive_eigenvalues, mxREAL);
      A = mxGetPr(mxA);
      
      // Create P
      mxP = mxCreateDoubleMatrix(num_rows_A, num_positive_eigenvalues, mxREAL);
      P = mxGetPr(mxP);
      
      k = 1; // counts across the diagonal of the new A
      for(col_counter =1 ; col_counter<= num_rows_A; col_counter++)
	{
	  if (S[MATRIX_INDEX(num_rows_A, col_counter, col_counter)] > eigenvalues_tolerance)
	    {
	      // This column's corresponding eigenvalue is > 0, so set
	      // this diagonal element of A to the corresponding
	      // eigenvalue.
	      A[MATRIX_INDEX(num_positive_eigenvalues, k, k)] =
		S[MATRIX_INDEX(num_rows_A, col_counter, col_counter)];
	      
	      // This column's corresponding eigenvalue is > 0, so
	      // copy this column into P.
	      for(row_counter = 1; row_counter <= num_rows_A; row_counter++)
		{
		  P[MATRIX_INDEX(num_rows_A, row_counter, k)] =
		    V[MATRIX_INDEX(num_rows_A, row_counter, col_counter)];
		}
	      k++;
	    }
	}
      
      // 2. Compute S
      
      // The eigenvalues S are now the same as the diagonal of A
      // (i.e. S := A)
      mxDestroyArray(mxS);
      mxS = NULL;
      S = NULL;
      mxS = mxDuplicateArray(mxA);
      S = mxGetPr(mxS);
      
      //. Compute V
      
      // The eigenvector matrix is an identity matrix of the same size
      // as the new A
      mxDestroyArray(mxV);
      mxV = NULL;
      V = NULL;
      mxV = mxCreateDoubleMatrix(num_positive_eigenvalues, num_positive_eigenvalues, mxREAL);
      V = mxGetPr(mxV);
      for(k=1; k <= num_positive_eigenvalues; k++)
	{
	  V[MATRIX_INDEX(num_positive_eigenvalues, k, k)] = 1.0;
	}
      
      // Set any variables that need changing to work with the generic
      // pseudo-inverse code below.
      num_rows_A = num_positive_eigenvalues;
      
#ifdef DEBUGGING
      mexPrintf("A is\n");
      display_matrix(mxA);
      mexPrintf("S is\n");
      display_matrix(mxS);
      mexPrintf("V is\n");
      display_matrix(mxV);
      mexPrintf("P is\n");
      display_matrix(mxP);
#endif
    }
  
  temporary_result = num_rows_A * max_S * eps; // temporary_result holds the tolerance
  for(k = 1; k <= num_rows_A; k++)
    {
      if(S[MATRIX_INDEX(num_rows_A, k, k)] > temporary_result)
	{
	  r++;
	}
    }
  
  plhs[0] = mxCreateDoubleMatrix(num_rows_A, num_rows_A, mxREAL);
  
  if (0 != r)
    {
      // This will point to the memory location of the matrix element
      // currently being considered
      double* pointer_to_current_result_element = NULL;
      
      ret_array = mxGetPr(plhs[0]);
      for(row_counter = 1; row_counter <= num_rows_A; row_counter++)
	{
	  // The inverse of a symmetric matrix is symmetric, so only
	  // compute the upper triangle, and copy to the lower
	  // triangle to save some computations.
	  for(col_counter = row_counter; col_counter <= num_rows_A; col_counter++)
	    {
	      // Compute the upper triangle value.
	      pointer_to_current_result_element =
		&ret_array[MATRIX_INDEX(num_rows_A, row_counter, col_counter)];
	      for(k = 1; k <= r; k++)
		{
		  *pointer_to_current_result_element +=
		    (V[MATRIX_INDEX(num_rows_A, row_counter, k)] * 
		     V[MATRIX_INDEX(num_rows_A, col_counter, k)] / S[MATRIX_INDEX(num_rows_A, k, k)]);
		}
	      // This doesn't HAVE to be done on the DIAGONAL, but
	      // doing it anyway is very cheap and we avoid an if
	      // statement here.
	      ret_array[MATRIX_INDEX(num_rows_A, col_counter, row_counter)] =
		*pointer_to_current_result_element;
	    }
	}
    }
  
  // if we did a PCA on the input matrix, project our pseudo-inverse
  // back into the original space
  if (num_positive_eigenvalues < orig_num_rows_A)
    {
      // Let A# be our inverse. We want to compute P * A# * P^T
      mxArray* mxTemp1 = NULL;
      mxArray* mxTemp2 = NULL;
      
      // Compute P * A#
      mxTemp1 = MatrixMatrixMult(mxP, plhs[0]);
      // Compute mxTemp1 * P^T
      mxTemp2 = MatrixTransposeMult(mxTemp1, mxP);
      
      // Copy mxTemp3 into plhs[0]
      mxDestroyArray(plhs[0]);
      plhs[0] = NULL;
      plhs[0] = mxDuplicateArray(mxTemp2);
      
#ifdef DEBUGGING
      mexPrintf("mxTemp1 is\n");
      display_matrix(mxTemp1);
      mexPrintf("mxTemp2 is\n");
      display_matrix(mxTemp2);
#endif
      
      // Clean up the memory allocated manually for the PCA case
      mxDestroyArray(mxP);
      mxDestroyArray(mxTemp1);
      mxDestroyArray(mxTemp2);
      
    }
  
  // Clean up the memory allocated manually (apparently more efficient
  // than letting Matlab do it for us).
  mxDestroyArray(mxA);
  mxDestroyArray(mxV);
  mxDestroyArray(mxS);
  
  return;
}
