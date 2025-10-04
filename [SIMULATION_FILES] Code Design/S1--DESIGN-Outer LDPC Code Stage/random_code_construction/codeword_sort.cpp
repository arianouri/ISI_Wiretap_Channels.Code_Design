#include "mex.h"
#include <vector>
#include <algorithm>

// MEX entry function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for proper number of input and output arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidNumInputs", "Two inputs required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:maxlhs", "Too many output arguments.");
    }

    // Get the inputs
    double *codeword_sys = mxGetPr(prhs[0]); // Input 1: codeword_sys
    double *permDex = mxGetPr(prhs[1]);      // Input 2: permDex

    // Get the length of the inputs
    mwSize N = mxGetM(prhs[0]); // Length of codeword_sys (assumed to be a column vector)
    if (N != mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:dimensionMismatch", "Inputs must have the same length.");
    }

    // Create output matrix
    plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    double *codeword = mxGetPr(plhs[0]);

    // Create a vector for permDex sorting
    std::vector<std::pair<double, int>> permDex_pairs(N);
    for (mwSize i = 0; i < N; i++) {
        permDex_pairs[i] = {permDex[i], i}; // Store value and index
    }

    // Sort permDex and retrieve true indices
    std::sort(permDex_pairs.begin(), permDex_pairs.end());

    // Populate the codeword output
    for (mwSize i = 0; i < N; i++) {
        codeword[i] = codeword_sys[permDex_pairs[i].second];
    }
}