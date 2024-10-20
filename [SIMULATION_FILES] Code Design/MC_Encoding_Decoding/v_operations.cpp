#include "mex.h"
#include <vector>
#include <algorithm>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Input validation
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidNumInputs", "Two inputs required.");
    }
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotRealMatrix", "MP_MAT must be a real matrix.");
    }
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) != 2 || mxGetM(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotRealVector", "init_LLR must be a real vector.");
    }
    
    // Retrieve inputs
    double *MP_MAT = mxGetPr(prhs[0]);
    size_t numRows = mxGetM(prhs[0]);
    size_t numCols = mxGetN(prhs[0]);
    double *init_LLR = mxGetPr(prhs[1]);
    size_t N = mxGetN(prhs[1]);  // Length of init_LLR (as it's a row vector)

    if (numCols < 4) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidMP_MAT", "MP_MAT must have at least 4 columns.");
    }

    // Output 'MP_MAT' (modified)
    plhs[0] = mxDuplicateArray(prhs[0]);  // Create a copy of the input MP_MAT for modification
    double *out_MP_MAT = mxGetPr(plhs[0]);

    // Output 'tx_est' (vector)
    plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);
    double *tx_est = mxGetPr(plhs[1]);

    // Temporary storage for each subMP_MAT (no need to create full array)
    std::vector<size_t> subRows;

    // Main loop
    for (size_t vDex = 0; vDex < N; vDex++) {
        // Clear subRows and find rows where MP_MAT(:,2) == vDex+1
        subRows.clear();
        for (size_t i = 0; i < numRows; i++) {
            if (static_cast<size_t>(out_MP_MAT[i + numRows]) == vDex + 1) {
                subRows.push_back(i);
            }
        }

        // For each row in subMP_MAT, update the third column
        for (size_t subDex = 0; subDex < subRows.size(); subDex++) {
            size_t rowIdx = subRows[subDex];

            // Update the third column with init_LLR(vDex) + sum of other rows' 4th column values
            double sumVal = init_LLR[vDex];
            for (size_t i = 0; i < subRows.size(); i++) {
                if (i != subDex) {
                    sumVal += out_MP_MAT[subRows[i] + 3 * numRows];
                }
            }
            out_MP_MAT[rowIdx + 2 * numRows] = sumVal;  // Update the 3rd column
        }

        // Update tx_est(vDex) with the sum of the third column
        double sumEst = 0.0;
        for (size_t subDex = 0; subDex < subRows.size(); subDex++) {
            sumEst += out_MP_MAT[subRows[subDex] + 2 * numRows];  // 3rd column
        }
        tx_est[vDex] = sumEst;
    }
}