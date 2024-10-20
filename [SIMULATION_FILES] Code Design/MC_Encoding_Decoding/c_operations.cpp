#include "mex.h"
#include <cmath>
#include <vector>

// Function to implement check_phi
double check_phi(double x) {
    if (x >= 33) {
        x = 33;
    }
    return std::log((std::exp(x) + 1) / (std::exp(x) - 1));
}

// MEX gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Input validation
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidNumInputs", "Two inputs required.");
    }
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 2) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:inputNotRealMatrix", "MP_MAT must be a real matrix.");
    }

    // Retrieve inputs
    double *MP_MAT = mxGetPr(prhs[0]);
    size_t numRows = mxGetM(prhs[0]);
    size_t numCols = mxGetN(prhs[0]);
    int M = static_cast<int>(mxGetScalar(prhs[1])); // M is the second input

    if (numCols < 4) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidMP_MAT", "MP_MAT must have at least 4 columns.");
    }

    // Output 'MP_MAT' (only the fourth column will be modified)
    plhs[0] = mxDuplicateArray(prhs[0]);  // Create a copy of the input MP_MAT for modification
    double *out_MP_MAT = mxGetPr(plhs[0]);

    // Temporary storage for each subMP_MAT (no need to create full array)
    std::vector<size_t> subRows;

    // Main loop over cDex
    for (int cDex = 1; cDex <= M; cDex++) {
        // Clear subRows and find rows where MP_MAT(:,1) == cDex
        subRows.clear();
        for (size_t i = 0; i < numRows; i++) {
            if (static_cast<int>(MP_MAT[i]) == cDex) {
                subRows.push_back(i);
            }
        }

        // Process each row in subMP_MAT for updating the fourth column only
        for (size_t subDex = 0; subDex < subRows.size(); subDex++) {
            double product_alpha = 1.0;
            double sum_beta_phi = 0.0;

            for (size_t i = 0; i < subRows.size(); i++) {
                if (i != subDex) {
                    size_t idx = subRows[i];
                    // Product of the second column (alpha, sign of the third column)
                    product_alpha *= std::signbit(MP_MAT[idx + 2 * numRows]) ? -1.0 : 1.0;

                    // Sum of check_phi of the third column (beta, absolute value)
                    sum_beta_phi += check_phi(std::abs(MP_MAT[idx + 2 * numRows]));
                }
            }

            // Set the fourth column for subDex, keeping the rest of MP_MAT unchanged
            out_MP_MAT[subRows[subDex] + 3 * numRows] = product_alpha * check_phi(sum_beta_phi);
        }
    }
}