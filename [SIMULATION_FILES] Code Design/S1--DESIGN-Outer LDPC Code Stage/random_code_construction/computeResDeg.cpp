#include "mex.h"
#include <algorithm>
#include <vector>

// Helper function to find minimum value and index
std::pair<double, mwSize> findMin(const std::vector<double> &ResDeg) {
    auto minIt = std::min_element(ResDeg.begin(), ResDeg.end());
    return std::make_pair(*minIt, std::distance(ResDeg.begin(), minIt));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Input arguments: H, t, N, K, g, VecCol, VecRow
    double *H = mxGetPr(prhs[0]);
    mwSize t = (mwSize)mxGetScalar(prhs[1]);
    mwSize N = (mwSize)mxGetScalar(prhs[2]);
    mwSize K = (mwSize)mxGetScalar(prhs[3]);
    mwSize g = (mwSize)mxGetScalar(prhs[4]);

    double *VecCol = mxGetPr(prhs[5]);
    double *VecRow = mxGetPr(prhs[6]);

    mwSize rows = N - K - g;
    mwSize cols = N - t;

    // Allocate ResDeg array
    std::vector<double> ResDeg(cols, 0);

    // Compute ResDeg: Sum of rows of Hresedual
    for (mwSize j = 0; j < cols; ++j) {
        for (mwSize i = 0; i < rows; ++i) {
            ResDeg[j] += H[i + j * rows];  // Summing each column of Hresedual
        }
    }

    // Find minimum value and its index
    std::pair<double, mwSize> minPair = findMin(ResDeg);
    double minval = minPair.first;
    mwSize minind = minPair.second;

    // Handle case where minval == 0
    mwSize cnt = 0;
    std::vector<double> ResDegTmp = ResDeg;
    while (minval == 0 && cnt < cols) {
        cnt++;
        ResDegTmp[minind] = 999;  // Set to large number to skip
        minPair = findMin(ResDegTmp);
        minval = minPair.first;
        minind = minPair.second;
    }

    // Perform row-column swap based on algorithm (this would depend on specific logic)
    // Swap VecCol and VecRow based on the logic
    std::swap(VecCol[t], VecCol[minind + t]);
    std::swap(VecRow[t], VecRow[minind + t]);

    // Return the updated values back to MATLAB
    plhs[0] = mxCreateDoubleScalar(minval);  // Return minval
    plhs[1] = mxCreateDoubleScalar(minind);  // Return minind

    plhs[2] = mxCreateDoubleMatrix(1, cols, mxREAL);  // Return updated VecCol
    double *VecColOut = mxGetPr(plhs[2]);
    std::copy(VecCol, VecCol + cols, VecColOut);

    plhs[3] = mxCreateDoubleMatrix(1, rows, mxREAL);  // Return updated VecRow
    double *VecRowOut = mxGetPr(plhs[3]);
    std::copy(VecRow, VecRow + rows, VecRowOut);

    plhs[4] = mxCreateDoubleMatrix(rows, cols, mxREAL);  // Return updated H matrix
    double *HOut = mxGetPr(plhs[4]);
    std::copy(H, H + rows * cols, HOut);
}