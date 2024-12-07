#include "mex.h"
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>

// Function entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Input validation: check that there's exactly one input
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidNumInputs", "One input required.");
    }

    // Check if the input is a logical matrix (binary)
    if (!mxIsLogical(prhs[0])) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidInputType", "Input must be a binary (logical) matrix.");
    }

    // Get input matrix H
    bool *H = mxGetLogicals(prhs[0]);  // Using bool* for binary matrix
    mwSize M = mxGetM(prhs[0]);        // Rows
    mwSize N = mxGetN(prhs[0]);        // Columns

    // Matrix size check: N should be greater than M (full row rank assumption)
    if (N <= M) {
        mexErrMsgIdAndTxt("MATLAB:mexFunction:invalidMatrixSize", "The matrix should have more columns than rows.");
    }

    int K = N - M;  // Calculate K (assuming H is full row rank)

    // Prepare output: Copy input matrix H into output
    plhs[0] = mxDuplicateArray(prhs[0]); // Output matrix H_UpTriAng (copy of H)
    bool *H_UpTriAng = mxGetLogicals(plhs[0]);  // Get logical array of output

    // Vectors for tracking row and column operations
    std::vector<int> VecRow(M);
    std::vector<int> VecCol(N);
    for (int i = 0; i < M; ++i) {
        VecRow[i] = i + 1;  // 1-based indexing for MATLAB compatibility
    }
    for (int j = 0; j < N; ++j) {
        VecCol[j] = j + 1;  // 1-based indexing
    }

    // Initialize variables
    int t = 0;
    int g = 0;

    // Random seed initialization for random selections
    std::srand(static_cast<unsigned int>(std::time(0)));

    // Main loop (t < N - K - g)
    while (t < N - K - g) {
        int cnt = 0;
        std::vector<int> ResDeg(N - t, 0);  // Residual degrees for each column

        // Calculate residual degrees
        for (int j = t; j < N; ++j) {
            int Sum = 0;
            for (int i = t; i < N - K - g; ++i) {
                Sum += H[i * N + j];  // H is a binary matrix
            }
            ResDeg[j - t] = Sum;  // Store sum as residual degree
        }

        // Find the column with the minimum residual degree
        auto minResDegIter = std::min_element(ResDeg.begin(), ResDeg.end());
        int minval = *minResDegIter;
        int minind = static_cast<int>(std::distance(ResDeg.begin(), minResDegIter));

        if (minval == 1) {  // Handle min degree == 1 case
            int col = t + minind + 1;  // MATLAB-style 1-based index
            std::vector<int> row;

            // Find rows with H(r, col) == 1
            for (int r = t; r < N - K - g; ++r) {
                if (H[r * N + col - 1] == 1) {
                    row.push_back(r + 1);  // Add row index (1-based)
                }
            }

            // Swap columns in VecCol and H
            std::swap(VecCol[t], VecCol[col - 1]);
            for (int r = 0; r < M; ++r) {
                std::swap(H_UpTriAng[r * N + t], H_UpTriAng[r * N + col - 1]);
            }

            // Swap rows in VecRow and H
            std::swap(VecRow[t], VecRow[row[0] - 1]);
            for (int j = 0; j < N; ++j) {
                std::swap(H_UpTriAng[t * N + j], H_UpTriAng[(row[0] - 1) * N + j]);
            }

            ++t;
        } else {  // Handle min degree > 1 case
            int col = t + minind + 1;  // MATLAB-style 1-based index
            std::vector<int> row;

            // Find rows with H(r, col) == 1
            for (int r = t; r < N - K - g; ++r) {
                if (H[r * N + col - 1] == 1) {
                    row.push_back(r + 1);  // Add row index (1-based)
                }
            }

            // Swap columns in VecCol and H
            std::swap(VecCol[t], VecCol[col - 1]);
            for (int r = 0; r < M; ++r) {
                std::swap(H_UpTriAng[r * N + t], H_UpTriAng[r * N + col - 1]);
            }

            // Swap rows in VecRow and H
            std::swap(VecRow[t], VecRow[row[0] - 1]);
            for (int j = 0; j < N; ++j) {
                std::swap(H_UpTriAng[t * N + j], H_UpTriAng[(row[0] - 1) * N + j]);
            }

            // Move rows to the bottom of the matrix
            for (int k = minval; k >= 2; --k) {
                std::swap(H_UpTriAng[(row[k - 1] - 1) * N + t], H_UpTriAng[(M - g + k - minval - 1) * N + t]);
                std::swap(VecRow[row[k - 1] - 1], VecRow[M - g + k - minval - 1]);
            }

            ++t;
            g += minval - 1;
        }
    }
}