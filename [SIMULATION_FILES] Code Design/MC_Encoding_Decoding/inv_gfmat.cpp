#include "mex.h"
#include <vector>
#include <algorithm>

// Function to calculate the inverse of a binary matrix
void inv_binary_matrix(std::vector<std::vector<bool>> &A, std::vector<std::vector<bool>> &inv_A, int row_A, int col_A) {
    int rank_A = row_A;  // Assuming full rank, since gfrank(A,2) functionality is not implemented here
    std::vector<std::vector<bool>> Aug_mtx(row_A, std::vector<bool>(2 * col_A, 0));
    
    // Initialize Augmented Matrix with Identity matrix
    for (int i = 0; i < row_A; i++) {
        for (int j = 0; j < col_A; j++) {
            Aug_mtx[i][j] = A[i][j];
        }
        Aug_mtx[i][col_A + i] = 1;  // Identity matrix part
    }
    
    if (rank_A == row_A) {
        // Gaussian elimination to convert to row echelon form
        for (int i = 0; i < col_A; i++) {
            // If diagonal element is 0, swap with a row below
            if (Aug_mtx[i][i] == 0) {
                for (int j = i + 1; j < row_A; j++) {
                    if (Aug_mtx[j][i] == 1) {
                        std::swap(Aug_mtx[i], Aug_mtx[j]);
                        break;
                    }
                }
            }
            std::vector<bool> temp_vec = Aug_mtx[i];
            // Eliminate below the diagonal
            for (int j = i + 1; j < row_A; j++) {
                if (Aug_mtx[j][i] == 1) {
                    for (int k = 0; k < 2 * col_A; k++) {
                        Aug_mtx[j][k] = Aug_mtx[j][k] ^ temp_vec[k];  // XOR for binary addition (mod 2)
                    }
                }
            }
        }

        // Back substitution to convert to reduced row echelon form
        for (int i = col_A - 1; i >= 0; i--) {
            std::vector<bool> temp_vec = Aug_mtx[i];
            for (int j = 0; j < i; j++) {
                if (Aug_mtx[j][i] == 1) {
                    for (int k = 0; k < 2 * col_A; k++) {
                        Aug_mtx[j][k] = Aug_mtx[j][k] ^ temp_vec[k];  // XOR for binary addition (mod 2)
                    }
                }
            }
        }

        // Extract the inverse from the augmented matrix
        for (int i = 0; i < row_A; i++) {
            for (int j = 0; j < col_A; j++) {
                inv_A[i][j] = Aug_mtx[i][j + col_A];
            }
        }
    } else {
        mexErrMsgIdAndTxt("inv_binary_matrix:invalidMatrix", "The matrix does not have an inverse!");
    }
}

// The entry point for the MEX function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("inv_binary_matrix:invalidNumInputs", "One input required.");
    }

    // Check if the input is a logical matrix (binary matrix)
    if (!mxIsLogical(prhs[0])) {
        mexErrMsgIdAndTxt("inv_binary_matrix:inputNotLogical", "Input matrix must be a binary logical matrix.");
    }

    // Get the input matrix dimensions
    int row_A = mxGetM(prhs[0]);
    int col_A = mxGetN(prhs[0]);

    // Create input matrix (logical, binary)
    std::vector<std::vector<bool>> A(row_A, std::vector<bool>(col_A, 0));
    mxLogical *input_ptr = mxGetLogicals(prhs[0]);

    for (int i = 0; i < row_A; i++) {
        for (int j = 0; j < col_A; j++) {
            A[i][j] = input_ptr[i + j * row_A];  // Logical (binary) values from input
        }
    }

    // Prepare the output matrix
    plhs[0] = mxCreateLogicalMatrix(row_A, col_A);  // Logical matrix for output (binary)
    std::vector<std::vector<bool>> inv_A(row_A, std::vector<bool>(col_A, 0));

    // Call the matrix inversion function
    inv_binary_matrix(A, inv_A, row_A, col_A);

    // Copy the result back to MATLAB (logical, binary matrix)
    mxLogical *output_ptr = mxGetLogicals(plhs[0]);
    for (int i = 0; i < row_A; i++) {
        for (int j = 0; j < col_A; j++) {
            output_ptr[i + j * row_A] = inv_A[i][j];
        }
    }
}