#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if (a.rows == b.rows && a.cols == b.cols) {
        Matrix result = create_matrix(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++) {
            for (int j = 0; j < a.cols; j++) {
                result.data[i][j] = a.data[i][j] + b.data[i][j];
            }
        }
        return result;
    }
    else{
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    if (a.rows == b.rows && a.cols == b.cols) {
        Matrix result = create_matrix(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++) {
            for (int j = 0; j < a.cols; j++) {
                result.data[i][j] = a.data[i][j] - b.data[i][j];
            }
        }
        return result;
    }
    else{
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    if (a.cols == b.rows) {
        Matrix result = create_matrix(a.rows, b.cols);
        for (int i = 0; i < a.rows; i++) {
            for (int j = 0; j < b.cols; j++) {
                result.data[i][j] = 0;
                for (int k = 0; k < a.cols; k++){
                    result.data[i][j] += a.data[i][k] * b.data[k][j];
                }
            }
        }
        return result;
    }
    else{
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
}

Matrix scale_matrix(Matrix a, double k)
{
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = k * a.data[i][j];
        }
    }
    return result;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix result = create_matrix(a.cols, a.rows);
    for (int i = 0; i < a.cols; i++) {
        for (int j = 0; j < a.rows; j++) {
            result.data[i][j] = a.data[j][i];
        }
    }
    return result;
}

double det_matrix(Matrix a)
{
    double result;
    if (a.rows == a.cols) {
        if (a.rows == 1) {
            result = a.data[0][0]; 
        }
        else{
            result = 0;
            int flag = 1;
            for (int col = 0; col < a.cols; col++) {
                Matrix sub;
                sub.rows = a.rows - 1;
                sub.cols = a.cols - 1;
                for (int i = 1; i < a.rows; i++) {
                    for (int j = 0; j < a.cols; j++) {
                        if (j < col) {
                            sub.data[i-1][j] = a.data[i][j];
                        } else if (j > col) {
                            sub.data[i-1][j-1] = a.data[i][j];
                        }
                    }
                }
                result += flag * a.data[0][col] * det_matrix(sub);
                flag = -flag;
            }
        }
        return result;
    }
    else {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
}

Matrix inv_matrix(Matrix a)
{
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    else if (det_matrix(a) == 0){
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    else {
        double det_a = det_matrix(a);
        Matrix cofactor = create_matrix(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++) {
            for (int j = 0; j < a.cols; j++) {
                Matrix sub = create_matrix(a.rows - 1, a.cols - 1);
                int sub_row = 0;
                for (int si = 0; si < a.rows; si++) {
                    if (si == i) continue;
                    int sub_col = 0;
                    for (int sj = 0; sj < a.cols; sj++) {
                        if (sj == j) continue;
                        sub.data[sub_row][sub_col] = a.data[si][sj];
                        sub_col++;
                    }
                    sub_row++;
                }
                cofactor.data[i][j] = pow(-1, i + j) * det_matrix(sub);
            }
        }
        Matrix adjugate = transpose_matrix(cofactor);
        return scale_matrix(adjugate, 1.0 / det_a);
    }
}

int rank_matrix(Matrix a)
{
    int rank = 0;
    int rows = a.rows, cols = a.cols;
    int min_dim = rows < cols ? rows : cols;
    Matrix tmp = create_matrix(rows, cols);
    for(int i=0; i<rows; i++)
        for(int j=0; j<cols; j++)
            tmp.data[i][j] = a.data[i][j];
    for(int col=0, row=0; col<cols && row<rows; col++) {
        int pivot = row;
        for(int i=row; i<rows; i++)
            if(fabs(tmp.data[i][col]) > fabs(tmp.data[pivot][col]))
                pivot = i;
        if(tmp.data[pivot][col] == 0) 
            continue;
        if(pivot != row) {
            for(int j=col; j<cols; j++) {
                double temp = tmp.data[row][j];
                tmp.data[row][j] = tmp.data[pivot][j];
                tmp.data[pivot][j] = temp;
            }
        }
        for(int i=row+1; i<rows; i++) {
            double factor = tmp.data[i][col] / tmp.data[row][col];
            for(int j=col; j<cols; j++)
                tmp.data[i][j] -= factor * tmp.data[row][j];
        }
        row++;
    }
    for(int i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            if(tmp.data[i][j] != 0) {
                rank++;
                break;
            }
        }
    }
    return rank;
}

double trace_matrix(Matrix a)
{
    double trace = 0;
    if (a.cols == a.rows) {
        for (int i = 0; i < a.cols; i++) {
            trace += a.data[i][i];
        }
        return trace;
    }
    else {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}