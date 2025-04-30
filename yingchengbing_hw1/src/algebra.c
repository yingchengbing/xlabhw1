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

Matrix add_matrix(Matrix a, Matrix b)//矩阵加法
{
    int i, j;
    if (a.rows == b.rows && a.cols == b.cols)
    {
        Matrix result = create_matrix(a.rows, a.cols);
        for (i = 0; i < a.rows; i++)
        {
            for (j = 0; j < a.cols; j++)
            {
                result.data[i][j] = a.data[i][j] + b.data[i][j];
            }
        }
        return result;
    }
    else
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
    }
}

Matrix sub_matrix(Matrix a, Matrix b)//矩阵减法
{
    int i, j;
    if (a.rows == b.rows && a.cols == b.cols)
    {
        Matrix result = create_matrix(a.rows, a.cols);
        for (i = 0; i < a.rows; i++)
        {
            for (j = 0; j < a.cols; j++)
            {
                result.data[i][j] = a.data[i][j] - b.data[i][j];
            }
        }
        return result;
    }
    else
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
    }
}

Matrix mul_matrix(Matrix a, Matrix b)//矩阵乘法
{
    int i, j, k;
    if (a.cols == b.rows)
    {
        Matrix result = create_matrix(a.rows, b.cols);
        for (i = 0; i < a.rows; i++)
        {
            for (j = 0; j < b.cols; j++)
            {
                result.data[i][j] = 0;
                for (k = 0; k < a.cols; k++)
                {
                    result.data[i][j] += a.data[i][k] * b.data[k][j];
                }
            }
        }
        return result;
    }
    else
    {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
    }
    
}

Matrix scale_matrix(Matrix a, double k)//矩阵数乘
{
    int i, j;
        Matrix result = create_matrix(a.rows, a.cols);
        for (i = 0; i < a.rows; i++)
        {
            for (j = 0; j < a.cols; j++)
            {
                result.data[i][j] = k*a.data[i][j];
            }
        }
        return result;
}

Matrix transpose_matrix(Matrix a)// 矩阵转置
{
    int i, j;
    Matrix result = create_matrix(a.cols, a.rows); // 转置后行列互换
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            result.data[j][i] = a.data[i][j]; // 行列互换
        }
    }
    return result;
}

double det_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }

    int n = a.rows;

    // 基础情况：1x1矩阵
    if (n == 1)
    {
        return a.data[0][0];
    }

    // 基础情况：2x2矩阵
    if (n == 2)
    {
        return a.data[0][0] * a.data[1][1] - a.data[0][1] * a.data[1][0];
    }

    // 递归计算行列式
    double determinant = 0;
    int col,i,j;
    // 对第一行的每个元素进行展开
    for (col = 0; col < n; col++)
    {
        // 创建子矩阵
        Matrix submatrix = create_matrix(n - 1, n - 1);
        for (i = 1; i < n; i++) // 从第1行开始
        {
            int sub_col = 0;
            for (j = 0; j < n; j++)
            {
                if (j == col) continue; // 跳过当前列
                submatrix.data[i - 1][sub_col] = a.data[i][j];
                sub_col++;
            }
        }

        // 递归计算子矩阵的行列式
        double sub_det = det_matrix(submatrix);

        // 加减交替
        determinant += (col % 2 == 0 ? 1 : -1) * a.data[0][col] * sub_det;
    }

    return determinant;

}

Matrix inv_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }

    int n = a.rows;
    Matrix result = create_matrix(n, n);
    Matrix augmented = create_matrix(n, 2 * n);

    // 初始化增广矩阵 [a | I]
    int i,j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            augmented.data[i][j] = a.data[i][j];
            augmented.data[i][j + n] = (i == j) ? 1 : 0; // 单位矩阵
        }
    }

    // 高斯消元
    int k;
    for (i = 0; i < n; i++)
    {
        // 寻找主元
        if (augmented.data[i][i] == 0)
        {
            int swap_row = -1;
            for (k = i + 1; k < n; k++)
            {
                if (augmented.data[k][i] != 0)
                {
                    swap_row = k;
                    break;
                }
            }
            if (swap_row == -1)
            {
                printf("Error: The matrix is singular.\n");
                return create_matrix(0, 0);
            }
            // 交换行
            for (j = 0; j < 2 * n; j++)
            {
                double temp = augmented.data[i][j];
                augmented.data[i][j] = augmented.data[swap_row][j];
                augmented.data[swap_row][j] = temp;
            }
        }

        // 将主元归一化
        double pivot = augmented.data[i][i];
        for (j = 0; j < 2 * n; j++)
        {
            augmented.data[i][j] /= pivot;
        }

        // 消去其他行的当前列
        for (k = 0; k < n; k++)
        {
            if (k != i)
            {
                double factor = augmented.data[k][i];
                for (j = 0; j < 2 * n; j++)
                {
                    augmented.data[k][j] -= factor * augmented.data[i][j];
                }
            }
        }
    }

    // 提取逆矩阵
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            result.data[i][j] = augmented.data[i][j + n];
        }
    }

    return result;
}

int rank_matrix(Matrix a)
{
    int rank = 0;
    int rows = a.rows;
    int cols = a.cols;
    double temp[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
    int i, j;

    // 复制矩阵数据到临时矩阵
    for (i= 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            temp[i][j] = a.data[i][j];
        }
    }

    // 高斯消元
    int col, row;
    for (col = 0, row = 0; col < cols && row < rows; col++)
    {
        // 寻找当前列的主元
        int pivot = row;
        for (i = row + 1; i < rows; i++)
        {
            if (fabs(temp[i][col]) > fabs(temp[pivot][col]))
            {
                pivot = i;
            }
        }

        // 如果主元为0，跳过当前列
        if (fabs(temp[pivot][col]) < 1e-10)
        {
            continue;
        }

        // 交换当前行和主元行
        for (j = 0; j < cols; j++)
        {
            double t = temp[row][j];
            temp[row][j] = temp[pivot][j];
            temp[pivot][j] = t;
        }

        // 将主元归一化
        double pivot_value = temp[row][col];
        for (j = 0; j < cols; j++)
        {
            temp[row][j] /= pivot_value;
        }

        // 消去其他行的当前列
        for (i = 0; i < rows; i++)
        {
            if (i != row)
            {
                double factor = temp[i][col];
                for (j = 0; j < cols; j++)
                {
                    temp[i][j] -= factor * temp[row][j];
                }
            }
        }

        // 处理下一行
        row++;
    }

    // 统计非零行的数量
    for (i = 0; i < rows; i++)
    {
        int non_zero = 0;
        for (j = 0; j < cols; j++)
        {
            if (fabs(temp[i][j]) > 1e-10)
            {
                non_zero = 1;
                break;
            }
        }
        rank += non_zero;
    }

    return rank;
}

double trace_matrix(Matrix a)
{
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }

    double trace = 0;
    int i;
    for (i = 0; i < a.rows; i++)
    {
        trace += a.data[i][i];
    }

    return trace;
}

void print_matrix(Matrix a)
{
    int i, j;
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}