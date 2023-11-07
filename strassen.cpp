#include <iostream>
#include <algorithm>
#include "strassen.h"

void add(matrix_complex &A, matrix_complex &B, matrix_complex &result_matrix, int size){
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            result_matrix[i][j] = A[i][j] + B[i][j];
        }
    }
}

void sub(matrix_complex &A, matrix_complex &B, matrix_complex &result_matrix, int size){
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            result_matrix[i][j] = A[i][j] - B[i][j];
        }
    }
}

double nextpowerof2(int x){
    return std::pow(2, int(std::ceil(std::log2(x))));
}

void display(matrix_complex &matrix, int m, int n){
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            if (j != 0){
                std::cout << "\t";
            }
            std::cout << matrix[i][j];
        }
        std::cout << std::endl;
    }
}

void strassen_util(matrix_complex &A, matrix_complex &B, matrix_complex &result_matrix, int size){
    if(size == 1){
        result_matrix[0][0] = A[0][0] * B[0][0];
        return;
    }
    int new_size = size / 2;
    std::vector<complex> z(new_size);

    matrix_complex a11(new_size, z), a12(new_size, z), a21(new_size, z), a22(new_size, z),
            b11(new_size, z), b12(new_size, z), b21(new_size, z), b22(new_size, z),
            c11(new_size, z), c12(new_size, z), c21(new_size, z), c22(new_size, z),
            m1(new_size, z), m2(new_size, z), m3(new_size, z), m4(new_size, z),
            m5(new_size, z), m6(new_size, z), m7(new_size, z),
            aResult(new_size, z), bResult(new_size, z);

    for(int i = 0; i < new_size; i++){
        for(int j = 0; j < new_size; j++){
            // Divide matrix into 4 sub-matrices
            a11[i][j] = A[i][j];
            a12[i][j] = A[i][j + new_size];
            a21[i][j] = A[i + new_size][j];
            a22[i][j] = A[i + new_size][j + new_size];

            b11[i][j] = B[i][j];
            b12[i][j] = B[i][j + new_size];
            b21[i][j] = B[i + new_size][j];
            b22[i][j] = B[i + new_size][j + new_size];
        }
    }
    add(a11, a22, aResult, new_size);     // a11 + a22
    add(b11, b22, bResult, new_size);    // b11 + b22
    strassen_util(aResult, bResult, m1, new_size);
    // m1 = (a11+a22) * (b11+b22)

    add(a21, a22, aResult, new_size); // a21 + a22
    strassen_util(aResult, b11, m2, new_size);
    // m2 = (a21+a22) * (b11)

    sub(b12, b22, bResult, new_size);      // b12 - b22
    strassen_util(a11, bResult, m3, new_size);
    // m3 = (a11) * (b12 - b22)

    sub(b21, b11, bResult, new_size);       // b21 - b11
    strassen_util(a22, bResult, m4, new_size);
    // m4 = (a22) * (b21 - b11)

    add(a11, a12, aResult, new_size);      // a11 + a12
    strassen_util(aResult, b22, m5, new_size);
    // m5 = (a11+a12) * (b22)

    sub(a21, a11, aResult, new_size);      // a21 - a11
    add(b11, b12, bResult, new_size);
    // m11 + b12
    strassen_util(aResult, bResult, m6, new_size);
    // m6 = (a21-a11) * (b11+b12)

    sub(a12, a22, aResult, new_size);      // a12 - a22
    add(b21, b22, bResult, new_size);      // b21 + b22

    strassen_util(aResult, bResult, m7, new_size);   // m7 = (a12-a22) * (b21+b22)

    // Calculate sub-matrices c_ij

    // c11 = m1 + m4 - m5 + m7
    add(m1, m4, aResult, new_size);       // m1 + m4
    add(aResult, m7, bResult, new_size);  // m1 + m4 + m7
    sub(bResult, m5, c11, new_size);

    // c12 = m3 + m5
    add(m3, m5, c12, new_size);

    // c21 = m2 + m4
    add(m2, m4, c21, new_size);

    // c22 = m1 + m3 - m2 + m6
    add(m1, m3, aResult, new_size);       // m1 + m3
    add(aResult, m6, bResult, new_size);  // m1 + m3 + m6
    sub(bResult, m2, c22, new_size);

    for(int i = 0; i < new_size; i++){
        for(int j = 0; j < new_size; j++){
            result_matrix[i][j] = c11[i][j];
            result_matrix[i][j + new_size] = c12[i][j];
            result_matrix[i + new_size][j] = c21[i][j];
            result_matrix[i + new_size][j + new_size] = c22[i][j];
        }
    }
}

void strassen(matrix_complex &A, matrix_complex &B){
    int rows_a = A.size();
    int columns_a = A[0].size();
    int rows_b = B.size();
    int columns_b = B[0].size();

    int x = std::max({rows_a, rows_b, columns_a, columns_b});
    double s = nextpowerof2(x);

    std::vector<complex> z(s);
    matrix_complex Aa(s, z), Bb(s, z), Cc(s, z);

    for (unsigned int i = 0; i < rows_a; i++)
    {
        for (unsigned int j = 0; j < columns_a; j++)
        {
            Aa[i][j] = A[i][j];
        }
    }
    for (unsigned int i = 0; i < rows_b; i++)
    {
        for (unsigned int j = 0; j < columns_b; j++)
        {
            Bb[i][j] = B[i][j];
        }
    }
    strassen_util(Aa, Bb, Cc, s);
    std::vector<complex> temp1(columns_b);
    matrix_complex C(rows_a, temp1);
    for (unsigned int i = 0; i < rows_a; i++)
    {
        for (unsigned int j = 0; j < columns_b; j++)
        {
            C[i][j] = Cc[i][j];
        }
    }
    display(C, rows_a, columns_b);
}