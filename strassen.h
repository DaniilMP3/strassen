#include <vector>
#include <complex>

typedef std::complex<long long> complex;
typedef std::vector<std::vector<complex>> matrix_complex;

void add(matrix_complex &A, matrix_complex &B, matrix_complex &result_matrix, int size);

void sub(matrix_complex &A, matrix_complex &B, matrix_complex &result_matrix, int size);

double nextpowerof2(int x);

void display(matrix_complex &matrix, int m, int n);

void strassen_util(matrix_complex &A, matrix_complex &B, matrix_complex &result_matrix, int size);

void strassen(matrix_complex &A, matrix_complex &B);