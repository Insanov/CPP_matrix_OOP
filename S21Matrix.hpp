#ifndef S21MATRIX
#define S21MATRIX

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <random>
#include <cmath>

void mult_vector_by_scalar(int size, double *arr, double scalar);
void rows_sub(double *minuend, double *subtrahend, int size);

class S21Matrix {
  private:
    int rows_, cols_;
    double **matrix_;
    bool is_equal_size(const S21Matrix &other);
    void row_inner_copy(double *src, int row);
    void row_outer_copy(double *dest, int row);
    void create_helper(S21Matrix &to, int row, int col);

  public:
    S21Matrix();
    ~S21Matrix();
    S21Matrix(int rows, int cols);
    S21Matrix(const S21Matrix &other);
    S21Matrix(S21Matrix &&other);

    bool operator==(const S21Matrix &other);
    S21Matrix &operator=(const S21Matrix &other);
    double &operator()(int row, int col) const;
    S21Matrix operator+(const S21Matrix& other) const;
    S21Matrix& operator+=(const S21Matrix& other);
    S21Matrix operator-(const S21Matrix& other) const;
    S21Matrix& operator-=(const S21Matrix& other);
    S21Matrix& operator*=(double number);
    S21Matrix operator*(const S21Matrix& other) const;
    S21Matrix& operator*=(const S21Matrix& other);
    S21Matrix operator*(double number) const;

    void printMatrix();
    void row_swap(int i, int j);
    int column_reduction(S21Matrix &copied, double *final_divisor);

    int getRows() const;
    int getCols() const;
    void setRows(int rows);
    void setCols(int cols);

    bool EqMatrix(const S21Matrix &other);
    void SumMatrix(const S21Matrix &other);
    void SubMatrix(const S21Matrix &other);
    void MulNumber(const double num);
    void MulMatrix(const S21Matrix &other);
    S21Matrix Transpose();
    S21Matrix CalcComplemets();
    double Determinant();
    S21Matrix InverseMatrix();
};

#endif // S21MATRIX
