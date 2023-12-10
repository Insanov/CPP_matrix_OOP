#include "S21Matrix.hpp"
#include <limits>
#include <random>
#include <stdexcept>

using namespace std;

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ < 0 || cols_ < 0)
    throw std::length_error("Too short");

  matrix_ = new double *[rows_] {};
  for (int i = 0; i < rows_; i++) {
    try {
      matrix_[i] = new double[cols_]{};
    } catch (...) {
      for (int j = 0; j < i; j++)
        delete[] matrix_[j];
      delete[] matrix_;
      throw;
    }
  }
}

bool S21Matrix::is_equal_size(const S21Matrix &other) {
  return (cols_ == other.cols_) && (rows_ == other.rows_);
}

int S21Matrix::getRows() const { return rows_; }

int S21Matrix::getCols() const { return cols_; }

void S21Matrix::setRows(int rows) { rows_ = rows; }

void S21Matrix::setCols(int cols) { cols_ = cols; }

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.getRows()), cols_(other.getCols()) {
  matrix_ = new double *[rows_] {};
  for (int i = 0; i < rows_; i++) {
    try {
      matrix_[i] = new double[cols_]{};
    } catch (...) {
      for (int j = 0; j < i; j++)
        delete[] matrix_[j];
      delete[] matrix_;
      throw;
    }
  }

  copy(other.matrix_,
       other.matrix_ + rows_ * sizeof(double *) + cols_ * sizeof(double),
       matrix_);
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (this != &other) {
    delete matrix_;
    copy(other.matrix_,
         other.matrix_ + rows_ * sizeof(double *) + cols_ * sizeof(double),
         this->matrix_);
  }

  return *this;
}

double &S21Matrix::operator()(int row, int col) const {
  if (row < 0 || col < 0 || rows_ < row || cols_ < col)
    throw std::out_of_range("Zhopa");
  return matrix_[row][col];
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  S21Matrix tmp{*this};
  tmp.SumMatrix(other);
  return tmp;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  S21Matrix tmp{*this};
  tmp.SumMatrix(other);
  return tmp;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(double number) {
  this->MulNumber(number);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  S21Matrix tmp{*this};
  tmp.MulMatrix(other);
  return tmp;
}

S21Matrix S21Matrix::operator*(double number) const {
  S21Matrix tmp{*this};
  tmp.MulNumber(number);
  return tmp;
}

bool S21Matrix::operator==(const S21Matrix &other) {
  return this->EqMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.getRows()), cols_(other.getCols()), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.setRows(0);
  other.setCols(0);
}

S21Matrix::~S21Matrix() { delete[] matrix_; }

void S21Matrix::printMatrix() {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++)
      cout << matrix_[i][j] << " ";
    cout << endl;
  }
}

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  if (!is_equal_size(other))
    return false;

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++)
      if (abs(matrix_[i][j] - other.matrix_[i][j]) < 1e7)
        return false;

  return true;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (!is_equal_size(other))
    return;

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++)
      matrix_[i][j] += other.matrix_[i][j];
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (!is_equal_size(other))
    return;

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++)
      matrix_[i][j] -= other.matrix_[i][j];
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++)
      matrix_[i][j] *= num;
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.rows_)
    return;

  S21Matrix res(rows_, other.cols_);
  for (int i = 0; i < res.rows_; i++)
    for (int j = 0; j < res.cols_; j++)
      for (int k = 0; k < cols_; k++)
        res.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];

  this->cols_ = other.cols_;
  *this = std::move(res);
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix res(cols_, rows_);

  for (int i = 0; i < res.rows_; i++)
    for (int j = 0; j < res.cols_; j++)
      res(i, j) = (*this)(j, i);

  return res;
}

void S21Matrix::row_inner_copy(double *src, int row) {
  for (int i = 0; i < cols_; i++)
    (*this)(row, i) = src[i];
}

void S21Matrix::row_outer_copy(double *dest, int row) {
  for (int i = 0; i < cols_; i++)
    dest[i] = (*this)(row, i);
}
void S21Matrix::row_swap(int i, int j) {
  if (i != j)
    for (int k = 0; k < cols_; k++)
      swap((*this)(i, k), (*this)(j, k));
}

void mult_vector_by_scalar(int size, double *arr, double scalar) {
  for (int i = 0; i < size; i++)
    arr[i] *= scalar;
}

void rows_sub(double *minuend, double *subtrahend, int size) {
  for (int i = 0; i < size; i++)
    minuend[i] -= subtrahend[i];
}

int S21Matrix::column_reduction(S21Matrix &copied, double *final_divisor) {
  for (int i = 0; i < this->cols_; i++) {
    if (!copied(i, i)) {
      int flag = 0;
      for (int k = i; k < this->rows_; k++)
        if (copied(k, k)) {
          copied.row_swap(i, k);
          *final_divisor *= -1;
          flag = 1;
        }
      if (!flag) {
        *final_divisor = 0;
        return 0;
      }
    }
    for (int j = i + 1; j < this->rows_; j++) {
      double *arr_h = new double[2 * this->cols_];
      double *arr_l = new double[2 * this->cols_];
      copied.row_outer_copy(arr_h, i);
      copied.row_outer_copy(arr_l, j);
      mult_vector_by_scalar(this->cols_, arr_h, copied(j, i));
      mult_vector_by_scalar(this->cols_, arr_l, copied(i, i));
      *final_divisor *= copied(i, i);
      rows_sub(arr_l, arr_h, this->cols_);
      copied.row_inner_copy(arr_l, j);
      delete[] arr_l;
      delete[] arr_h;
    }
  }

  return 1;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_)
    throw out_of_range("Zhopa");

  double result = 1;
  double final_divisor = 1;

  S21Matrix copied(rows_, cols_);
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++)
      copied(i, j) = (*this)(i, j);

  column_reduction(copied, &final_divisor);
  if (!final_divisor)
    return 0;

  for (int i = 0; i < cols_; i++)
    result *= copied(i, i);
  result /= final_divisor;

  return result;
}

S21Matrix S21Matrix::CalcComplemets() {
  if (rows_ != cols_)
    throw out_of_range("Zhopa");

  S21Matrix result(rows_, cols_);
  if (rows_ == 1 && cols_ == 1) {
    result(0, 0) = 1;
  } else {
    S21Matrix helper(rows_ - 1, cols_ - 1);
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++) {
        create_helper(helper, i, j);
        double det = helper.Determinant();
        result(i, j) = det * pow(-1, (i + j));
      }
    }
  }

  return result;
}

void S21Matrix::create_helper(S21Matrix &to, int row, int col) {
  int help_row = 0;
  for (int i = 0; i < this->rows_; ++i) {
    if (i == row) continue;
    int help_col = 0;
    for (int j = 0; j < this->cols_; j++) {
      if (j == col) continue;
      to(help_row, help_col) = (*this)(i, j);
      help_col++;
    }
    help_row++;
  }
}

S21Matrix S21Matrix::InverseMatrix() {
  if (rows_ != cols_)
    throw out_of_range("Zhopa");

  if (cols_ == 1 && rows_ == 1) {
    S21Matrix result(1, 1);
    result(0, 0) = 1 / (*this)(0, 0);
    return result;
  } else {
    S21Matrix complements = CalcComplemets();
    S21Matrix result = Transpose();
    double det = Determinant();
    result.MulNumber(1.0 / det);
    return result;
  }
}

void matrixFill(S21Matrix *src) {
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> distribution(0, 100);
  for (int i = 0; i < src->getRows(); i++)
    for (int j = 0; j < src->getCols(); j++)
      (*src)(i, j) = distribution(gen);
}

int main() {
  S21Matrix a = S21Matrix(3, 3);
  matrixFill(&a);
  a.printMatrix();
  double res = a.Determinant();
  cout << res << endl;
  S21Matrix resM = a.InverseMatrix();
  resM.printMatrix();

  return 0;
}
