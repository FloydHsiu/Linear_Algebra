#include "matrix.h"


Matrix::Matrix()
{
	this->RowSize = 0;
	this->ColumnSize = 0;
}

Matrix::Matrix(int m, int n, std::vector<std::vector<double> > mtx_data)
{
	this->RowSize = m;
	this->ColumnSize = n;
	this->data.clear();
	this->data.assign(mtx_data.begin(), mtx_data.end());
}

Matrix::~Matrix() {

}

int Matrix::get_RowIndex() 
{
	return this->RowSize;
}

int Matrix::get_ColumnIndex() 
{
	return this->ColumnSize;
}

double Matrix::get_data(int m, int n) throw (MatrixException)
{
	if (m >= this->RowSize || n >= this->ColumnSize) {
		throw MatrixException("Matrix::get_data", 'A');
	}
	return this->data[m][n];
}

void Matrix::setData(int m, int n, std::vector<std::vector<double> > mtx_data)
{
	this->RowSize = m;
	this->ColumnSize = n;
	this->data.clear();
	this->data.assign(mtx_data.begin(), mtx_data.end());
}

void Matrix::show()
{
	for (int i = 0; i < this->get_RowIndex(); i++) {
		for (int j = 0; j < this->get_ColumnIndex(); j++) {
			std::cout << this->data[i][j] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << "," << std::endl;
}

//essential

Matrix Matrix::transpose(Matrix mtx)
{
	int _m = mtx.get_RowIndex();
	int _n = mtx.get_ColumnIndex();
	std::vector<std::vector<double> > temp(_n, std::vector<double>(_m));
	for (int i = 0; i < _m; i++) {
		for (int j = 0; j < _n; j++) {
			temp[j][i] = mtx.get_data(i, j);
		}
	}
	return Matrix(_n, _m, temp);
}

Matrix Matrix::reduceRowColumn(int m, int n) throw (MatrixException)
{
	if (!(m >= 0 && n >= 0 && m < this->RowSize && n < this->ColumnSize)) throw MatrixException("Matrix::reduceRowColumn", 'A');
	std::vector<std::vector<double> >tmp(this->RowSize - 1, std::vector<double>(this->ColumnSize - 1));
	int k = 0;
	for (int i = 0; i < this->RowSize - 1; i++) {
		int l = 0;
		if (m == k) k++;//當執行到要刪除列時，自動+1以跳過該列
		for (int j = 0; j < this->ColumnSize - 1; j++) {
			if (l == n) l++;//當執行到要刪除行時，自動+1以跳過該行
			tmp[i][j] = this->data[k][l ];
			l++;//往下一行
		}
		k++;//往下一列
	}
	return Matrix(this->RowSize-1, this->ColumnSize-1, tmp);
}

double Matrix::multiDiagonal()//對角線相乘
{
	double result = 1;
	for (int i = 0; i < this->RowSize; i++) {
		result *=this->data[i][i];
	}
	return result;
}

Matrix Matrix::operator+(Matrix &mtx) throw (MatrixException)
{
	if (this->RowSize != mtx.RowSize || this->ColumnSize != mtx.ColumnSize) {
		//check if m and n are the same
		throw MatrixException("Matrix::operator+", 'D');
	}
	std::vector<std::vector<double> > temp(mtx.RowSize, std::vector<double>(mtx.ColumnSize));
	for (int i = 0; i < mtx.RowSize; i++) {
		for (int j = 0; j < mtx.ColumnSize; j++) {
			temp[i][j] = mtx.data[i][j] + this->data[i][j];
		}
	}
	return Matrix(mtx.RowSize, mtx.ColumnSize, temp);
}

Matrix Matrix::operator-(Matrix & mtx) throw(MatrixException)
{
	if (this->RowSize != mtx.RowSize || this->ColumnSize != mtx.ColumnSize) {
		//check if m and n are the same
		throw MatrixException("Matrix::operator-", 'D');
	}
	std::vector<std::vector<double> > temp(mtx.RowSize, std::vector<double>(mtx.ColumnSize));
	for (int i = 0; i < mtx.RowSize; i++) {
		for (int j = 0; j < mtx.ColumnSize; j++) {
			temp[i][j] = this->data[i][j] - mtx.data[i][j];
		}
	}
	return Matrix(mtx.RowSize, mtx.ColumnSize, temp);
}

Matrix Matrix::operator*(Matrix &mtx)
{
	if (this->ColumnSize != mtx.RowSize ) {
		//check if m and n are the same
		throw MatrixException("Matrix::operator*", 'M');
	}
	std::vector<std::vector<double> > temp(this->RowSize, std::vector<double>(mtx.ColumnSize));//A[m*n] * B[n*o] = C[m*o]
	for (int i = 0; i < this->RowSize; i++) {
		for (int j = 0; j < mtx.ColumnSize; j++) {
			temp[i][j] = 0;
			for (int k = 0; k < this->ColumnSize; k++) {
				temp[i][j] += this->data[i][k] * mtx.data[k][j];
			}
		}
	}
	return Matrix(this->RowSize, mtx.ColumnSize, temp);
}

Matrix Matrix::operator*(double constant)
{
	std::vector<std::vector<double> > temp(this->RowSize, std::vector<double>(this->ColumnSize));
	for (int i = 0; i < this->RowSize; i++) {
		for (int j = 0; j < this->ColumnSize; j++) {
			temp[i][j] =  constant * this->data[i][j];
		}
	}
	return Matrix(this->RowSize, this->ColumnSize, temp);
}

//complex

std::vector<Matrix> Matrix::LU() throw (MatrixException)
{
	if (this->RowSize != this->ColumnSize) throw MatrixException("Matrix::LU", 'S');
	std::vector<std::vector<double> > L (this->RowSize, std::vector<double>(this->ColumnSize));
	std::vector<std::vector<double> > U (this->RowSize, std::vector<double>(this->ColumnSize));
	//Initialize U,L and let diagonal of L be 1
	for (int i = 0; i < this->RowSize; i++) {
		for (int j = 0; j < this->ColumnSize; j++) {
			if (i == j)L[i][j] = 1;
			else L[i][j] = 0;
			U[i][j] = 0;
		}
	}
	//
	if( this->data[0][0] != 0 ) U[0][0] = this->data[0][0] / L[0][0];
	else throw MatrixException("Matrix::LU", 'L');
	//First row of U, first column of L
	for (int j = 1; j < this->RowSize; j++) {
		U[0][j] = this->data[0][j] / L[0][0];
		L[j][0] = this->data[j][0] / U[0][0];
	}
	//
	for (int i = 1; i < this->RowSize - 1; i++) {
		//
		double tmp = 0;
		for (int k = 0; k < i; k++){
			tmp += L[i][k] * U[k][i];
		}
		double a_sub_tmp = this->data[i][i] - tmp;
		if (std::abs(a_sub_tmp) < 1e-8) throw MatrixException("Matrix::LU", 'L');
		else U[i][i] = a_sub_tmp / L[i][i];
		//
		for (int j = i + 1; j < this->RowSize; j++) {
			double tmp_u = 0;
			double tmp_l = 0;
			for (int k = 0; k < i; k++) {
				tmp_u += L[i][k] * U[k][j];
				tmp_l += L[j][k] * U[k][i];
			}
			U[i][j] = (this->data[i][j] - tmp_u) / L[i][i];
			L[j][i] = (this->data[j][i] - tmp_l) / U[i][i];
		}
	}
	//
	double tmp_nn = 0;
	for (int k = 0; k < this->RowSize-1; k++) {
		tmp_nn += L[this->RowSize-1][k] * U[k][this->RowSize-1];
	}
	U[this->RowSize-1][this->RowSize-1] = this->data[this->RowSize-1][this->RowSize-1] - tmp_nn;
	std::vector<Matrix> L_U;
	L_U.push_back(Matrix(this->RowSize, this->RowSize, L));
	L_U.push_back(Matrix(this->RowSize, this->RowSize, U));
	return L_U;
}

double Matrix::cofactor(int m, int n)
{
	Matrix tmp = this->reduceRowColumn(m,n);
	return pow(-1, (m+n))*(tmp.Det_Gauss());
}

Matrix Matrix::adjoint() throw(MatrixException)//先算矩陣的cofactor矩陣，再做轉置
{
	if (this->RowSize != this->ColumnSize) throw MatrixException("Matrix::adjoint", 'S');
	std::vector<std::vector<double> > tmp(this->RowSize, std::vector<double>(this->ColumnSize));
	for (int i = 0; i < this->RowSize; i++) {
		for (int j = 0; j < this->ColumnSize; j++) {
			tmp[i][j] = this->cofactor(i, j);
		}
	}
	return Matrix::transpose(Matrix(this->RowSize, this->ColumnSize, tmp));
}

Matrix Matrix::inverse() throw(MatrixException)
{
	if (this->Det_Gauss() == 0) throw MatrixException("Matrix::inverse", 'I');
	try {
		return this->adjoint()*(1 / this->Det_Gauss());
	}
	catch (MatrixException e) {
		throw MatrixException("Matrix::inverse", 'S');
	}
}

Matrix Matrix::gauss()
{
	std::vector<std::vector<double> >tmp;
	tmp.assign(this->data.begin(), this->data.end());
	int n = this->RowSize;
	//待修正，針對非n*n矩陣
	for (int i = 0; i < n; i++) {
		// Search for maximum in this column
		double maxEl = abs(tmp[i][i]);
		int maxRow = i;
		for (int k = i + 1; k<n; k++) {
			if (abs(tmp[k][i]) > maxEl) {
				maxEl = abs(tmp[k][i]);
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k < n; k++) {
			double temp = tmp[maxRow][k];
			tmp[maxRow][k] = tmp[i][k];
			tmp[i][k] = temp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k < n; k++) {
			double c = -tmp[k][i] / tmp[i][i];
			if (tmp[i][i] == 0) c = 1;//當此值為0時，因為有取絕對值最大值的關係，tmp[i][i]下面行全為0，所以不須乘以倍率
			for (int j = i; j < n ; j++) {
				if (i == j) {
					tmp[k][j] = 0;
				}
				else {
					tmp[k][j] += c * tmp[i][j];
				}
			}
		}
	}
	return Matrix(n, n ,tmp);
}

double Matrix::Det_Gauss()
{
	std::vector<std::vector<double> >tmp;
	tmp.assign(this->data.begin(), this->data.end());
	int n = this->RowSize;
	int swap_time = 0;
	for (int i = 0; i < n; i++) {
		// Search for maximum in this column
		double maxEl = abs(tmp[i][i]);
		int maxRow = i;
		for (int k = i + 1; k<n; k++) {
			if (abs(tmp[k][i]) > maxEl) {
				maxEl = abs(tmp[k][i]);
				maxRow = k;
			}
		}

		if (maxRow != i) swap_time++;//每做列交換一次，值就乘以-1

		// Swap maximum row with current row (column by column)
		for (int k = i; k < n; k++) {
			double temp = tmp[maxRow][k];
			tmp[maxRow][k] = tmp[i][k];
			tmp[i][k] = temp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k < n; k++) {
			double c = -tmp[k][i] / tmp[i][i];
			if (tmp[i][i] == 0) c = 1;//當此值為0時，因為有取絕對值最大值的關係，tmp[i][i]下面行全為0，所以不須乘以倍率
			for (int j = i; j < n; j++) {
				if (i == j) {
					tmp[k][j] = 0;
				}
				else {
					tmp[k][j] += c * tmp[i][j];
				}
			}
		}
	}
	return pow(-1, swap_time) * Matrix(n, n, tmp).multiDiagonal();
}

double Matrix::Det_RecursiveAndVector()
{
	double result  = Matrix::DetRecursive(*this);
	return result;
}

double Matrix::Det_RecursiveAndArray()
{
	int m = this->RowSize;
	int n = this->ColumnSize;
	return DetRecursive(m, n, StdVectorArrayTransform::parseToArray(m ,n , this->data));
}

int Matrix::rank()
{
	Matrix tmp = this->gauss();//先做高斯消去法
	int counter = 0;
	for (int i = 0; i < tmp.ColumnSize; i++) {
		if (std::abs(tmp.data[i][i]) < 1e-8) counter++;//計算pivot為零的個數
	}
	return this->ColumnSize - counter;
}

Matrix Matrix::solveLinearEquation(Matrix A, Matrix B) throw (MatrixException)// AX = B --> X = inverse(A)*B
{
	if (A.ColumnSize != A.RowSize || A.RowSize != B.RowSize) throw MatrixException("Matrix::solveLinearEquation", 'M');
	Matrix InverA = A.inverse();
	return InverA*B;
}

Matrix Matrix::LeastSquare(Matrix A, Matrix B) throw(MatrixException)//tanspose(A)AX = tanspose(A)B  --> X = inverse(tansepose(A)A)transpose(A)B 
{
	if (A.RowSize != B.RowSize) throw MatrixException("Matrix::LeastSquare", 'M');
	Matrix AtA = Matrix::transpose(A) * A;
	Matrix inverAtA = AtA.inverse();
	return inverAtA*Matrix::transpose(A)*B;
}

double Matrix::DetRecursive(Matrix &mtx) 
{
	int m = mtx.get_RowIndex(), n = mtx.get_ColumnIndex();
	double result = 0;
	if (m == 1 && n == 1) {
		//std::cout << mtx.get_data(0, 0) << std::endl;
		result = mtx.get_data(0, 0);
		return result;//簡化到matrix中只有一個element就return該值
	}
	for (int i = 0; i < m; i++) {
		//重複取新的matrix的第一列去簡化
		Matrix tmp = mtx.reduceRowColumn(0, i);
		result = result +  mtx.get_data(0,i) * std::pow(-1, i) * Matrix::DetRecursive(tmp);
	}
	return result;
}

double Matrix::DetRecursive(int m, int n, double * data)
{
	double result = 0;
	if (m == 1 && n == 1) {
		//std::cout << mtx.get_data(0, 0) << std::endl;
		result = data[0];
		//std::cout << result << std::endl;
		return result;//簡化到matrix中只有一個element就return該值
	}
	for (int i = 0; i < m; i++) {
		//重複取新的matrix的第一列去簡化
		double *tmp = StdVectorArrayTransform::ArrayReduceRowColumn(0, i, m, n, data);
		result = result +data[i] * std::pow(-1, i) * Matrix::DetRecursive((m-1), (n-1), tmp);
	}
	return result;
}

Matrix operator*(double constant, Matrix &mtx) {
	std::vector<std::vector<double> > temp(mtx.RowSize, std::vector<double>(mtx.ColumnSize));
	for (int i = 0; i < mtx.RowSize; i++) {
		for (int j = 0; j < mtx.ColumnSize; j++) {
			temp[i][j] = constant * mtx.data[i][j];
		}
	}
	return Matrix(mtx.RowSize, mtx.ColumnSize, temp);
}
