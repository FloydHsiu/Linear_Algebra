#pragma once
#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <exception>
#include <string>
#include "StdVectorArrayTransform.h"

class MatrixException
{
public:
	MatrixException() {
		this->msg = "Error";
	};
	MatrixException(std::string department, char type) {
		this->department = department;
		switch (type) {
		case 'A'://A
		case 'a':
			this->msg = "Access Matrix position is out of range.";
			break;
		case 'D'://D
		case 'd':
			this->msg = "Two different size Matrix error.";
			break;
		case 'M'://M
		case 'm':
			this->msg = "Matrix multiple Matrix error.";
			break;
		case 'S'://S
		case 's':
			this->msg = "No a square Matrix error.";
			break;
		case 'L':
		case 'l':
			this->msg = "Unable to do LU decomposition error.";
			break;
		case 'I':
		case 'i':
			this->msg = "Determinant is zero error.";
			break;
		default:
			this->msg = "Error";

		}
	};
	void log() {
		std::cout << msg << std::endl;
	};

private:
	std::string msg;
	std::string department;
};

class Matrix
{
public:
	Matrix();
	Matrix(int m, int n, std::vector<std::vector<double> > mtx_data);
	~Matrix();
	int get_RowIndex();
	int get_ColumnIndex();
	double get_data(int m, int n) throw (MatrixException);
	void setData(int m, int n, std::vector<std::vector<double> > mtx_data);
	void show();
	//essential
	static Matrix transpose(Matrix mtx);
	Matrix reduceRowColumn(int m, int n) throw (MatrixException);//消除某一列行
	double multiDiagonal();
	//operator overloading
	Matrix operator+(Matrix &mtx) throw (MatrixException);
	Matrix operator-(Matrix &mtx) throw(MatrixException);
	Matrix operator*(Matrix &mtx) throw (MatrixException);
	Matrix operator*(double constant);
	friend Matrix operator*(double constant, Matrix &mtx);
	//complex
	std::vector<Matrix> LU() throw (MatrixException);
	double cofactor(int m, int n);
	Matrix adjoint() throw(MatrixException);
	Matrix inverse() throw(MatrixException);
	Matrix gauss();//高斯消去法
	double Det_Gauss();
	double Det_RecursiveAndVector();
	double Det_RecursiveAndArray();
	int rank();
	static Matrix solveLinearEquation(Matrix A, Matrix B) throw (MatrixException);// AX = B --> X = inverse(A)*B
	static Matrix LeastSquare(Matrix A, Matrix B) throw (MatrixException);//tanspose(A)AX = tanspose(A)B  --> X = inverse(tansepose(A)A)transpose(A)B 

private:
	int RowSize;
	int ColumnSize;
	std::vector<std::vector<double> > data;
	static double DetRecursive(Matrix &mtx);
	static double DetRecursive(int m, int n, double *data);
};
