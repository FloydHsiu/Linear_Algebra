#pragma once
#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <exception>
#include "StdVectorArrayTransform.h"

class MatrixException
{
public:
	MatrixException() {
		this->msg = "Error";
	};
	void log() {
		std::cout << msg << std::endl;
	};

private:
	char *msg;
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
	Matrix adjoint();
	Matrix inverse() throw(MatrixException);
	Matrix gauss();//高斯消去法
	double Det_Gauss();
	double Det_RecursiveAndVector();
	double Det_RecursiveAndArray();
	int rank();
	static Matrix solveLinearEquation(Matrix A, Matrix B) throw (MatrixException);// AX = B --> X = inverse(A)*B

private:
	int RowIndex;
	int ColumnIndex;
	std::vector<std::vector<double> > data;
	static double DetRecursive(Matrix &mtx);
	static double DetRecursive(int m, int n, double *data);
};
