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
	//operator overloading
	Matrix operator+(Matrix &mtx) throw (MatrixException);
	Matrix operator*(Matrix &mtx) throw (MatrixException);
	Matrix operator*(double constant);
	friend Matrix operator*(double constant, Matrix &mtx);
	//complex
	std::vector<Matrix> LU() throw (MatrixException);
	Matrix reduce(int m, int n) throw (MatrixException);//消除某一列行
	double determinant();
	double determinant_A();
	double multiDiagonal();
	Matrix inverse();
	Matrix gauss();

private:
	int RowIndex;
	int ColumnIndex;
	double **data_a;
	std::vector<std::vector<double> > data;
	void setData_a(int m, int n, double **mtx_data);
	static double determinant_recursive(Matrix &mtx);
	static double determinant_recursive(int m, int n, double *data);
};
