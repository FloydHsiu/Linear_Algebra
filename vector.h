#pragma once

#include <iostream>
#include <vector>
#include <cmath>

class VectorException
{
public:
	VectorException() {
		this->msg = "Error";
	};
	char* log() {
		return msg;
	};
private:
	char * msg;
};

class Vector
{
public:
	Vector();
	Vector(int len, std::vector<double> vtr_data);
	int get_m();
	double get_data(int m) throw (VectorException);
	void setData(int len, std::vector<double> vtr_data);
	void show();
	//operator overloading
	Vector operator+(Vector vtr) throw (VectorException);
	double operator*(Vector vtr) throw (VectorException);
	Vector operator*(double constant);
	friend Vector operator*(double constant, Vector vtr);
	//essential
	double length();
	Vector normalization();
	double angle(Vector &vtr) throw(VectorException);
	double area(Vector&vtr) throw(VectorException);
	bool isParallel(Vector &vtr) throw (VectorException);
	bool isOrthogonal(Vector &vtr) throw (VectorException);

private:
	std::vector<double> data;
	int _m;
};
