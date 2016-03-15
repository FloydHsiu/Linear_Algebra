#include "vector.h"

#define PI  3.141592653589793238463

Vector::Vector()
{
	this->_m = 0;
}

Vector::Vector(int len, std::vector<double> vtr_data)
{
	this->_m = len;
	this->data.assign(vtr_data.begin(), vtr_data.end());
}

int Vector::get_m()
{
	return this->_m;
}

double Vector::get_data(int m)  throw (VectorException)
{
	if (m >= this->_m) {
		throw VectorException();
	}
	return this->data[m];
}

void Vector::setData(int len, std::vector<double> vtr_data)
{
	this->_m = len;
	this->data.assign(vtr_data.begin(), vtr_data.end());
}

void Vector::show()
{
	for (int i = 0; i < this->_m; i++) {
		std::cout << this->data[i] << " ";
	}
	std::cout << std::endl;
}

Vector Vector::operator+(Vector vtr) throw(VectorException)
{
	if (this->_m != vtr._m) {
		throw VectorException();
	}
	std::vector<double> temp(_m);
	for (int i = 0; i < this->_m; i++) {
		temp[i] = this->data[i] + vtr.data[i];
	}
	return Vector(this->_m, temp);
}

double Vector::operator*(Vector vtr) throw(VectorException)
{
	if (this->_m != vtr._m) {
		throw VectorException();
	}
	double result = 0;
	for (int i = 0; i < this->_m; i++) {
		result += this->data[i] * vtr.data[i];
	}
	return result;
}

Vector Vector::operator*(double constant)
{
	std::vector<double> temp(this->_m);
	for (int i = 0; i < this->_m; i++) {
		temp[i] = this->data[i] * constant;
	}
	return Vector(this->_m, temp);
}

Vector operator*(double constant, Vector vtr)
{
	std::vector<double> temp(vtr._m);
	for (int i = 0; i < vtr._m; i++) {
		temp[i] = vtr.data[i] * constant;
	}
	return Vector(vtr._m, temp);
}

double Vector::length()
{
	double result;
	//do dot with self
	result = (*this) * (*this);
	return sqrt(result);
}

Vector Vector::normalization()
{
	return ((1 / this->length()) * (*this));
}

double Vector::angle(Vector & vtr) throw(VectorException)
{
	double cos_ = (*this * vtr) / (this->length()*vtr.length());
	return (acos(cos_)/PI) * 180;
}

double Vector::area(Vector & vtr) throw(VectorException)
{
	return 0.5 * this->length() * vtr.length() * sin(this->angle(vtr) / 180 * PI);
}

bool Vector::isParallel(Vector & vtr) throw (VectorException)
{
	if (this->_m != vtr._m) {
		throw VectorException();
	}
	if (this->angle(vtr) < 1e-6) return true;
	else return false;
}

bool Vector::isOrthogonal(Vector & vtr) throw(VectorException)
{
	if ((*this) * vtr < 1e-6) return true;
	else return false;
}


