#include "vector.h"

#define PI  3.141592653589793238463

Vector::Vector()
{
	this->Size = 0;
}

Vector::Vector(int len, std::vector<double> vtr_data)
{
	this->Size = len;
	this->data.assign(vtr_data.begin(), vtr_data.end());
}

int Vector::get_m()
{
	return this->Size;
}

double Vector::get_data(int m)  throw (VectorException)
{
	if (m >= this->Size) {
		throw VectorException();
	}
	return this->data[m];
}

void Vector::setData(int len, std::vector<double> vtr_data)
{
	this->Size = len;
	this->data.assign(vtr_data.begin(), vtr_data.end());
}

void Vector::show()
{
	for (int i = 0; i < this->Size; i++) {
		std::cout << this->data[i] << " ";
	}
	std::cout << std::endl;
}

Vector Vector::operator+(Vector vtr) throw(VectorException)
{
	if (this->Size != vtr.Size) {
		throw VectorException();
	}
	std::vector<double> temp(Size);
	for (int i = 0; i < this->Size; i++) {
		temp[i] = this->data[i] + vtr.data[i];
	}
	return Vector(this->Size, temp);
}

double Vector::operator*(Vector vtr) throw(VectorException)
{
	if (this->Size != vtr.Size) {
		throw VectorException();
	}
	double result = 0;
	for (int i = 0; i < this->Size; i++) {
		result += this->data[i] * vtr.data[i];
	}
	return result;
}

Vector Vector::operator*(double constant)
{
	std::vector<double> temp(this->Size);
	for (int i = 0; i < this->Size; i++) {
		temp[i] = this->data[i] * constant;
	}
	return Vector(this->Size, temp);
}

Vector operator*(double constant, Vector vtr)
{
	std::vector<double> temp(vtr.Size);
	for (int i = 0; i < vtr.Size; i++) {
		temp[i] = vtr.data[i] * constant;
	}
	return Vector(vtr.Size, temp);
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
	if (this->Size != vtr.Size) {
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

double Vector::Component(Vector A, Vector B) throw (VectorException)
{
	if (A.Size != B.Size) throw VectorException();
	return (A * B) / B.length();
}

Vector Vector::Projection(Vector A, Vector B)
{
	return ((A * A) / (B * B)) * B;
}

Matrix Vector::parsetoMatrix(std::vector<Vector> Vectors, char type='R') throw (VectorException)
{
	int SizeofVector = Vectors[0].Size;
	for (int i = 1; i < Vectors.size(); i++) {
		if (SizeofVector != Vectors[i].Size) throw VectorException();
	}
	if (type == 'R') {
		std::vector<std::vector<double> > tmp(Vectors.size(), std::vector<double>(SizeofVector));
		for (int i = 0; i < Vectors.size(); i++) {
			for (int j = 0; j < SizeofVector; j++) {
				tmp[i][j] = Vectors[i].data[j];
			}
		}
		return Matrix(Vectors.size(), SizeofVector, tmp);
	}
	else if (type == 'C') {
		std::vector<std::vector<double> > tmp(SizeofVector, std::vector<double>(Vectors.size()));
		for (int i = 0; i < Vectors.size(); i++) {
			for (int j = 0; j < SizeofVector; j++) {
				tmp[j][i] = Vectors[i].data[j];
			}
		}
		return Matrix(SizeofVector, Vectors.size(), tmp);
	}
	else {
		throw VectorException();
	}
}






