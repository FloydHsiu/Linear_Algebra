#pragma once
#include "vector.h"
#include "matrix.h"
#include<vector>
#include<iostream>
#include<fstream>
#include<string>


class DataManager 
{
private:
	std::vector<Vector> Vectors;
	std::vector<Matrix> Matrixs;
	int VectorVariableIndex;
	int MatrixVariableindex;
	std::string FileName;
	std::vector<std::string> VariableNameList;

public:
	DataManager();
	bool LoadData();
	void SetFileName(std::string FileName);
	Vector GetVector(int index);
	Matrix GetMatrix(int index);
};