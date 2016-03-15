#pragma once
#include <vector>
#include <iostream>

class StdVectorArrayTransform 
{
public:
	static double* parseToArray(int m, int n, std::vector<std::vector<double> > data);
	static double* ArrayReduceRowColumn(int reducedRow, int reducedColumn, int m, int n, double* data);
};