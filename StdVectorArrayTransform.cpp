#include "StdVectorArrayTransform.h"

double* StdVectorArrayTransform::parseToArray(int m, int n, std::vector<std::vector<double>> data)
{
	double *result = new double[m*n];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			result[i*m + j] = data[i][j];
		}
	}
	return result;
}

double * StdVectorArrayTransform::ArrayReduceRowColumn(int reducedRow, int reducedColumn, int m, int n, double * data)
{
	double *result = new double[(m - 1)*(n - 1)];
	int k = 0;
	for (int i = 0; i < m-1; i++) {
		int l = 0;
		if (i == reducedRow) k++;
		for (int j = 0; j < n-1; j++) {
			if (j == reducedColumn) l++;
			result[i*(m - 1) + j] = data[k*m + l];
			l++;
		}
		k++;
	}
	return result;
}
