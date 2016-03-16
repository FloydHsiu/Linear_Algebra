#pragma once
#include <iostream>
#include "matrix.h"
#include "vector.h"
#include "DataManager.h"
#include <vector>
#include <string.h>
#include <fstream>


int main() {
	DataManager dataManager;
	dataManager.SetFileName("C:/Users/Floyd/Documents/Visual Studio 2015/Projects/Project1/TestData/Matrix/M5.txt");
	//dataManager.SetFileName("/Users/Floyd/Downloads/project1/TestData/Matrix");
	dataManager.LoadData();
	Matrix::solveLinearEquation(dataManager.GetMatrix(0), dataManager.GetMatrix(1)).show();
	Matrix::solveLinearEquation(dataManager.GetMatrix(2), dataManager.GetMatrix(3)).show();
	Matrix::solveLinearEquation(dataManager.GetMatrix(4), dataManager.GetMatrix(5)).show();
	/*int m = 0;
	Matrix tmp = dataManager.GetMatrix(1);
	for (int i = 0; i < m; i++) {
		tmp = tmp.reduceRowColumn(0, 0);
	}
	tmp.show();*/
	
	//std::cout << dataManager.GetMatrix(3).determinant() << std::endl;
	//std::cout << dataManager.GetMatrix(1).get_m();
}