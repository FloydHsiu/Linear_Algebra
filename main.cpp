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
	dataManager.SetFileName("C:/Users/Floyd/Documents/Visual Studio 2015/Projects/Project1/TestData/Matrix/m6.txt");
	//dataManager.SetFileName("/Users/Floyd/Downloads/project1/TestData/Matrix");
	dataManager.LoadData();
	int m = 0;
	Matrix tmp = dataManager.GetMatrix(1);
	for (int i = 0; i < m; i++) {
		tmp = tmp.reduce(0, 0);
	}
	tmp.show();
	std::cout << tmp.gauss().multiDiagonal() << std::endl;
	std::cout << dataManager.GetMatrix(0).gauss().multiDiagonal() << std::endl;
	std::cout << dataManager.GetMatrix(2).gauss().multiDiagonal() << std::endl;
	//std::cout << dataManager.GetMatrix(3).determinant() << std::endl;
	//std::cout << dataManager.GetMatrix(1).get_m();
}