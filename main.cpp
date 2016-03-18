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
	dataManager.SetFileName("C:/Users/Floyd/Documents/Visual Studio 2015/Projects/Project1/TestData/Vector/V14.txt");
	//dataManager.SetFileName("/Users/Floyd/Downloads/project1/TestData/Matrix");
	dataManager.LoadData();
	Vector::Cross(dataManager.GetVector(0), dataManager.GetVector(1)).show();
	Vector::Cross(dataManager.GetVector(2), dataManager.GetVector(3)).show();
	/*int m = 0;
	Matrix tmp = dataManager.GetMatrix(1);
	for (int i = 0; i < m; i++) {
		tmp = tmp.reduceRowColumn(0, 0);
	}
	tmp.show();*/
	
	//std::cout << dataManager.GetMatrix(3).determinant() << std::endl;
	//std::cout << dataManager.GetMatrix(1).get_m();
}