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
	//dataManager.GetMatrix(2).show();
	//dataManager.GetMatrix(3).show();
	/*std::fstream fin;
	std::string temp;
	std::string FileName = "C:/Users/Floyd/Documents/Visual Studio 2015/Projects/Project1/TestData/Matrix/m6.txt";
	//開啟檔案，傳入open函數的參數有兩個，欲開起的檔案名稱，開啟檔案的模式參數(這邊std::ios::in為讀取(輸入)狀態)
	fin.open(FileName, std::ios::in);
	for (int i = 0; i < 20; i++) {
		fin >> temp;
		std::cout <<i << std::endl <<temp << std::endl;
	}*/
}