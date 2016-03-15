#include "DataManager.h"
#include "matrix.h"
#include "vector.h"

DataManager::DataManager()
{
	this->VectorVariableIndex = 0;
	this->MatrixVariableindex = 0;
}

bool DataManager::LoadData()
{
	std::fstream fin;
	//開啟檔案，傳入open函數的參數有兩個，欲開起的檔案名稱，開啟檔案的模式參數(這邊std::ios::in為讀取(輸入)狀態)
	fin.open(FileName, std::ios::in);
	this->Matrixs.clear();
	this->Vectors.clear();
	//讀取失敗回傳false
	if (!fin)
	{
		return false;
	}
	else
	{
		//標記當前讀取向量ID
		int currentLoadVectorID = 0;
		//定義向量資料暫存變數
		std::vector<double> tempVectorData;
		//定義讀取檔案字串暫存變數
		std::string tempString;
		//從檔案讀取字串，解析掉向量總數
		fin >> tempString;
		int totalinput = std::stoi(tempString);
		//建立變數名稱清單
		this->VariableNameList = std::vector<std::string>(totalinput);
		std::string InputType;

		//執行讀檔迴圈，並在讀到檔案結尾時結束
		while (currentLoadVectorID < totalinput)
		{
			fin >> InputType;
			if (InputType[0] == 'V') {
				//加入變數名稱清單
				//VariableNameList[currentLoadVectorID] = "$v" + std::to_string(currentLoadVectorID);
				//建立足夠數量的Vector
				//this->Vectors = std::vector<Vector>(totalinput);
				//取得Vector的長度
				fin >> tempString;
				int len = std::stoi(tempString);
				//建立暫存data的std::vector<double>
				std::vector<double> vtr_data(len);
				//讀取所有資料
				for (int i = 0; i < len; i++) {
					fin >> tempString;
					vtr_data[i] = std::stod(tempString);
				}
				Vectors.push_back(Vector(len, vtr_data));
				currentLoadVectorID++;
			}
			else if (InputType[0] == 'M') {
				//加入變數名稱清單
				//VariableNameList[currentLoadVectorID] = "$m" + std::to_string(currentLoadVectorID);
				//建立足夠數量的Vector
				//this->Matrixs = std::vector<Matrix>(totalinput);
				//取得Matrix的row和column
				fin >> tempString;
				int m = std::stoi(tempString);
				fin >> tempString;
				int n = std::stoi(tempString);
				//建立暫存matrix data std::vector<std::vector<double> >
				std::vector<std::vector<double> > mtx_data(m, std::vector<double>(n));
				//讀取所有資料
				for (int i = 0; i < m; i++) {
					for (int j = 0; j < n; j++) {
						fin >> tempString;
						mtx_data[i][j] = std::stod(tempString);
					}
				}
				this->Matrixs.push_back(Matrix(m, n, mtx_data));
				currentLoadVectorID++;
			}
			else {
			}
		}
		return true;
	}
}

void DataManager::SetFileName(std::string FileName)
{
	this->FileName = FileName;
}

Vector DataManager::GetVector(int index)
{
	return this->Vectors[index];
}

Matrix DataManager::GetMatrix(int index)
{
	return this->Matrixs[index];
}
