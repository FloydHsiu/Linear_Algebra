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
	//�}���ɮסA�ǤJopen��ƪ��ѼƦ���ӡA���}�_���ɮצW�١A�}���ɮת��Ҧ��Ѽ�(�o��std::ios::in��Ū��(��J)���A)
	fin.open(FileName, std::ios::in);
	this->Matrixs.clear();
	this->Vectors.clear();
	//Ū�����Ѧ^��false
	if (!fin)
	{
		return false;
	}
	else
	{
		//�аO��eŪ���V�qID
		int currentLoadVectorID = 0;
		//�w�q�V�q��ƼȦs�ܼ�
		std::vector<double> tempVectorData;
		//�w�qŪ���ɮצr��Ȧs�ܼ�
		std::string tempString;
		//�q�ɮ�Ū���r��A�ѪR���V�q�`��
		fin >> tempString;
		int totalinput = std::stoi(tempString);
		//�إ��ܼƦW�ٲM��
		this->VariableNameList = std::vector<std::string>(totalinput);
		std::string InputType;

		//����Ū�ɰj��A�æbŪ���ɮ׵����ɵ���
		while (currentLoadVectorID < totalinput)
		{
			fin >> InputType;
			if (InputType[0] == 'V') {
				//�[�J�ܼƦW�ٲM��
				//VariableNameList[currentLoadVectorID] = "$v" + std::to_string(currentLoadVectorID);
				//�إߨ����ƶq��Vector
				//this->Vectors = std::vector<Vector>(totalinput);
				//���oVector������
				fin >> tempString;
				int len = std::stoi(tempString);
				//�إ߼Ȧsdata��std::vector<double>
				std::vector<double> vtr_data(len);
				//Ū���Ҧ����
				for (int i = 0; i < len; i++) {
					fin >> tempString;
					vtr_data[i] = std::stod(tempString);
				}
				Vectors.push_back(Vector(len, vtr_data));
				currentLoadVectorID++;
			}
			else if (InputType[0] == 'M') {
				//�[�J�ܼƦW�ٲM��
				//VariableNameList[currentLoadVectorID] = "$m" + std::to_string(currentLoadVectorID);
				//�إߨ����ƶq��Vector
				//this->Matrixs = std::vector<Matrix>(totalinput);
				//���oMatrix��row�Mcolumn
				fin >> tempString;
				int m = std::stoi(tempString);
				fin >> tempString;
				int n = std::stoi(tempString);
				//�إ߼Ȧsmatrix data std::vector<std::vector<double> >
				std::vector<std::vector<double> > mtx_data(m, std::vector<double>(n));
				//Ū���Ҧ����
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
