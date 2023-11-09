#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

#include "ManageData.h"

ManageData::ManageData(){
	_ncols = 0; _nrows = 0; _nf = 0;
	_time = 0; _nmin = 0;
}

ManageData::ManageData(string name){
	_name = name;
	_get_params();
	_number_columns();
	_number_rows();
}

void ManageData::_get_params(){
	ifstream file;
	file.open(_name);

	string line;
	stringstream t;
	getline(file, line);
	t << line;
	t >> _file_name;


	getline(file, line);

	stringstream s;
	s << line;
	s >> _time;
	s >> _nf;
	s >> _nmin;

	file.close();
}

void ManageData::_number_columns(){
	// Return the number of columns of the data

	_ncols = 0;

	ifstream file;
	file.open(_name);

	string line;
	getline(file, line);
	getline(file, line); // jump the first 2 lines since it is only used for saving time
	getline(file, line);

	stringstream s;
	s << line;
	float value; 
	while(s >> value) _ncols++;

	file.close();
}

void ManageData::_number_rows(){
	// Return the number of rows of the data

	_nrows = 0;

	ifstream file;
	file.open(_name);

	string line;
	while (getline(file, line)) _nrows++;

	file.close();
}

float ManageData::operator()(int row, int col){

	ifstream file;
	file.open(_name);

	string line;
	stringstream s;
	float value, aux;
	int i = 0;
	int j = 0;

	row+=2;

	if (row < _nrows && col < _ncols && row >= 1 && col >= 0){
		while (file){
			getline(file, line);
			if (i == row){
				s << line;
				while(s >> aux){
					if (j == col) value = aux;
					j++;
				}
			}
			i++;
		}
	}
	else{
		cout << "Index out of range" << endl;
		exit(1);
	}

	file.close();

	return value;
}
