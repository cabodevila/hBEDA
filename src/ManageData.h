/*
 * Reads data from file (in a very unefficient way)
*/

#define MANAGE_DATA_H

#include <string>

using namespace std;

class ManageData{

private:
	string _name;
	unsigned int _ncols, _nrows, _nf, _nmin;
	float _time;
	string _file_name;

	void _number_columns();
	void _number_rows();
	void _get_params();

public:
	ManageData();
	ManageData(string);

	int number_columns(){return _ncols;};
	int number_rows(){return _nrows-2;};
	float get_time(){return _time;};
	float get_nf(){return _nf;};
	float get_nmin(){return _nmin;};
	string get_name(){return _file_name;};


	float operator()(int, int); // Arguments: row, column

};