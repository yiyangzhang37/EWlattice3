#ifndef EW_HELPER_H
#define EW_HELPER_H

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<system_error>
#include<algorithm>
#include<map>

#ifdef _WIN32
const std::string CONFIG_PATH = ".\\config.ini";
#else
const std::string CONFIG_PATH = "./config.ini";
#endif

std::string ReadConfigIni(const std::string& item,
	const std::string& file_path = CONFIG_PATH);



class DataTable {
public:
	DataTable();
	DataTable(const std::vector<std::string>& keys);
	~DataTable();
	//add a value to the current table.
	//It assumes the previous data are all in order, and will assume itself in the right position.
	void initialize(const std::vector<std::string>& keys);
	void append_value(const std::string& key, const double val) const;
	void append_value(const int index, const double val) const;
	void modify_value(const std::string& key, const double val, const int row);
	void save_table(const std::string& fname) const;

	//get the index of the key in the this->columns_
	//returns -1 if the key does not exist in the this->columns_
	int get_column_index(const std::string& key) const;
	const std::vector<double>& get_column(const std::string& key) const;

	bool is_table_complete() const;
private:
	std::vector< std::string > columns_;
	std::vector< std::vector<double>* > val_ptr_;
};



class Parameters {
public:
	Parameters();
	void add(const std::string& key, const double val);
	void add(const std::string& key, const int val);
	void add(const std::string& key, const std::string& val);
	template<class Type>
	void add(const std::string& key, const Type val);
	void set_delimiter(const std::string& delimiter);
	void save(const std::string& fname) const;
private:
	std::vector< std::pair<std::string, std::string> > params_;
	std::string delimiter_ = "=";
};

template<class Type>
void Parameters::add(const std::string& key, const Type val) {
	add(key, std::to_string(val));
	return;
}


#endif