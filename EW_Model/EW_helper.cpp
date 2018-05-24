#include "EW_helper.h"
#include <sstream>
#include <iomanip>

std::string ReadConfigIni(const std::string& key, const std::string& file_path) {
	std::ifstream fin;
	fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	std::string val = "";
	std::string curr_line;
	try {
		fin.open(file_path, std::ifstream::in);
	}
	catch (std::system_error& e) {
		std::cerr << e.code().message() << std::endl;
	}
	if (fin.is_open()) {
		while (getline(fin, curr_line))
		{
			int pos = curr_line.find(" = ");
			std::string curr_key = curr_line.substr(0, pos);
			if (curr_key == key) {
				return curr_line.substr(pos + 3);
			}
		}
		fin.close();
	}
	else {
		std::cerr << "config.ini DOES NOT EXIST." << std::endl;
		getchar();
		exit(EXIT_FAILURE);
	}
	return val;
}

/*DataTable class*/

//public functions
DataTable::DataTable()
{}

DataTable::DataTable(const std::vector<std::string>& keys)
{
	initialize(keys);
}

DataTable::~DataTable(){
	for(auto& item : val_ptr_){
		delete item;
	}
}

void DataTable::initialize(const std::vector<std::string>& keys) {	
	columns_ = keys;
	val_ptr_.resize(columns_.size());
	for(auto& item : val_ptr_){
		item = new std::vector<double>();
	}
	return;
}

void DataTable::append_value(const std::string& key, const double val) const {
	const int col = get_column_index(key);
	if (col != -1) {
		this->val_ptr_[col]->push_back(val);
	}
	else {
		std::cerr << "Column does not exist." << std::endl;
	}
	return;
}

void DataTable::append_value(const int index, const double val) const {
	this->val_ptr_[index]->push_back(val);
	return;
}

void DataTable::save_table(const std::string& fname) const {
	std::ofstream fout;
	fout.open(fname, std::ios::trunc);
	if (fout.is_open()) {
		//const std::string first = "DataCount";
		//fout << first;
		for (auto& iter : this->columns_) {
			fout << iter << ",";
		}
		fout << std::endl;

		if (is_table_complete()) {
			const auto cols_size = val_ptr_.size();
			const auto common_size = val_ptr_[0]->size();
			for (auto count = 0; count < common_size; ++count) {
				for (auto l = 0; l < cols_size; ++l) {
					fout << (*val_ptr_[l])[count] << ",";
				}
				fout << std::endl;
			}
		}
		else{
			std::cerr << "Table not complete." << std::endl;
		}
		fout.close();
	}
	else {
		std::cerr << "ERROR IN OPENING FILE. (CLASS DataTable)" << std::endl;
	}
	return;
}

const std::vector<double>& DataTable::get_column(const std::string& key) const {
	auto col_index = this->get_column_index(key);
	if (col_index == -1) {
		std::cerr << "Key does not exist." << std::endl;
		return *(this->val_ptr_[0]);
	}
	else {
		return *(this->val_ptr_[col_index]);
	}
}

int DataTable::get_column_index(const std::string& key) const {
	auto iter = std::find(columns_.begin(), columns_.end(), key);
	if (iter != columns_.end()) {
		return static_cast<int>(iter - columns_.begin());
	}
	else {
		return -1;
	}
}

bool DataTable::is_table_complete() const {
	const auto common_size = (*val_ptr_[0]).size();
	for (auto& iter : val_ptr_) {
		if (iter->size() != common_size) return false;
	}
	return true;
}

/*Parameters class*/

Parameters::Parameters()
	:
	params_()
{}

void Parameters::add(const std::string& key, const double val, const bool is_scientific) {
	std::stringstream ss;
	if(is_scientific)
		ss << std::scientific << std::setprecision(6) << val; 
	else 
		ss << std::setprecision(6) << val;
	add(key, ss.str());
	return;
}

void Parameters::add(const std::string& key, const int val){
	add(key, std::to_string(val));
	return;
}

void Parameters::add(const std::string& key, const std::string& val){
	params_.push_back(std::make_pair(key, val));
	return;
}

void Parameters::set_delimiter(const std::string& delimiter){
	this->delimiter_ = delimiter;
	return;
}

void Parameters::save(const std::string& fname) const {
	std::ofstream fout;
	fout.open(fname, std::ios::trunc);
	if (fout.is_open()) {
		for (auto& item : this->params_) {
			fout << item.first 
				<< " " << this->delimiter_ << " "
				<< item.second << std::endl;
		}
		fout.close();
	}
	else {
		std::cerr << "ERROR IN OPENING FILE. (CLASS Parameters)" << std::endl;
	}
	return;
}