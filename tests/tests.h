#ifndef TEST_PARASITE_H
#define TEST_PARASITE_H

#include "../ParaSite/ParaSite.h"
#include <iostream>
#include <string>

using namespace MPI_Wrapper;
using namespace ParaSite;

//MPIWrapper2D class tests
int Test_Build_Parallel2D(const int n_rows, const int n_cols);
int Test_Build_2DGrid(const int n_rows, const int n_cols);
int Test_Broadcast(const int n_rows, const int n_cols);
int Test_Sum(const int n_rows, const int n_cols);
int Test_Max(const int n_rows, const int n_cols);
int Test_Min(const int n_rows, const int n_cols);
int Test_SendRecv(const int n_rows, const int n_cols);
int Test_SendUp(const int n_rows, const int n_cols);
int Test_SendDown(const int n_rows, const int n_cols);

int Test_SerialOperation();

//Lattice class tests
int Test_Build_Lattice3D(
				const IndexType size[3], 
                const int halo, 
                const int node_size[2], 
                const int grid_loc[2]);

int Test_Lattice3D_Global_Coordinate_Transformation();
int Test_Local_Visible_Loop_Single_Process(
	const int rounds = 1, 
	const bool use_hash_table = false);
int Test_Local_Visible_Loop_Pseudo_Multi_Processes();

//MPI-based tests: Site class
int Test_Local_Visible_Loop(const int n_rows, const int n_cols);
int Test_Halo_Loop(const int n_rows, const int n_cols);
int Test_Site_Move(const int n_rows, const int n_cols);
int Test_Set_Get_Coordinates(const int n_rows, const int n_cols);

//MPI-based tests: Field class
int Test_UpdateHalo(const int n_rows, const int n_cols);

//HDF5 tests
int Test_HDF5_WriteDataset();
int Test_ParaSiteFieldWrite(const int n_rows, const int n_cols);

//helper functions -------------------
/*
compare whether theo_val == exp_val
*/
template<class Type>
int compare_results(const Type theo_val, const Type exp_val, 
					const std::string& name, 
					const bool output = true) {
	if (theo_val == exp_val) {
		if(output){
			std::cout << name << " match. value = " << theo_val 
				<< "." << std::endl;
		}
		return 1;
	}
	else {
		if(output){
			std::cout << name << " UNMATCH. Theoretical value = " << theo_val
				<< " ; Experimental value = " << exp_val << "." << std::endl;
		}
		return 0;
	}
}

template<class Type>
int compare_results(const Type* theo_val, const Type* exp_val, const int len, const std::string& name){
	auto status = 1;
	for(auto i = 0; i < len; ++i){
		status = (theo_val[i] == exp_val[i]) && status;
	}
	if(status){
		std::cout << name << " match." << std::endl;
	}
	else{
		std::cout << name << " UNMATCH. {";
		for(auto i = 0; i<len; ++i){
			std::cout << theo_val[i] << " : " << exp_val[i] <<", ";
		}
		std::cout << "}" << std::endl;
	}
	return status;
}

template<class Type>
void print_coordinate(const Type* c, const int len, 
		const std::string& name, 
		const bool end_line = true){
	std::cout << name <<": (";
	for(auto i = 0; i < len; ++i){
		std::cout << c[i] << ",";
	}
	std::cout << ") ";
	if(end_line) std::cout << std::endl;
}

#endif // !TEST_MPI_WRAPPER_2D
