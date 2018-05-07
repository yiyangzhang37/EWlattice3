#ifndef TEST_MPI_WRAPPER_2D
#define TEST_MPI_WRAPPER_2D

#include "MPIWrapper2D.h"
#include "LATfield3_Lattice.h"
#include <iostream>
#include <string>

using namespace MPI_Wrapper;
using namespace LATfield3;

int Test_Build_Parallel2D(const int n_rows, const int n_cols);
int Test_Build_2DGrid(const int n_rows, const int n_cols);
int Test_Broadcast(const int n_rows, const int n_cols);
int Test_Sum(const int n_rows, const int n_cols);
int Test_Max(const int n_rows, const int n_cols);
int Test_Min(const int n_rows, const int n_cols);
int Test_SendRecv(const int n_rows, const int n_cols);
int Test_SendUp(const int n_rows, const int n_cols);
int Test_SendDown(const int n_rows, const int n_cols);

int Test_Build_Lattice3D(
				const IndexType size[3], 
                const int halo, 
                const int node_size[2], 
                const int grid_loc[2]);

int Test_Lattice3D_Global_Coordinate_Transformation();
int Test_Local_Visible_Loop_Single_Process();
int Test_Local_Visible_Loop_Pseudo_Multi_Processes();

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

#endif // !TEST_MPI_WRAPPER_2D
