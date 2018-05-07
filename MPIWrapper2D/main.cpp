#include "MPIWrapper2D.h"
#include "Test_MPIWrapper2D.h"
#include "LATfield3_Lattice.h"
#include <iostream>
#include <thread>
#include <chrono>

int main(int argc, char** argv) {
	Parallel_Init();

	int n_rows = 1, n_cols = 1;
	for (int i = 1; i < argc; i++) {
		if (argv[i][0] != '-')
			continue;
		switch (argv[i][1]) {
		case 'r':
			n_rows = atoi(argv[++i]);
			break;
		case 'c':
			n_cols = atoi(argv[++i]);
			break;
		}
	}
	int world_rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int status = 1;
	/*
	status = Test_Build_Parallel2D(n_rows, n_cols) && status;

	status = Test_Build_2DGrid(n_rows, n_cols) && status;

	status = Test_Broadcast(n_rows, n_cols) && status;

	status = Test_Sum(n_rows, n_cols) && status;

	status = Test_SendRecv(n_rows, n_cols) && status;

	status = Test_SendUp(n_rows, n_cols) && status;

	status = Test_SendDown(n_rows, n_cols) && status;
	*/
	unsigned int lat_size[] = {11,16,23};
	int node_size[] = {1,1};
	int loc[] = {0,0};
	status = Test_Build_Lattice3D(lat_size, 2, node_size, loc) && status;

	status = Test_Lattice3D_Global_Coordinate_Transformation() && status;

	status = Test_Local_Visible_Loop_Single_Process() && status;

	status = Test_Local_Visible_Loop_Pseudo_Multi_Processes() && status;

	// summerize results
	int overall_status = 1;
	MPI_Reduce(&status, &overall_status, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank == 0){
		if(overall_status){
			std::cout << "TEST SUCCESS." << std::endl;
		}
		else{
			std::cout << "TEST FAIL." << std::endl;
		}
	}
	Parallel_Finalize();
	return !status;
} 