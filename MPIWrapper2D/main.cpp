//#include "./tests/Test_MPIWrapper2D.h"
#include <iostream>
#include <thread>
#include <chrono>
#include <array>
#include "EW_parameter.h"
#include "EW_helper.h"

#include <hdf5_hl.h>

#include "./ParaSite/ParaSite.h"

#include "EW_Model.h"
#include "EW_Observation.h"
//#include "EW_BubbleNucl.h"

using namespace MPI_Wrapper;
using namespace ParaSite;
using namespace Electroweak;
using namespace HDF5_Wrapper;

int main(int argc, char** argv) {
	Parallel_Init();
	{
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
		/*
		Parallel2D parallel(n_rows, n_cols, MPI_COMM_WORLD);
		GridIndexType node_size[] = { n_rows, n_cols };
		GridIndexType grid_rank[] = { parallel.get_grid_rank()[0], parallel.get_grid_rank()[1] };

		GridIndexType grid_loc[2];
		transform_gridrank_to_gridloc(grid_rank, grid_loc);
		IndexType lat_size[3] = {12,24,36};
		Lattice<3> lat(lat_size, 2, node_size, grid_loc);
		*/
		/*
		ElectroweakEvolution<DIM> bubble(lat, parallel, "test");
		bubble.RecordParameters();
		bubble.SaveParameters("test_param.txt");
		ElectroweakObserver<DIM> obs(bubble);
		for (int t = 0; t < 5; ++t) {
			obs.Measure();
			bubble.UpdateFields();
			bubble.EvolveInterior_KS();
			bubble.TimeAdvance();
		}
		obs.SaveDataTable("test_datatable.txt");*/
		
		//Parallel2D parallel(MPI_COMM_WORLD);
		//parallel.InitializeGrid(n_rows, n_cols);
		/*
		Parallel2D parallel(n_rows, n_cols, MPI_COMM_WORLD);
		GridIndexType node_size[] = { n_rows, n_cols };
		GridIndexType grid_rank[] = { parallel.get_grid_rank()[0], parallel.get_grid_rank()[1] };

		GridIndexType grid_loc[2];
		transform_gridrank_to_gridloc(grid_rank, grid_loc);
		IndexType lat_size[3] = {12,24,36};
		Lattice<3> lat(lat_size, 2, node_size, grid_loc);
		Site<3> x(lat);
		Field<int, 3> f(lat, 1, parallel);
		for(x.first(); x.test(); x.next()){
			f(x, 0) = x.get_index();
		}
		f.update_halo();
		f.write("f.h5");*/
		SU2vector a;
		std::cout<<sizeof(a)<<std::endl;

		

		/*
		HDF5Wrapper h5("test1.h5");
		int data[100];
		std::generate_n(data, 100, [n = 0]() mutable {return static_cast<int>(n++);});
		hsize_t size[2] = {6,6};
		hsize_t mem_size[2] = {10,10};
		for(auto b = 0; b < 3; ++b){
			hsize_t file_offset[2] = {static_cast<hsize_t>(2*b),0};
			hsize_t mem_offset[2] = {static_cast<hsize_t>(2*b+1),3};
			hsize_t block_size[2] = {2,6};
			h5.SaveSingleDatasetFile<int,2>(
				"test.h5", 
				"int_data", 
				size, 
				file_offset, 
				block_size,
				mem_size,
				mem_offset,
				block_size, 
				data, b == 0);
		}
		*/
	}
	Parallel_Finalize();
	return 0;
} 