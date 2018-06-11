//#include "./tests/Test_MPIWrapper2D.h"
#include <iostream>
#include <chrono>

#include "./ParaSite/ParaSite.h"

#include "./EW_Model/EW_examples.h"

#include "./tests/test_fft.h"

using namespace MPI_Wrapper;
using namespace ParaSite;
//using namespace Electroweak;
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
		
		EW_Nucl_NonRand(n_rows, n_cols);
		//EW_Random_Nucl(n_rows, n_cols);
		//test_FFT3D_r2cc2r();
	}
	Parallel_Finalize();
	return 0;
} 