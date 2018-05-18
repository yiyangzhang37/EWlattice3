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

#include "EW_examples.h"

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

		EW_BaseModel_SymmetricInit(n_rows, n_cols);

	}
	Parallel_Finalize();
	return 0;
} 