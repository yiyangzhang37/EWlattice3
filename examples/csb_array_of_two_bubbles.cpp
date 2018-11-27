#include <iostream>

#include "./ParaSite/ParaSite.h"

#include "./EW_Model/EW_examples.h"

using namespace MPI_Wrapper;
using namespace ParaSite;
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

		using namespace EW_BubbleNucleation;
        Parallel2D parallel(n_rows, n_cols, MPI_COMM_WORLD);
        GridIndexType node_size[] = { n_rows, n_cols };
        GridIndexType grid_rank[] = { parallel.get_grid_rank()[0], parallel.get_grid_rank()[1] };

        GridIndexType grid_loc[2];
        transform_gridrank_to_gridloc(grid_rank, grid_loc);

        Lattice<DIM> lat(nSize, halo, node_size, grid_loc);
        std::string id = ReadConfigIni("RUN_ID");
        CSBubble<DIM> bubble(lat, parallel, id);
        bubble.RecordParameters();
        bubble.SaveParameters(id+"_param.txt");
        
        NucleationObserver<DIM> obs(bubble);
        obs.SetObservables(ObserverFlags::OBS_EnergyAllParts |
                            ObserverFlags::OBS_CSNumber |
                            ObserverFlags::OBS_MagneticEnergy |
                            ObserverFlags::OBS_HiggsMagnitude2 |
                            ObserverFlags::OBS_MinHiggsMagnitude2 |
                            ObserverFlags::OBS_HiggsWinding |
                            ObserverFlags::OBS_NewBubbleCount,

                            ObserverFlags::OBS_TotalEnergy | 
                            ObserverFlags::OBS_MagneticEnergy |
                            ObserverFlags::OBS_HiggsMagnitude2);
        
        bubble.InitializeSymmetricPhase();

        SU2vector phi1_hat, phi2_hat;
        phi1_hat(0) = Cmplx(0, 0);
        phi1_hat(1) = Cmplx(1, 0);

        //phi2_hat(0) = Cmplx(1, 0)/sqrt(2.0);
        //phi2_hat(1) = Cmplx(0, -1)/sqrt(2.0);
        phi2_hat = phi1_hat;

        for(auto i = 0; i <= Ntimesteps; ++i){
            if(DensityDataCalcFreq == 0 || bubble.get_time_step() % DensityDataCalcFreq == 0){
                obs.Measure();
            }

            bubble.UpdateFields();

            bubble.ArrayOfTwoBubblesTest_WithWinding(1, 1, phi1_hat, phi2_hat);
        
            bubble.EvolveInterior_RadialDamping();
            obs.SaveDensityData(id + "_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);
            obs.SaveDataTable(id+"_dtable.txt", 50);

            bubble.TimeAdvance();
        }
        obs.SaveDataTable(id+"_dtable.txt");
        bubble.ConcludeEvolution(id+"_param.txt");
		
	}
	Parallel_Finalize();
	return 0;
}