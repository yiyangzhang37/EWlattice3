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
        BubbleNucleation<DIM> bubble(lat, parallel, id);
        bubble.RecordParameters();
        bubble.SaveParameters(id+"_param.txt");
        
        NucleationObserver<DIM> obs(bubble);
        obs.SetObservables(ObserverFlags::OBS_EnergyAllParts |
                            ObserverFlags::OBS_MagneticEnergy |
                            ObserverFlags::OBS_HiggsMagnitude2 |
                            ObserverFlags::OBS_MinHiggsMagnitude2 |
                            ObserverFlags::OBS_NewBubbleCount,

                            ObserverFlags::OBS_TotalEnergy | 
                            ObserverFlags::OBS_MagneticEnergy |
                            ObserverFlags::OBS_HiggsMagnitude2);
    
        NucleationObserver<DIM> bfield(bubble);
        bfield.SetObservables(ObserverFlags::OBS_MagneticField, 
                                ObserverFlags::OBS_MagneticField);

        bubble.InitializeSymmetricPhase();

        for(auto i = 0; i <= Ntimesteps; ++i){
            obs.Measure();

            bubble.UpdateFields();

            /*nucleation is inserted here*/
            if( ! bubble.isHiggsAllInBrokenPhase(obs) ){
                bubble.RandomBubbleNucleation();
            }
        
            bubble.EvolveInterior_RadialDamping();
            obs.SaveDensityData(id + "_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);
            obs.SaveDataTable(id+"_dtable.txt", 50);

            /*early stop*/
            if( bubble.CheckEarlyStop(obs) ) {
                obs.Measure();
                obs.SaveDensityData(id + "_den_es" + ".h5");
                bfield.Measure();
                bfield.SaveDensityData(id + "_bfield.h5");
                break;
            }
            bubble.TimeAdvance();
        }
        obs.SaveDataTable(id+"_dtable.txt");
        bubble.ConcludeEvolution(id+"_param.txt");
		
	}
	Parallel_Finalize();
	return 0;
}