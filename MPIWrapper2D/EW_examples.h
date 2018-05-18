#ifndef EW_EXAMPLES_H
#define EW_EXAMPLES_H

#include "EW_Base.h"

    using namespace Electroweak;

void EW_BaseModel_SymmetricInit(const int n_rows, const int n_cols){

    Parallel2D parallel(n_rows, n_cols, MPI_COMM_WORLD);
    GridIndexType node_size[] = { n_rows, n_cols };
    GridIndexType grid_rank[] = { parallel.get_grid_rank()[0], parallel.get_grid_rank()[1] };

    GridIndexType grid_loc[2];
	transform_gridrank_to_gridloc(grid_rank, grid_loc);

    Lattice<DIM> lat(nSize, halo, node_size, grid_loc);

    ElectroweakEvolution<DIM> bubble(lat, parallel, "base_model");
    bubble.RecordParameters();
    bubble.SaveParameters("base_model_param.txt");
    bubble.TrivialInitialCondition();

    ElectroweakObserver<DIM> obs(bubble);
    obs.SetObservables(ObserverFlags::OBS_EnergyAllParts |
                        ObserverFlags::OBS_CSNumber |
                        ObserverFlags::OBS_HiggsMagnitude2,

                        ObserverFlags::OBS_TotalEnergy | 
                        ObserverFlags::OBS_CSNumber |
                        ObserverFlags::OBS_HiggsMagnitude2);
    
    //Evolve 10 steps.
    for(auto i = 0; i < Ntimesteps; ++i){
        obs.Measure();
        bubble.UpdateFields();
        bubble.EvolveInterior_KS();
        bubble.TimeAdvance();
        obs.SaveDensityData("base_model_den_" + std::to_string(i) + ".h5");
    }
    obs.SaveDataTable("base_model_dtable.txt");
    return;
}



#endif