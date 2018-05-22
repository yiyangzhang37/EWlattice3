#include "EW_examples.h"


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


void EW_Nucl_OneBubble(const int n_rows, const int n_cols){
    using namespace EW_BubbleNucleation;
    Parallel2D parallel(n_rows, n_cols, MPI_COMM_WORLD);
    GridIndexType node_size[] = { n_rows, n_cols };
    GridIndexType grid_rank[] = { parallel.get_grid_rank()[0], parallel.get_grid_rank()[1] };

    GridIndexType grid_loc[2];
	transform_gridrank_to_gridloc(grid_rank, grid_loc);

    Lattice<DIM> lat(nSize, halo, node_size, grid_loc);
    BubbleNucleation<DIM> bubble(lat, parallel, "one_bubble");
    bubble.RecordParameters();
    bubble.RecordCustomParameters();
    bubble.SaveParameters("one_bubble_param.txt");
    
    NucleationObserver<DIM> obs(bubble);
    obs.SetObservables(ObserverFlags::OBS_EnergyAllParts |
                        ObserverFlags::OBS_MagneticEnergy |
                        ObserverFlags::OBS_HiggsMagnitude2 |
                        ObserverFlags::OBS_MinHiggsMagnitude2 |
                        ObserverFlags::OBS_NewBubbleCount,

                        ObserverFlags::OBS_TotalEnergy | 
                        ObserverFlags::OBS_MagneticEnergy |
                        ObserverFlags::OBS_HiggsMagnitude2);


    bubble.InitializeSymmetricPhase();

    for(auto i = 0; i <= Ntimesteps; ++i){
        obs.ExtendMeasure();
        bubble.UpdateFields();

        /*nucleation is inserted here*/
        bubble.OneBubbleTest();

        bubble.EvolveInterior_KS();
        obs.SaveDensityData("one_bubble_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);

        bubble.TimeAdvance();
    }
    obs.SaveDataTable("one_bubble_dtable.txt");
    bubble.ConcludeEvolution("one_bubble_param.txt");
    return;
}


void EW_Nucl_TwoBubbles(const int n_rows, const int n_cols){
    using namespace EW_BubbleNucleation;
    Parallel2D parallel(n_rows, n_cols, MPI_COMM_WORLD);
    GridIndexType node_size[] = { n_rows, n_cols };
    GridIndexType grid_rank[] = { parallel.get_grid_rank()[0], parallel.get_grid_rank()[1] };

    GridIndexType grid_loc[2];
	transform_gridrank_to_gridloc(grid_rank, grid_loc);

    Lattice<DIM> lat(nSize, halo, node_size, grid_loc);
    std::string id = "two_bubbles";
    BubbleNucleation<DIM> bubble(lat, parallel, id);
    bubble.RecordParameters();
    bubble.RecordCustomParameters();
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


    bubble.InitializeSymmetricPhase();

    for(auto i = 0; i <= Ntimesteps; ++i){
        obs.ExtendMeasure();
        bubble.UpdateFields();

        /*nucleation is inserted here*/
        bubble.TwoBubblesTest(BUBBLES_HALF_SEP);

        bubble.EvolveInterior_KS();
        obs.SaveDensityData(id+"_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);

        bubble.TimeAdvance();
    }
    obs.SaveDataTable(id+"_dtable.txt");
    bubble.ConcludeEvolution(id+"_param.txt");
    return;
}

void EW_Random_Nucl(const int n_rows, const int n_cols){
    using namespace EW_BubbleNucleation;
    Parallel2D parallel(n_rows, n_cols, MPI_COMM_WORLD);
    GridIndexType node_size[] = { n_rows, n_cols };
    GridIndexType grid_rank[] = { parallel.get_grid_rank()[0], parallel.get_grid_rank()[1] };

    GridIndexType grid_loc[2];
	transform_gridrank_to_gridloc(grid_rank, grid_loc);

    Lattice<DIM> lat(nSize, halo, node_size, grid_loc);
    std::string id = "random_nucl";
    BubbleNucleation<DIM> bubble(lat, parallel, id);
    bubble.RecordParameters();
    bubble.RecordCustomParameters();
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


    bubble.InitializeSymmetricPhase();

    for(auto i = 0; i <= Ntimesteps; ++i){
        obs.ExtendMeasure();
        bubble.UpdateFields();

        /*nucleation is inserted here*/
        bubble.RandomBubbleNucleation();
       
        bubble.EvolveInterior_RadialDamping();
        obs.SaveDensityData(id+"_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);

        bubble.TimeAdvance();
    }
    obs.SaveDataTable(id+"_dtable.txt");
    bubble.ConcludeEvolution(id+"_param.txt");
    return;
}