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
    std::string id = ReadConfigIni("RUN_ID");
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
        bubble.OneBubbleTest();

        bubble.EvolveInterior_KS();
        obs.SaveDensityData(id + "_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);

        bubble.TimeAdvance();
    }
    obs.SaveDataTable(id + "_dtable.txt");
    bubble.ConcludeEvolution(id + "_param.txt");
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
    std::string id = ReadConfigIni("RUN_ID");
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

void EW_Nucl_NonRand(const int n_rows, const int n_cols) {
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
 
    NucleationObserver<DIM> bfield(bubble);
    bfield.SetObservables(ObserverFlags::OBS_MagneticField, 
                            ObserverFlags::OBS_MagneticField);

    bubble.InitializeSymmetricPhase();

    for(auto i = 0; i <= Ntimesteps; ++i){
        obs.ExtendMeasure();

        bubble.UpdateFields();
        bubble.NonRandomTest(BUBBLES_HALF_SEP);
       
        bubble.EvolveInterior_RadialDamping();
        obs.SaveDensityData(id + "_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);
	    obs.SaveDataTable(id + "_dtable.txt", 50);

        if( (i - RECORD_TIME_OFFSET >= 0) && (i - RECORD_TIME_OFFSET) % 20 == 0 ){
            bfield.Measure();
            bfield.SaveDensityData(id + "_bfield_" + std::to_string(i) + ".h5");
        }

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
    std::string id = ReadConfigIni("RUN_ID");
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
 
    NucleationObserver<DIM> bfield(bubble);
    bfield.SetObservables(ObserverFlags::OBS_MagneticField, 
                            ObserverFlags::OBS_MagneticField);

    bubble.InitializeSymmetricPhase();

    for(auto i = 0; i <= Ntimesteps; ++i){
        obs.ExtendMeasure();

        bubble.UpdateFields();

        /*nucleation is inserted here*/
        //if( ! bubble.CheckHiggsAllInBrokenPhase(obs) ){
            bubble.RandomBubbleNucleation();
        //}
       
        bubble.EvolveInterior_RadialDamping();
        obs.SaveDensityData(id + "_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);
	    obs.SaveDataTable(id+"_dtable.txt", 50);

        /*early stop*/
        if( bubble.CheckEarlyStop(obs) ) {
            obs.ExtendMeasure();
            obs.SaveDensityData(id + "_den_es" + ".h5");
            bfield.Measure();
            bfield.SaveDensityData(id + "_bfield.h5");
            break;
        }

        bubble.TimeAdvance();
    }
    obs.SaveDataTable(id+"_dtable.txt");
    bubble.ConcludeEvolution(id+"_param.txt");
    return;
}

void EW_Random_Nucl_FFT(const int n_rows, const int n_cols){
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
 
    NucleationObserver<DIM> bfield(bubble);
    bfield.SetObservables(ObserverFlags::OBS_MagneticField, 
                            ObserverFlags::OBS_MagneticField);

    bubble.InitializeSymmetricPhase();

    for(auto i = 0; i <= Ntimesteps; ++i){
        obs.ExtendMeasure();

        bubble.UpdateFields();

        /*nucleation is inserted here*/
        //if( ! bubble.CheckHiggsAllInBrokenPhase(obs) ){
            bubble.RandomBubbleNucleation();
        //}
       
        bubble.EvolveInterior_RadialDamping();
        obs.SaveDensityData(id + "_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);
	    obs.SaveDataTable(id+"_dtable.txt", 50);

        if(bubble.get_time_step() % 20 == 0){
            bfield.Measure();
            const auto& bfield_data = bfield.get_density_data();
        }

        /*early stop*/
        if( bubble.CheckEarlyStop(obs) ) {
            obs.ExtendMeasure();
            obs.SaveDensityData(id + "_den_es" + ".h5");
            bfield.Measure();
            bfield.SaveDensityData(id + "_bfield.h5");
            break;
        }

        bubble.TimeAdvance();
    }
    obs.SaveDataTable(id+"_dtable.txt");
    bubble.ConcludeEvolution(id+"_param.txt");
    return;
}