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
        obs.Measure();
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
    //lat.make_all_index_tables();
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


    bubble.InitializeSymmetricPhase();

    for(auto i = 0; i <= Ntimesteps; ++i){
        obs.Measure();
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
        bubble.NonRandomTest(BUBBLES_HALF_SEP);
       
        bubble.EvolveInterior_RadialDamping();
        obs.SaveDensityData(id + "_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);
	    obs.SaveDataTable(id + "_dtable.txt", 50);

        if( (i - RECORD_TIME_OFFSET >= 0) && (i - RECORD_TIME_OFFSET) % 50 == 0 ){
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
    return;
}

void EW_Random_Nucl_LongTimeSpectrum(const int n_rows, const int n_cols){
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
        if(DensityDataCalcFreq == 0 || bubble.get_time_step() % DensityDataCalcFreq == 0){
            obs.Measure();
        }

        bubble.UpdateFields();

        /*nucleation is inserted here*/
        if( ! bubble.isHiggsAllInBrokenPhase(obs) ){
            bubble.RandomBubbleNucleation();
        }
       
        bubble.EvolveInterior_RadialDamping();
        obs.SaveDensityData(id + "_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);
	    obs.SaveDataTable(id+"_dtable.txt", 50);

        /*save bfield*/
        if( BFIELD_SAVE_FREQ != 0 && bubble.get_time_step() % BFIELD_SAVE_FREQ == 0 ) {
            bfield.Measure();
            bfield.SaveDensityData(id + "_bfield_" + std::to_string(bubble.get_time_step()) + ".h5");
        }

        if(bubble.get_time_step() % 10000 == 0 && bubble.get_time_step() > 0 ){
            bubble.SaveFields(std::to_string(bubble.get_time_step()) + ".h5");
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
        //if( ! bubble.CheckHiggsAllInBrokenPhase(obs) ){
            bubble.RandomBubbleNucleation();
        //}
       
        bubble.EvolveInterior_RadialDamping();
        obs.SaveDensityData(id + "_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);
	    obs.SaveDataTable(id+"_dtable.txt", 50);

        if(bubble.get_time_step() % 20 == 0){
            bfield.Measure();
            //const auto& bfield_data = bfield.get_density_data();
        }

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
    return;
}

void EW_CSB_OneBubble(const int n_rows, const int n_cols, const std::string& method){
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
                        ObserverFlags::OBS_GaussConstraint |
                        ObserverFlags::OBS_HiggsWinding |
                        ObserverFlags::OBS_NewBubbleCount,

                        ObserverFlags::OBS_TotalEnergy |
                        ObserverFlags::OBS_CSNumber | 
                        ObserverFlags::OBS_MagneticEnergy |
                        ObserverFlags::OBS_HiggsMagnitude2);
    
    bubble.InitializeSymmetricPhase();
    //bubble.InitPureGauge(3, DX*10);

    SU2vector phi_hat;
    phi_hat(0) = Cmplx(0, 0);
    phi_hat(1) = Cmplx(1, 0);

    for(auto i = 0; i <= Ntimesteps; ++i){
        obs.Measure();

        bubble.UpdateFields();
        if(method == "unperturbed"){
            bubble.OneBubbleTest_WithWinding(0, phi_hat);
        } else if(method == "perturbed"){
            bubble.OneBubbleTest_WithWinding_Perturbed(0, phi_hat);
        } else {
            std::cout << "method error." << std::endl;
        }
       
        bubble.EvolveInterior_RadialDamping();
        obs.SaveDensityData(id + "_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);
	    obs.SaveDataTable(id+"_dtable.txt", 5);

        bubble.TimeAdvance();
    }
    obs.SaveDataTable(id+"_dtable.txt");
    bubble.ConcludeEvolution(id+"_param.txt");
    return;
}

void EW_CSB_TwoBubbles(const int n_rows, const int n_cols, const std::string& method){
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
                        ObserverFlags::OBS_GaussConstraint |
                        ObserverFlags::OBS_HiggsWinding |
                        ObserverFlags::OBS_NewBubbleCount,

                        ObserverFlags::OBS_TotalEnergy | 
                        ObserverFlags::OBS_MagneticEnergy |
                        ObserverFlags::OBS_HiggsMagnitude2);
    
    bubble.InitializeSymmetricPhase();

    SU2vector phi1_hat, phi2_hat;
    phi1_hat(0) = Cmplx(0, 0);
    phi1_hat(1) = Cmplx(1, 0);

    phi2_hat(0) = Cmplx(0, 0);
    phi2_hat(1) = Cmplx(1, 0);

    for(auto i = 0; i <= Ntimesteps; ++i){
        if(DensityDataCalcFreq == 0 || bubble.get_time_step() % DensityDataCalcFreq == 0){
            obs.Measure();
        }

        bubble.UpdateFields();

        if(method == "unperturbed"){
            bubble.TwoBubblesTest_WithWinding(0, 0, phi1_hat, phi2_hat);
        } else if(method == "perturbed"){
            bubble.TwoBubblesTest_WithWinding_Perturbed(0, 0, phi1_hat, phi2_hat);
        } else{
            std::cout << "method error." << std::endl;
        }
        bubble.EvolveInterior_RadialDamping();
        obs.SaveDensityData(id + "_den_" + std::to_string(i) + ".h5", DensityDataSaveFreq);
	    obs.SaveDataTable(id+"_dtable.txt", 50);
        //if(bubble.get_time_step() % 100 == 0){
        //    bubble.get_phi().write("phi_"+std::to_string(bubble.get_time_step())+".dat");
        //}
        bubble.TimeAdvance();
    }
    obs.SaveDataTable(id+"_dtable.txt");
    bubble.ConcludeEvolution(id+"_param.txt");
    return;
}

void EW_CSB_ArrayOfTwoBubbles(const int n_rows, const int n_cols){
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
    return;
}