#ifndef EW_OBSERVATION_H
#define EW_OBSERVATION_H

#include "EW_Model.h"
#include <algorithm>
#include <functional>

namespace Electroweak{

    using namespace ParaSite;

	typedef unsigned int FlagType;

	namespace ObserverFlags {

		FlagType OBS_TotalEnergy = 1 << 0;      	// TotalEnergy
		FlagType OBS_CSNumber = 1 << 1;         	// CSNumber
		FlagType OBS_MagneticEnergy = 1 << 2;   	// MagneticEnergy
		FlagType OBS_MagneticHelicity = 1 << 3;     // MagneticHelicity
		FlagType OBS_ElectricEnergy = 1 << 4;       // ElectricEnergy
		FlagType OBS_GaussConstraint = 1 << 5;      // GaussConstraint
		FlagType OBS_HiggsMagnitude2 = 1 << 6;      // HiggsMagnitude2
		FlagType OBS_MinHiggsMagnitude2 = 1 << 7;   // MinHiggsMagnitude2
		FlagType OBS_HiggsWinding = 1 << 8;         // HiggsWinding
		FlagType OBS_EnergyKinetic = 1 << 9;
		FlagType OBS_EnergyGradient = 1 << 10;
		FlagType OBS_EnergyPotential = 1 << 11;
		FlagType OBS_EnergyU1 = 1 << 12;
		FlagType OBS_EnergySU2 = 1 << 13;

		//Combined observations
		FlagType OBS_EnergyAllParts = OBS_TotalEnergy | OBS_EnergyKinetic | OBS_EnergyGradient
						| OBS_EnergyPotential | OBS_EnergyU1 | OBS_EnergySU2;
		
		//The flags can be extended.
		//Custom flags should be start at 1 << 21.
	}
	
	template<int DIM>
    class ElectroweakObserver{
    public:
        ElectroweakObserver(const ElectroweakEvolution<DIM>& evo);
        ~ElectroweakObserver() = default;
		
		/*
		set flags for data_table observations and density_data observations.
		density_data_flags should be a subset of data_table_flags.
		This function is declared as virtual and it should be re-implemented in derived class.
		*/
		virtual void SetObservables(
			const FlagType data_table_flags,
			const FlagType density_data_flags);
		/*
		measure all the quantities determined by the flags.
		*/
		void Measure();
		/*
		Measure more quantities, including those defined in the derived class.
		This function should be re-implemented in the derived class.
		The form should be:
		void ExtendMeasure(){
			this->Measure();
			auto time_step = this->evo_.get_time_step();
			if (this->data_table_flags_ & OBS_YourCustomFlag) {
				this->CalcCustomQuantity(time_step);
			}
			...
			return;
		}
		*/
		virtual void ExtendMeasure() {;} 

		void CalcEnergy(const int time_step);
		void CalcCSNumber(const int time_step);
		void CalcMagneticEnergy(const int time_step);
		void CalcMagneticHelicity(const int time_step);
		void CalcElectricEnergy(const int time_step);
		//CalcGaussConstraint computes the Hamiltonian of Gauss constraints.
		void CalcGaussConstraint(const int time_step);
		//CalcHiggsMagnitude2 outputs the average HiggsMagnitude^2 to Data Table.
		void CalcHiggsMagnitude2(const int time_step);
		void CalcMinHiggsMagnitude2(const int time_step);
		void CalcHiggsWinding(const int time_step);

		void CalcSimple(const int time_step, 
			const std::string& data_table_key,
			const std::string& density_key,
			const std::function<Real(const Site<DIM>&, const int)> sc_quantity);

		/*
		Save the data table.
		save operation will be performed on the root process.
		*/
		void SaveDataTable(const std::string& file_name, 
						const int freq = 1,
						const int starting_time_step = 0) const {
			if(this->evo_.get_parallel_object().is_root()){
				auto t = this->evo_.get_time_step();
				if( ((t - starting_time_step) >= 0) && ((t - starting_time_step)%freq ==0) ){
					this->data_table_.save_table(file_name);
				}
			}
		}

		/*
		save the density data for a specific time step.
		Parallel save procedure has been taken care of in the Field class.
		*/
		void SaveDensityData(const std::string& file_name, 
							const int freq = 1, 
							const int starting_time_step = 0) const;

    protected:
		const ElectroweakEvolution<DIM>& evo_;
		const Field< SU2vector, DIM >& phi_;
		const Field< SU2vector, DIM >& pi_;
		const Field< SU2matrix, DIM >& U_;
		const Field< SU2matrix, DIM >& F_;
		const Field< U1matrix, DIM >& V_;
		const Field< U1matrix, DIM>& E_;
		const Parallel2D& parallel_;

		FlagType data_table_flags_;
		FlagType density_data_flags_;

		/*total quantities*/
		std::vector<std::string> data_table_names_;
        DataTable data_table_;

		/*density quantities*/
		std::vector<std::string> density_names_;
        Field< RealSave, DIM > density_data_;

		/*
		The following two functions are declared as virtual,
		and should be re-implemented in derived class,
		with init_name_vector() replaced by init_extend_name_vector().
		*/
		virtual void init_data_table();
		virtual void init_density_data();

		void init_name_vector(
			const FlagType flags,
			std::vector<std::string>& names) const;
		/*
		This function should be re-implemented in derived class to take the position of init_name_vector.
		The function should take the form:
		init_extend_name_vector(...){
			this->init_name_vector(...);
			if (flags & ObserverFlags::YourCustomFlag) {
				names.push_back("YourCustomId");
			}
			...
			return;
		}
		*/
		virtual void init_extend_name_vector(
					const FlagType flags,
					std::vector<std::string>& names) const {;}

		int find_index(const std::vector<std::string>& names, const std::string& key) const;

		/*single cell quantities*/
		/*
		Normal single cell quantities should accept two parameters, in the form:
		Real single_cell_func_(const Site<DIM>& x, const int T) const;
		The return value of the function does not include the DX^3 factor.
		*/
		Real sc_GaussU1(const Site<DIM>& x, const int T) const;
		Real sc_GaussSU2(const Site<DIM>& x, const int T,
			const int su2dir) const;
		Real sc_HamiltonianGC(const Site<DIM>& x, const int T) const;
		Real sc_MagneticEnergy(const Site<DIM>& x, const int T) const;
		Real sc_ElectricEnergy(const Site<DIM>& x, const int T) const;
		Real sc_HiggsWinding(const Site<DIM>& x, const int T) const;

    };


	template<int DIM>
	ElectroweakObserver<DIM>::ElectroweakObserver(const ElectroweakEvolution<DIM>& evo)
		:
		evo_(evo),
		phi_(evo.get_phi()),
		pi_(evo.get_pi()),
		U_(evo.get_U()),
		F_(evo.get_F()),
		V_(evo.get_V()),
		E_(evo.get_E()),
		parallel_(evo.get_parallel_object()),
		density_data_(evo.get_lattice(), 1, 1, parallel_)
	{
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::SetObservables(
		const FlagType data_table_flags,
		const FlagType density_data_flags) {
		this->data_table_flags_ = data_table_flags;
		this->density_data_flags_ = density_data_flags;
		init_data_table();
		init_density_data();
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::init_data_table() {
		std::vector<std::string> data_table_names;
		this->init_name_vector(this->data_table_flags_, data_table_names);
		this->data_table_.initialize(data_table_names);
		this->data_table_names_ = data_table_names;
		return;
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::init_density_data() {
		std::vector<std::string> density_data_names;
		this->init_name_vector(this->density_data_flags_, density_data_names);
		this->density_names_ = density_data_names;
		this->density_data_.reinit(this->density_names_.size(), 1);
		return;
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::init_name_vector(
		const FlagType flags,
		std::vector<std::string>& names) const {
		names.clear();
		if (flags & ObserverFlags::OBS_TotalEnergy) {
			names.push_back("TotalEnergy");
		}
		if (flags & ObserverFlags::OBS_CSNumber) {
			names.push_back("CSNumber");
		}
		if (flags & ObserverFlags::OBS_MagneticEnergy) {
			names.push_back("MagneticEnergy");
		}
		if (flags & ObserverFlags::OBS_MagneticHelicity) {
			names.push_back("MagneticHelicity");
		}
		if (flags & ObserverFlags::OBS_ElectricEnergy) {
			names.push_back("ElectricEnergy");
		}
		if (flags & ObserverFlags::OBS_GaussConstraint) {
			names.push_back("GaussConstraint");
		}
		if (flags & ObserverFlags::OBS_HiggsMagnitude2) {
			names.push_back("HiggsMagnitude2");
		}
		if (flags & ObserverFlags::OBS_MinHiggsMagnitude2){
			names.push_back("MinHiggsMagnitude2");
		}
		if (flags & ObserverFlags::OBS_HiggsWinding) {
			names.push_back("HiggsWinding");
		}
		if (flags & ObserverFlags::OBS_EnergyKinetic) {
			names.push_back("E_Kinetic");
		}
		if (flags & ObserverFlags::OBS_EnergyGradient) {
			names.push_back("E_Grad");
		}
		if (flags & ObserverFlags::OBS_EnergyPotential) {
			names.push_back("E_Pot");
		}
		if (flags & ObserverFlags::OBS_EnergyU1) {
			names.push_back("E_U1");
		}
		if (flags & ObserverFlags::OBS_EnergySU2) {
			names.push_back("E_SU2");
		}
		return;
	}

	template<int DIM>
	int ElectroweakObserver<DIM>::find_index(
		const std::vector<std::string>& names,
		const std::string& key) const {
		auto loc = std::find(names.cbegin(), names.cend(), key);
		if (loc == names.end()) return -1;
		else return loc - names.begin();
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::CalcEnergy(const int time_step) {
		//in this function, the energy is calculated as the current step, not the updated step.
		//The problem of this function is, when initial condition is not accurate, the calculated energy density may be not positive-definite.
		auto nowTime = time_step % CYCLE;
		auto nextTime = (time_step + 1) % CYCLE;

		Site<DIM> x(this->evo_.get_lattice());

		auto dt_idx_total_energy = find_index(this->data_table_names_, "TotalEnergy");
		auto dt_idx_kinetic = find_index(this->data_table_names_, "E_Kinetic");
		auto dt_idx_grad = find_index(this->data_table_names_, "E_Grad");
		auto dt_idx_pot = find_index(this->data_table_names_, "E_Pot");
		auto dt_idx_u1 = find_index(this->data_table_names_, "E_U1");
		auto dt_idx_su2 = find_index(this->data_table_names_, "E_SU2");
		auto den_idx_total_energy = find_index(this->density_names_, "TotalEnergy");
		
		auto total_energy = 0.0;
		auto total_e_kinetic = 0.0;
		auto total_e_grad = 0.0;
		auto total_e_pot = 0.0;
		auto total_e_u1 = 0.0;
		auto total_e_su2 = 0.0;

		for (x.first(); x.test(); x.next()) {
			if (this->evo_.is_boundary(x)) {
				if (den_idx_total_energy != -1) {
					this->density_data_(x, den_idx_total_energy) = 0.0;
				}
			}
			else {
				Real eKinetic = 0, eGrad = 0, ePotential = 0, eU1 = 0;
				Cmplx eSU2(0, 0);
				eKinetic = (0.5*(pi_(x, nowTime) + pi_(x, nextTime))).squaredNorm();
				ePotential = Vphi(x, phi_, nowTime);
				
				for (int i = 0; i<DIM; ++i) {
					eGrad += 1.0 / DX2*(U_(x, i, nowTime)*V_(x, i, nowTime)*phi_(x + i, nowTime) - phi_(x, nowTime)).squaredNorm();
					eU1 += 4.0 / (pow(DT*DX*gp, 2)) * (1.0 - 1.0 / UNITY_E * 0.5*(E_(x, i, nowTime) + E_(x, i, nextTime)).real());
					eSU2 += 2.0 / (pow(DT*DX*g, 2)) * (2.0 - 1.0 / UNITY_F * 0.5*(F_(x, i, nowTime) + F_(x, i, nextTime)).trace());
					for (int j = 0; j < DIM; ++j) {
						eU1 += 2.0 / (DX4*pow(gp, 2)) * (1.0 - Plaquette(x, V_, i, j, nowTime).real());
						eSU2 += 1.0 / (DX4*pow(g, 2)) * (2.0 - Plaquette(x, U_, i, j, nowTime).trace());
					}
				}
				if (den_idx_total_energy != -1) {
					this->density_data_(x, den_idx_total_energy) = 
						static_cast<RealSave> ((eKinetic + eGrad + ePotential + eU1 + eSU2.real()));
				}
				total_energy += (eKinetic + ePotential + eGrad + eU1 + eSU2.real())*DX3;
				total_e_kinetic += eKinetic * DX3;
				total_e_grad += eGrad * DX3;
				total_e_pot += ePotential * DX3;
				total_e_u1 += eU1 * DX3;
				total_e_su2 += eSU2.real() * DX3;

			}
		}
		this->parallel_.Sum(total_energy);
		this->parallel_.Sum(total_e_kinetic);
		this->parallel_.Sum(total_e_grad);
		this->parallel_.Sum(total_e_pot);
		this->parallel_.Sum(total_e_u1);
		this->parallel_.Sum(total_e_su2);
		if (dt_idx_total_energy != -1) this->data_table_.append_value(dt_idx_total_energy, total_energy);
		if (dt_idx_kinetic != -1) this->data_table_.append_value(dt_idx_kinetic, total_e_kinetic);
		if (dt_idx_grad != -1) this->data_table_.append_value(dt_idx_grad, total_e_grad);
		if (dt_idx_pot != -1) this->data_table_.append_value(dt_idx_pot, total_e_pot);
		if (dt_idx_u1 != -1) this->data_table_.append_value(dt_idx_u1, total_e_u1);
		if (dt_idx_su2 != -1) this->data_table_.append_value(dt_idx_su2, total_e_su2);
		return;
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::CalcCSNumber(const int time_step){
		auto updateTime = time_step % CYCLE;
		Site<DIM> x(this->evo_.get_lattice());

		auto dt_idx_cs = find_index(this->data_table_names_, "CSNumber");
		auto den_idx_cs = find_index(this->density_names_, "CSNumber");

		Real totalCSnumber = 0;

		const auto pre_factor = FermiFamily / (32.0*PI*PI) * DX3;
		const auto g3 = 2.0*g*g*g / 3.0;

		for (x.first(); x.test(); x.next()) {
			Real tempSave = 0.0;
			if (this->evo_.is_boundary(x)) {
				if (den_idx_cs != -1) this->density_data_(x, den_idx_cs) = 0.0;
				continue;
			}
			else {
				const Real SU2W[DIM][SU2DIM] = {
					{ SU2_W_SITE(x, U_, X_AXIS, 0, updateTime), SU2_W_SITE(x, U_, X_AXIS, 1, updateTime), SU2_W_SITE(x, U_, X_AXIS, 2, updateTime) },
					{ SU2_W_SITE(x, U_, Y_AXIS, 0, updateTime), SU2_W_SITE(x, U_, Y_AXIS, 1, updateTime), SU2_W_SITE(x, U_, Y_AXIS, 2, updateTime) },
					{ SU2_W_SITE(x, U_, Z_AXIS, 0, updateTime), SU2_W_SITE(x, U_, Z_AXIS, 1, updateTime), SU2_W_SITE(x, U_, Z_AXIS, 2, updateTime) } };
				const Real U1Y[DIM] = {
					U1_Y_SITE(x, V_, X_AXIS, updateTime),
					U1_Y_SITE(x, V_, Y_AXIS, updateTime),
					U1_Y_SITE(x, V_, Z_AXIS, updateTime)
				};
				Real U1B[DIM][DIM];
				Real SU2B[DIM][DIM][SU2DIM];
				for (int i = 0; i < DIM; ++i) {
					for (int j = i; j < DIM; ++j) {
						U1B[i][j] = U1_B_SITE(x, V_, i, j, updateTime);
						U1B[j][i] = -U1B[i][j];
						for (int a = 0; a < SU2DIM; ++a) {
							SU2B[i][j][a] = SU2_B_SITE(x, U_, i, j, a, updateTime);
							SU2B[j][i][a] = -SU2B[i][j][a];
						}
					}
				}
				for (int i = 0; i<DIM; ++i) {
					int j = (i + 1) % DIM;
					int k = (i + 2) % DIM;
					tempSave -= gp*gp*(U1B[i][j] * U1Y[k] - U1B[i][k] * U1Y[j]);
					for (int a = 0; a<SU2DIM; ++a) {
						int b = (a + 1) % DIM;
						int c = (a + 2) % DIM;
						tempSave += g*g*(SU2B[i][j][a] * SU2W[k][a] - SU2B[i][k][a] * SU2W[j][a]);
						tempSave -= g3*SU2W[i][a] * (SU2W[j][b] * SU2W[k][c] - SU2W[j][c] * SU2W[k][b]);
					}
				}
				tempSave *= pre_factor;
				if(den_idx_cs != -1) this->density_data_(x, den_idx_cs) = static_cast<RealSave>(tempSave);
				totalCSnumber += tempSave;
			}
		}
		this->parallel_.Sum(totalCSnumber);
		if(dt_idx_cs != -1) this->data_table_.append_value(dt_idx_cs, totalCSnumber);
		return;
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::CalcMagneticEnergy(const int time_step){
		this->CalcSimple(time_step, "MagneticEnergy", "MagneticEnergy",
					[this](const Site<DIM>& x, const int time_step)
					{return this->sc_MagneticEnergy(x, time_step);} );
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::CalcMagneticHelicity(const int time_step){
		const auto T = time_step % CYCLE;
		Site<DIM> x(this->evo_.get_lattice());

		auto dt_idx_helicity = find_index(this->data_table_names_, "MagneticHelicity");
		auto den_idx_helicity = find_index(this->density_names_, "MagneticHelicity");

		Real totalHelicity = 0;
		for (x.first(); x.test(); x.next()) {
			Real tempSave = 0;
			if (this->evo_.is_boundary(x)) {
				if (den_idx_helicity != -1) this->density_data_(x, den_idx_helicity) = 0.0;
				continue;
			}
			else {
					for (int i = 0; i < DIM; ++i) {
						int j = (i + 1) % DIM;
						int k = (i + 2) % DIM;
						tempSave += EM_A_SITE(x, T, i, phi_, U_, V_) *EM_F_SITE(x, T, j, k, phi_, U_, V_);
					}
					tempSave *= DX3;
				totalHelicity += tempSave;
			}
			if (den_idx_helicity != -1) 
				this->density_data_(x, den_idx_helicity) = static_cast<RealSave>(tempSave);
		}
		this->parallel_.Sum(totalHelicity);
		if(dt_idx_helicity != -1) 
			this->data_table_.append_value(dt_idx_helicity, totalHelicity);
		return;
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::CalcElectricEnergy(const int time_step) {
		this->CalcSimple(time_step, "ElectricEnergy", "ElectricEnergy", 
						std::bind(&ElectroweakObserver::sc_ElectricEnergy, this, 
						std::placeholders::_1, std::placeholders::_2) );
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::CalcGaussConstraint(const int time_step) {
		this->CalcSimple(time_step, "GaussConstraint", "GaussConstraint",
						std::bind(&ElectroweakObserver::sc_HamiltonianGC, this, 
						std::placeholders::_1, std::placeholders::_2) );
	}


	template<int DIM>
	void ElectroweakObserver<DIM>::CalcHiggsMagnitude2(const int time_step){
		const auto T = time_step % CYCLE;
		Site<DIM> x(this->evo_.get_lattice());

		auto dt_idx_h2 = find_index(this->data_table_names_, "HiggsMagnitude2");
		auto den_idx_h2 = find_index(this->density_names_, "HiggsMagnitude2");

		Real totalOutput = 0;
		for (x.first(); x.test(); x.next()) {
			if (this->evo_.is_boundary(x)) {
				if (den_idx_h2 != -1) this->density_data_(x, den_idx_h2) = 0.0;
				continue;
			}
			else {
				Real tempSave = phi_(x, T).squaredNorm();
				//if (tempSave < VECN_THRESHOLD) totalOutput += 1.0;
				totalOutput += tempSave;
				if (den_idx_h2 != -1) 
					this->density_data_(x, den_idx_h2) = static_cast<RealSave>(tempSave);
			}
		}
		this->parallel_.Sum(totalOutput);
		auto ave = totalOutput / this->evo_.get_lattice().get_visible_global_sites();
		if(dt_idx_h2 != -1) this->data_table_.append_value(dt_idx_h2, ave);
		return;
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::CalcMinHiggsMagnitude2(const int time_step){
		const auto nowTime = time_step % CYCLE;
		Site<DIM> x(this->evo_.get_lattice());

		auto dt_idx_minh2 = find_index(this->data_table_names_, "MinHiggsMagnitude2");
		auto den_idx_minh2 = find_index(this->density_names_, "MinHiggsMagnitude2");

		/*This function does not support store density info yet*/
		if(den_idx_minh2 != -1) 
			std::cerr << "Density names error. (MinHiggsMagnitude2)" << std::endl;

		auto min = v2;
		for (x.first(); x.test(); x.next()) {
			if (this->evo_.is_boundary(x)) {
				continue;
			}
			else {
				auto phi2 = phi_(x, nowTime).squaredNorm();
				min = phi2 < min ? phi2 : min;
			}
		}
		this->parallel_.Min(min);
		if(dt_idx_minh2 != -1) this->data_table_.append_value(dt_idx_minh2, min);
		return;
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::CalcHiggsWinding(const int time_step){
		this->CalcSimple(time_step, "HiggsWinding", "HiggsWinding", 
					std::bind(&ElectroweakObserver::sc_HiggsWinding, this, 
					 std::placeholders::_1, std::placeholders::_2) );
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::CalcSimple(const int time_step, 
			const std::string& data_table_key,
			const std::string& density_key,
			const std::function<Real(const Site<DIM>&, const int)> sc_quantity) {
		const auto T = time_step % CYCLE;
		Site<DIM> x(this->evo_.get_lattice());

		auto dt_idx = find_index(this->data_table_names_, data_table_key);
		auto den_idx = find_index(this->density_names_, density_key);

		Real total = 0;
		for (x.first(); x.test(); x.next()) {
			if (this->evo_.is_boundary(x)) {
				if (den_idx != -1) this->density_data_(x, den_idx) = 0.0;
				continue;
			}
			else {
				Real tempSave = sc_quantity(x, T);
				total += tempSave * DX3;
				if (den_idx != -1) 
					this->density_data_(x, den_idx) = static_cast<RealSave>(tempSave);
			}
		}
		this->parallel_.Sum(total);
		if(dt_idx != -1)
			this->data_table_.append_value(dt_idx, total);
		return;
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::Measure() {
		auto time_step = this->evo_.get_time_step();
		if (this->data_table_flags_ & ObserverFlags::OBS_TotalEnergy) {
			this->CalcEnergy(time_step);
		}
		if (this->data_table_flags_ & ObserverFlags::OBS_CSNumber) {
			this->CalcCSNumber(time_step);
		}
		if (this->data_table_flags_ & ObserverFlags::OBS_MagneticHelicity) {
			this->CalcMagneticHelicity(time_step);
		}
		if (this->data_table_flags_ & ObserverFlags::OBS_MagneticEnergy) {
			this->CalcMagneticEnergy(time_step);
		}
		if (this->data_table_flags_ & ObserverFlags::OBS_ElectricEnergy) {
			this->CalcElectricEnergy(time_step);
		}
		if (this->data_table_flags_ & ObserverFlags::OBS_GaussConstraint) {
			this->CalcGaussConstraint(time_step);
		}
		if (this->data_table_flags_ & ObserverFlags::OBS_HiggsMagnitude2) {
			this->CalcHiggsMagnitude2(time_step);
		}
		if (this->data_table_flags_ & ObserverFlags::OBS_MinHiggsMagnitude2) {
			this->CalcMinHiggsMagnitude2(time_step);
		}
		if (this->data_table_flags_ & ObserverFlags::OBS_HiggsWinding) {
			this->CalcHiggsWinding(time_step);
		}
		return;
	}

	template<int DIM>
	Real ElectroweakObserver<DIM>::sc_GaussU1(const Site<DIM>& x, const int T) const {
		Real u1term = 0;
		for (int i = 0; i < DIM; ++i) {
			u1term += E_(x, i, T).imag() - E_(x - i, i, T).imag();
		}
		return u1term / DX - gp*(pi_(x, T).dot(phi_(x, T))).imag();
	}

	template<int DIM>
	Real ElectroweakObserver<DIM>::sc_GaussSU2(const Site<DIM>& x, const int T,
			const int su2dir) const {
		Real su2term = 0;
		for (int i = 0; i < DIM; ++i) {
			su2term += (I*Pauli[su2dir] * (F_(x, i, T) - U_(x - i, i, T).adjoint()*F_(x - i, i, T)*U_(x - i, i, T))).trace().real();
		}
		return su2term / DX - g*(pi_(x, T).dot(iPauli[su2dir] * phi_(x, T))).real();
	}

	template<int DIM>
	Real ElectroweakObserver<DIM>::sc_HamiltonianGC(const Site<DIM>& x, const int T) const {
		return 0.5*(pow(sc_GaussU1(x, T), 2) + pow(sc_GaussSU2(x, T, 0), 2) + pow(sc_GaussSU2(x, T, 1), 2) + pow(sc_GaussSU2(x, T, 2), 2));
	}

	template<int DIM>
	Real ElectroweakObserver<DIM>::sc_MagneticEnergy(const Site<DIM>& x, const int T) const {
		return 0.5*(pow(EM_F_SITE(x, T, 0, 1, phi_, U_, V_), 2)
				+ pow(EM_F_SITE(x, T, 1, 2, phi_, U_, V_), 2)
				+ pow(EM_F_SITE(x, T, 0, 2, phi_, U_, V_), 2));
	}


	template<int DIM>
	Real ElectroweakObserver<DIM>::sc_ElectricEnergy(const Site<DIM>& x, const int T) const {
		return 0.5*(
		pow(EM_E_SITE(x, T, 0, phi_, U_, V_, pi_, F_, E_), 2) +
		pow(EM_E_SITE(x, T, 1, phi_, U_, V_, pi_, F_, E_), 2) +
		pow(EM_E_SITE(x, T, 2, phi_, U_, V_, pi_, F_, E_), 2));
	}

	template<int DIM>
	Real ElectroweakObserver<DIM>::sc_HiggsWinding(const Site<DIM>& x, const int T) const{
		Cmplx hw = 0;
		const SU2matrix Uphi_p[DIM] = {
			SU2invphi(x + X_AXIS,T,phi_),
			SU2invphi(x + Y_AXIS,T,phi_),
			SU2invphi(x + Z_AXIS,T,phi_) };
		const SU2matrix Uphi_m[DIM] = {
			SU2invphi(x - X_AXIS,T,phi_),
			SU2invphi(x - Y_AXIS,T,phi_),
			SU2invphi(x - Z_AXIS,T,phi_) };
		auto Uphi = SU2invphi(x, T, phi_);
		for (int i = 0; i < DIM; ++i) {
			const int j = (i + 1) % DIM;
			const int k = (i + 2) % DIM;
			hw += ((Uphi_p[i] - Uphi_m[i]).adjoint()*(
				(Uphi_p[j] - Uphi_m[j])*(Uphi_p[k] - Uphi_m[k]).adjoint() -
				(Uphi_p[k] - Uphi_m[k])*(Uphi_p[j] - Uphi_m[j]).adjoint()
				)*Uphi).trace();
		}
		hw *= 1.0 / (24.0*8.0*PI*PI*DX3);
		if (std::isnan(hw.real())) {
			return 0.0;
		}
		else {
			//somehow there needs to be a minus sign.
			return hw.real();
		}
	}

	template<int DIM>
	void ElectroweakObserver<DIM>::SaveDensityData(const std::string& file_name, 
							const int freq, 
							const int starting_time_step) const {
		auto t = this->evo_.get_time_step();
		if( ((t - starting_time_step) >= 0) && ((t - starting_time_step)%freq == 0) ){
			this->density_data_.write(file_name, &(this->density_names_));
		}
		return;
	}



}

#endif