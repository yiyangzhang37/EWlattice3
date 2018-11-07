#ifndef EWMODEL_H
#define EWMODEL_H

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>

#include "../ParaSite/ParaSite.h"

#include "EW_parameter.h"
#include "EW_functions.h"

namespace HDF5_Wrapper {
	template<>
    inline hid_t get_H5_datatype<Electroweak::SU2vector>(){
		hsize_t arr_size = 4;
		return H5Tarray_create(
			get_H5_datatype<Electroweak::Real>(),  /*base type*/
			1, /*rank of the array*/
			&arr_size /*size of each array dimension*/
		);
	}

	template<>
    inline hid_t get_H5_datatype<Electroweak::SU2matrix>(){
		hsize_t arr_size = 8;
		return H5Tarray_create(
			get_H5_datatype<Electroweak::Real>(),
			1,
			&arr_size);
	}

	template<>
	inline hid_t get_H5_datatype<Electroweak::U1matrix>(){
		hsize_t arr_size = 2;
		return H5Tarray_create(
			get_H5_datatype<Electroweak::Real>(),
			1,
			&arr_size);
	}
}

/*
The ElectroweakEvolution class shall be the base class,
implementing the basic evolution and boundary condition.
*/

namespace Electroweak{

    using namespace ParaSite;
	using namespace MPI_Wrapper;

	constexpr int PHI_FIELD = 1 << 0;
	constexpr int PI_FIELD = 1 << 1;
	constexpr int U_FIELD = 1 << 2;
	constexpr int F_FIELD = 1 << 3;
	constexpr int V_FIELD = 1 << 4;
	constexpr int E_FIELD = 1 << 5;
	constexpr int ALL_FIELDS = PHI_FIELD | PI_FIELD | U_FIELD | F_FIELD | V_FIELD | E_FIELD;
    
    template<int DIM>
    class ElectroweakEvolution{

	protected:
		/*member variables*/
		const Lattice<DIM>& lat_;

		const Parallel2D& parallel_;

		Field< SU2vector, DIM > phi_;
		Field< SU2vector, DIM > pi_;
		Field< SU2matrix, DIM > U_;
		Field< SU2matrix, DIM > F_;
		Field< U1matrix, DIM > V_;
		Field< U1matrix, DIM> E_;

		int time_step_ = 0;

		Parameters param_;

		std::string id_ = "";
		time_t start_time_;
		time_t finish_time_;

		/*
		Here are the two common ways to define boundary region (or constraint region).
		periodic boundary is just no boundary at all.
		*/
		constexpr int test_periodic_boundary(const Site<DIM>& x) const {
			return 0;
		}
		int test_one_layer_boundary(const Site<DIM>& x) const {
			int isBD = 0;
			for (int i = 0; i < DIM; ++i) {
				isBD += (x.coord(i) == 0);
				isBD += 2 * (x.coord(i) == nSize[i] - 1);
			}
			return isBD;
		}

		void record_basic_parameters();

    public:
		ElectroweakEvolution(const Lattice<DIM>& lat, const Parallel2D& parallel, const std::string& id);
        virtual ~ElectroweakEvolution() {}

		virtual void RecordParameters();
		void SaveParameters(const std::string& name, const bool is_root_save_only = true) const;
		
		/*
		Measure the finish time and compute the time cost.
		Save all the parameters again, including the total time cost.
		*/
		void ConcludeEvolution(const std::string& name, const bool is_root_save_only = true);

        /* evolution functions */
		void TimeAdvance() { this->time_step_++; }

        void UpdateFields() const;
        void EvolveInterior_KS() const;
        void EvolveInterior_RadialDamping() const;
        void EvolveScalarOneSite_KS(const Site<DIM>& x, const int nowTime, const int futureTime) const;
        void EvolveScalarOneSite_RadialDamping(const Site<DIM>& x, const int nowTime, const int futureTime) const;
        void EvolveU1OneSite_KS(const Site<DIM>& x, const int component, const int nowTime, const int futureTime) const;
        void EvolveSU2OneSite_KS(const Site<DIM>& x, const int component, const int nowTime, const int futureTime) const;

		void SaveFields(const std::string& filename, const int field_type = ALL_FIELDS) const;

		/*
		The constraint region and boundary region can be modified 
		by specific problems.
		*/
		virtual const int is_constraint_region(const Site<DIM>& x) const {
			return test_periodic_boundary(x);
		}
		virtual const int is_boundary(const Site<DIM>& x) const {
			return test_periodic_boundary(x);
		}

		/* boundary conditions */

		/* trivial initial condition*/
		void TrivialInitialCondition() const;

		/* get fields for computations*/
		int get_time_step() const { return this->time_step_; }
		const Lattice<DIM>& get_lattice() const { return this->lat_; }
		const Parallel2D& get_parallel_object() const { return this->parallel_; }
		auto get_phi() const -> const decltype(phi_)& {return this->phi_; }
		auto get_pi() const -> const decltype(pi_)& {return this->pi_; }
		auto get_U() const -> const decltype(U_)& {return this->U_; }
		auto get_F() const -> const decltype(F_)& {return this->F_; }
		auto get_V() const -> const decltype(V_)& {return this->V_; }
		auto get_E() const -> const decltype(E_)& {return this->E_; }
    };

    template<int DIM>
	ElectroweakEvolution<DIM>::ElectroweakEvolution(
        const Lattice<DIM>& lat, 
        const Parallel2D& parallel, 
        const std::string& id)
        :
        lat_(lat),
        parallel_(parallel),
		phi_(lat, CYCLE, parallel),
		pi_(lat, CYCLE, parallel),
		U_(lat, DIM, CYCLE, parallel),
		F_(lat, DIM, CYCLE, parallel),
		V_(lat, DIM, CYCLE, parallel),
		E_(lat, DIM, CYCLE, parallel),
		param_(),
        id_(id),
		start_time_(std::time(nullptr)),
		finish_time_(std::time(nullptr))
    {       
    }

	template<int DIM>
	void ElectroweakEvolution<DIM>::RecordParameters() {
		this->record_basic_parameters();
	}

	template<int DIM>
	void ElectroweakEvolution<DIM>::record_basic_parameters(){
		std::stringstream ss1, ss2;
    	ss1 << std::ctime(&start_time_);
		param_.add("ID", this->id_);
		param_.add("Comment", ReadConfigIni("COMMENT"));
		param_.add("Start Time", ss1.str());

		param_.add("Dimension", DIM);
		param_.add("HaloLayer", halo);
		param_.add("LatticeSize(X)", nSize[0]);
		param_.add("LatticeSize(Y)", nSize[1]);
		param_.add("LatticeSize(Z)", nSize[2]);
		param_.add("PresetMaxTimeSteps", Ntimesteps);
		param_.add("DX", DX);
		param_.add("DT", DT);
		param_.add("HiggsCoupling(lambda)", lambda);
		param_.add("HiggsVEV(eta)", v);
		param_.add("SU2GaugeCoupling(g)", g);
		param_.add("WeinbergAngle(sin2w)", sin2w);
		param_.add("HiggsDamping", DAMPING_HIGGS);
		param_.add("GaugeDamping", DAMPING_GAUGE);
		param_.add("FermiFamily", FermiFamily);
		param_.add("MPI_GridSize", std::to_string(lat_.get_node_size()[0]) + "*" + std::to_string(lat_.get_node_size()[1]));
		return;
	}

	template<int DIM>
	void ElectroweakEvolution<DIM>::SaveParameters(
		const std::string& name,
		const bool is_root_save_only) const {
		if (is_root_save_only) {
			if(this->parallel_.is_root()) param_.save(name);
		}
		else {
			param_.save(name);
		}
		return;
	}

	template<int DIM>
	void ElectroweakEvolution<DIM>::ConcludeEvolution(
		const std::string& name, 
		const bool is_root_save_only) {
		this->finish_time_ = std::time(nullptr);
		std::stringstream ss2;
		ss2 << std::ctime(&finish_time_);
		param_.add("Finish Time", ss2.str());
		param_.add("Time Cost(s)", std::difftime(finish_time_, start_time_));
		this->SaveParameters(name, is_root_save_only);
		return;
	}


	template<int DIM>
	void ElectroweakEvolution<DIM>::UpdateFields() const {
		auto nowTime = this->time_step_ % CYCLE;
		auto futureTime = (this->time_step_ + 1) % CYCLE;
		Site<DIM> x(this->lat_);
		for (x.first(); x.test(); x.next()) {
			phi_(x, futureTime) = phi_(x, nowTime) + DT * pi_(x, nowTime);
			for (int i = 0; i<DIM; ++i) {
				V_(x, i, futureTime) = 1.0 / UNITY_E * E_(x, i, nowTime) * V_(x, i, nowTime);
				U_(x, i, futureTime) = 1.0 / UNITY_F * F_(x, i, nowTime) * U_(x, i, nowTime);
			}
		}
		phi_.update_halo();
		U_.update_halo();
		V_.update_halo();
		return;
	}

	template<int DIM>
	void ElectroweakEvolution<DIM>::EvolveInterior_KS() const {
		auto nowTime = this->time_step_ % CYCLE;
		auto futureTime = (this->time_step_ + 1) % CYCLE;
		Site<DIM> x(this->lat_);
		//update pi,F,E fields to t+3/2*dt time, according to evolution equations
		for (x.first(); x.test(); x.next()) {
			//test isBoundary
			auto isBD = this->is_constraint_region(x);
			if (isBD > 1) continue;//boundary: isBoundary>1
			else if (isBD == 1) {//a special case: one component of gauge fields is inside
				int dirInside = -1; //find the inside direction
				for (int i = 0; i < DIM; ++i) {
					if (x.coord(i) == 0) { dirInside = i; break; }
				}
				//evolve gauge fields as interior way
				EvolveU1OneSite_KS(x, dirInside, nowTime, futureTime);
				EvolveSU2OneSite_KS(x, dirInside, nowTime, futureTime);
			}
			//Evolve interior points
			else {
				//scalar field
				EvolveScalarOneSite_KS(x, nowTime, futureTime);
				//gauge fields
				for (int k = 0; k < DIM; ++k) {
					EvolveU1OneSite_KS(x, k, nowTime, futureTime);
					EvolveSU2OneSite_KS(x, k, nowTime, futureTime);
				}
			}
		}
		pi_.update_halo();
		F_.update_halo();
		E_.update_halo();
		return;
	}

	template<int DIM>
	void ElectroweakEvolution<DIM>::EvolveInterior_RadialDamping() const {
		int nowTime = this->time_step_ % CYCLE;
		int futureTime = (this->time_step_ + 1) % CYCLE;
		Site<DIM> x(this->lat_);
		//update pi,F,E fields to t+3/2*dt time, according to evolution equations
		for (x.first(); x.test(); x.next()) {
			//test isBoundary
			int isBD = this->is_constraint_region(x);
			if (isBD > 1) continue;//boundary: isBoundary>1
			else if (isBD == 1) {//a special case: one component of gauge fields is inside
				int dirInside = -1; //find the inside direction
				for (int i = 0; i < DIM; ++i) {
					if (x.coord(i) == 0) { dirInside = i; break; }
				}
				//evolve gauge fields as interior way
				EvolveU1OneSite_KS(x, dirInside, nowTime, futureTime);
				EvolveSU2OneSite_KS(x, dirInside, nowTime, futureTime);
			}
			//Evolve interior points
			else {
				//scalar field
				EvolveScalarOneSite_RadialDamping(x, nowTime, futureTime);
				//gauge fields
				for (int k = 0; k < DIM; ++k) {
					EvolveU1OneSite_KS(x, k, nowTime, futureTime);
					EvolveSU2OneSite_KS(x, k, nowTime, futureTime);
				}
			}
		}
		pi_.update_halo();
		F_.update_halo();
		E_.update_halo();
		return;
	}

	template<int DIM>
	void ElectroweakEvolution<DIM>::EvolveScalarOneSite_KS(
		const Site<DIM>& x, 
		const int nowTime, 
		const int futureTime) const {
		SU2vector gradientPart;
		gradientPart << 0, 0;
		for (int i = 0; i < DIM; ++i) {
			gradientPart += U_(x, i, futureTime)*V_(x, i, futureTime)*phi_(x + i, futureTime)
				- 2.0 * phi_(x, futureTime)
				+ U_(x - i, i, futureTime).adjoint() * conj(V_(x - i, i, futureTime)) * phi_(x - i, futureTime);
		}
		gradientPart /= DX2;
		pi_(x, futureTime) = (1.0 - DAMPING_HIGGS * DT)*pi_(x, nowTime) + DT * (gradientPart - DVDphi(x, phi_, futureTime));
		return;
	}

	template<int DIM>
	void ElectroweakEvolution<DIM>::EvolveScalarOneSite_RadialDamping(
		const Site<DIM>& x,
		const int nowTime,
		const int futureTime) const {
		SU2vector gradientPart;
		gradientPart << 0, 0;
		for (int i = 0; i<DIM; ++i) {
			gradientPart += U_(x, i, futureTime)*V_(x, i, futureTime)*phi_(x + i, futureTime)
				- 2.0 * phi_(x, futureTime)
				+ U_(x - i, i, futureTime).adjoint() * conj(V_(x - i, i, futureTime)) * phi_(x - i, futureTime);
		}
		gradientPart /= DX2;
		pi_(x, futureTime) = pi_(x, nowTime) + DT * (gradientPart - DVDphi(x, phi_, futureTime));
		auto phi_norm_now = phi_(x, nowTime).norm();
		auto phi_norm_future = phi_(x, futureTime).norm();
		//const auto phi_norm_now = sqrt(std::norm(phi(x, nowTime)(0)) + std::norm(phi(x, nowTime)(1)));
		//const auto phi_norm_future = sqrt(std::norm(phi(x, futureTime)(0)) + std::norm(phi(x, futureTime)(1)));
		if (phi_norm_future > 1e-6 && phi_norm_now > 1e-6) {
			pi_(x, futureTime) -= DAMPING_HIGGS * phi_(x, futureTime)*log(phi_norm_future / phi_norm_now);
		}
		return;
	}

	template<int DIM>
	void ElectroweakEvolution<DIM>::EvolveU1OneSite_KS(
		const Site<DIM>& x,
		const int component,
		const int nowTime,
		const int futureTime) const {
		Cmplx cmpImFutureE(0, 0);
		Real reFutureE = 0.0, imFutureE = 0.0;
		U1matrix plaq = PlaquetteSum4(x, V_, component, futureTime);
		cmpImFutureE = (1.0 - DAMPING_GAUGE * DT) * E_(x, component, nowTime).imag()
			+ DT * (gp / DX * (phi_(x + component, futureTime).dot(U_(x, component, futureTime).adjoint() * conj(V_(x, component, futureTime)) * phi_(x, futureTime))).imag() - 2.0 / (gp*DX3) * plaq.imag());
		if (cmpImFutureE.imag() > EVO_ERROR) { std::cerr << "U1 FIELD ERROR:COMPLEX VALUE" << std::endl; exit(EXIT_FAILURE); }
		else imFutureE = cmpImFutureE.real();
		if (pow(imFutureE, 2) <= pow(UNITY_E, 2)) reFutureE = sqrt(pow(UNITY_E, 2) - pow(imFutureE, 2));
		else { std::cerr << "U1 FIELD ERROR: NON-UNITARY." << std::endl; exit(EXIT_FAILURE); }
		E_(x, component, futureTime) = U1matrix(reFutureE, imFutureE);
		return;
	}

	template<int DIM>
	void ElectroweakEvolution<DIM>::EvolveSU2OneSite_KS(
		const Site<DIM>& x,
		const int component,
		const int nowTime,
		const int futureTime) const {
		std::vector<Cmplx> su2_Fcomp(SU2DIM, 0);
		Real su2_Fsum = 0.0;
		Real su2_F0 = 0.0;
		Cmplx plaqTrace;
		for (int a = 0; a<SU2DIM; ++a) {  //su2 index
			plaqTrace = (I*Pauli[a] * PlaquetteSum4(x, U_, component, futureTime)).trace();
			if (plaqTrace.imag() > EVO_ERROR) { std::cerr << "SU2 FIELD ERROR:COMPLEX VALUE" << std::endl; exit(EXIT_FAILURE); }
			su2_Fcomp[a] = (1.0 - DAMPING_GAUGE * DT) * (I * Pauli[a] * F_(x, component, nowTime)).trace()
				+ DT * (g / DX * (phi_(x + component, futureTime).dot(U_(x, component, futureTime).adjoint() * conj(V_(x, component, futureTime)) * I*Pauli[a] * phi_(x, futureTime))).real()
					- 1.0 / (g*DX3) * plaqTrace);
		}

		for (int a = 0; a<SU2DIM; ++a) {
			if (su2_Fcomp[a].imag() > EVO_ERROR) { std::cerr << "SU2 FIELD ERROR:COMPLEX VALUE" << std::endl; exit(EXIT_FAILURE); }
			su2_Fsum += pow(su2_Fcomp[a].real() / 2.0, 2);
		}
		if (su2_Fsum <= pow(UNITY_F, 2)) su2_F0 = sqrt(pow(UNITY_F, 2) - su2_Fsum);
		else {
			std::cerr << "SU2 FIELD ERROR: NON-UNITARY." << std::endl;
			std::cout << "( " << x.coord(0) << " , " << x.coord(1) << " , " << x.coord(2) << " )" << std::endl;
			std::cout << su2_Fsum << std::endl;
			getchar();  exit(EXIT_FAILURE);
		}
		F_(x, component, futureTime) = su2_F0 * Ident;
		for (int a = 0; a<SU2DIM; ++a) {
			F_(x, component, futureTime) += -I / 2.0 * Pauli[a] * su2_Fcomp[a].real();
		}
		return;
	}

	template<int DIM>
	void ElectroweakEvolution<DIM>::TrivialInitialCondition() const {
		Site<DIM> x(this->lat_);
		for (auto t = 0; t < CYCLE; ++t) {
			for (x.first(); x.test(); x.next()) {
				phi_(x, t) = SU2vector(0, 0);
				pi_(x, t) = SU2vector(0, 0);
				for (auto i = 0; i < DIM; ++i) {
					U_(x, i, t) = Ident;
					F_(x, i, t) = UNITY_F * Ident;
					V_(x, i, t) = Cmplx(1, 0);
					E_(x, i, t) = UNITY_E * Cmplx(1, 0);
				}
			}
		}
		this->phi_.update_halo();
		this->pi_.update_halo();
		this->U_.update_halo();
		this->F_.update_halo();
		this->V_.update_halo();
		this->E_.update_halo();
		return;
	}

	template<int DIM>
	void ElectroweakEvolution<DIM>::SaveFields(const std::string& filename, const int field_type) const{
		
		if(field_type & PHI_FIELD) this->phi_.write("phi_" + filename, nullptr, "phi");
		if(field_type & PI_FIELD) this->pi_.write("pi_" + filename, nullptr, "pi");
		if(field_type & U_FIELD) this->U_.write("U_" + filename, nullptr, "U");
		if(field_type & V_FIELD) this->V_.write("V_" + filename, nullptr, "V");
		if(field_type & F_FIELD) this->F_.write("F_" + filename, nullptr, "F");
		if(field_type & E_FIELD) this->E_.write("E_" + filename, nullptr, "E");
		return;
	}

}



#endif