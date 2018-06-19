#ifndef EW_BUBBLE_NUCL_H
#define EW_BUBBLE_NUCL_H

#include "EW_Base.h"
#include <random>
#include <type_traits>

namespace Electroweak{
	//extend the ObserverFlags
	namespace ObserverFlags {
		const FlagType OBS_NewBubbleCount = 1 << 31; //"NewBubbleCount"
	}
}

namespace EW_BubbleNucleation{

    using namespace Electroweak;

    /*The probability that one site can nucleate*/
	const Real NUCLEATION_PROB = std::stod(ReadConfigIni("NULEATION_PROBABILITY"));
	/*Higgs magnitude larger than this value cannot nucleate*/
	const Real NUCLEATION_LIMIT = 0.0001*v2;
	/*Hard cutoff for the nucleation area, should be > 3 * HIGGS_STD*/
	const unsigned int NUCLEATION_RADIUS_SITE = std::stoi(ReadConfigIni("BUBBLE_LAT_RADIUS"));
	const Real NUCLEATION_RADIUS_LEN = NUCLEATION_RADIUS_SITE*DX;
	
	const int BUBBLES_HALF_SEP = std::stoi(ReadConfigIni("BUBBLE_HALF_SEPARATION"));
	// time offset due to the starting point of bubble collision at different time 
	const int RECORD_TIME_OFFSET = std::stoi(ReadConfigIni("RECORD_TIME_OFFSET"));

	
	// The threshold value of Phi^2/eta^2.
	constexpr Real EARLY_STOP_LIMIT = 0.25; 
	// The mean average steps taken to calculate the early stopping value.
	constexpr int STOP_MEAN_STEPS = 10;

	template<int DIM>
    class BubbleNucleation : public ElectroweakEvolution<DIM> {
    public:
		BubbleNucleation(
			const Lattice<DIM>& lat, 
			const Parallel2D& parallel, 
			const std::string& id);
		BubbleNucleation() = default;

		void RecordParameters() override;
		
		/*Initial condition*/
		void InitializeSymmetricPhase() const;

		/*nulceation*/
		/*
		The overall control of random bubble nucleation.
		The update_halo() will be called at the end of this function for phi_ field.
		Return: new bubble count.
		*/
        int RandomBubbleNucleation();

		void OneBubbleTest() const;
        void TwoBubblesTest(const int sep) const;
        void NonRandomTest(const int half_sep);
		
		//get the new bubble count for the current time step.
		int GetNewBubbleCount() const {return this->new_bubbles_count_;}

		//check whether the simulation satisfies the early stopping condition.
		//return true if the condition is satisfied, and can perform early stop.
		//otherwise return false.
		bool CheckEarlyStop(const ElectroweakObserver<DIM>& obs) const;

		//check whether the Higgs field has all left the symmetric phase.
		//by looking at the MinHiggsMagnitude2
		//if yes, then skip the random nucleation.
		///This function purely for efficiency.
		bool CheckHiggsAllInBrokenPhase(const ElectroweakObserver<DIM>& obs) const;

    private:
		/*
		check whether each site is in symmtric phase (0) or not (1), for the time slice nowTime.
		The threshold is set by phi2_limit.
		The infomation will be gathered to the root process. to a single vector higgs_phase_info.
		*/
		void GatherHiggsPhaseInfo(
			const int nowTime,
			const double phi2_limit,
			std::vector<char>& higgs_phase_info) const;
		/*
		Reorder HiggsPhaseInfo to the correct global index order.
		*/
		void ReorderHiggsPhaseInfo(std::vector<char>& higgs_phase_info) const;
		
		/*
		Verify if the preselected positions are legible to nucelate bubbles.
		Whether it's in the symmetric region, or some of them are overlapped.
		*/
		void VerifyPickedPositions(
			const std::vector<char>& sorted_higgs_phase_info,
			const std::vector<IndexType>& preselected_positions,
			std::vector<IndexType>& picked_positions) const;
		
		/*
		Nucleate one bubble with the bubble center specified with global_coord,
		and the randomized Higgs angluar components specified by ha[4].
		*/
		void NucleateOneBubble_Exp(
			const int nowTime,
			const IndexType global_index, 
			const SU2vector& phi_hat) const;
		
		/*
		Wrapper
		*/
		void NucleateOneBubble(
			const int nowTime,
			const IndexType* global_coord, 
			const SU2vector& phi_hat) const{
			this->NucleateOneBubble_Exp(nowTime, 
					this->lat_.global_coord2index(global_coord), 
					phi_hat);
		}

		/*
		Given a global index, get the info about that site:
		a list of all the sites (global index) within the bubble.
		Will assume DIM == 3.
		*/
		void GetBubbleRegion(
			const IndexType global_index, 
			const int bubble_lat_radius,
			std::vector< IndexType >& region_list) const;
		
		/*
		Generate a random Higgs direction for the NucleateOneBubble_Exp(...) function to use.
		This function takes place at each process.
		The random seed is broadcasted to every process.
		*/
		void GenerateRandomHiggsComponents(SU2vector& phi_hat) const;

		/*
		stores the number of new bubbles generated at each time step.
		This number is already the sum over all the processes.
		*/
        int new_bubbles_count_ = 0;
	};


	template<int DIM>
	class NucleationObserver : public ElectroweakObserver<DIM> {
	public:
		NucleationObserver(const BubbleNucleation<DIM>& bubble);
		~NucleationObserver() = default;
		
		void Measure() override;

		void CalcNewBubbleCount(const int time_step);
	
	protected:

		const BubbleNucleation<DIM>& nucl_;

		void init_name_vector(
			const FlagType flags,
			std::vector<std::string>& names) const override;

	};

	template<int DIM>
	BubbleNucleation<DIM>::BubbleNucleation(
		const Lattice<DIM>& lat,
		const Parallel2D& parallel,
		const std::string& id)
		:
		ElectroweakEvolution<DIM>(lat, parallel, id) {

	}

	template<int DIM>
	void BubbleNucleation<DIM>::RecordParameters() {
		this->record_basic_parameters();
		this->param_.add("NucleationProbability", NUCLEATION_PROB, true);
		this->param_.add("NucleationLimit", NUCLEATION_LIMIT);
		this->param_.add("BubbleRadius(in units of DX)", NUCLEATION_RADIUS_SITE);
		this->param_.add("BubbleRadius", NUCLEATION_RADIUS_LEN);
		this->param_.add("BubbleFixedHalfSeparation(ifFixed)", BUBBLES_HALF_SEP);
		return;
	}

	template<int DIM>
	void BubbleNucleation<DIM>::InitializeSymmetricPhase() const {
		Site<DIM> x(this->lat_);
		for (auto t = 0; t < CYCLE; ++t) {
			for (x.first(); x.test(); x.next()) {
				this->phi_(x, t) = SU2vector(0, 0);
				this->pi_(x, t) = SU2vector(0, 0);
				for (auto i = 0; i < DIM; ++i) {
					this->U_(x, i, t) = Ident;
					this->F_(x, i, t) = UNITY_F * Ident;
					this->V_(x, i, t) = Cmplx(1, 0);
					this->E_(x, i, t) = UNITY_E * Cmplx(1, 0);
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
	bool BubbleNucleation<DIM>::CheckEarlyStop(const ElectroweakObserver<DIM>& obs) const {
		auto T = this->time_step_;
		if(T > STOP_MEAN_STEPS){
			const auto& dtable = obs.get_data_table();
			const auto& col_vals = dtable.get_column("MinHiggsMagnitude2");
			auto mean = std::accumulate(
						col_vals.end() - STOP_MEAN_STEPS, col_vals.end(), 0.0) / (v2 * STOP_MEAN_STEPS);
			if(mean >= EARLY_STOP_LIMIT) return true;
			else return false;
		} else {
			return false;
		}
	}

	template<int DIM>
	bool BubbleNucleation<DIM>::CheckHiggsAllInBrokenPhase(
		const ElectroweakObserver<DIM>& obs) const {
		auto T = this->time_step_;
		if(T > 50){ //only use this check after 50 timesteps
			const auto& dtable = obs.get_data_table();
			const auto& col_vals = dtable.get_column("MinHiggsMagnitude2");
			auto last_val = *col_vals.end();
			if(last_val >= 10 * NUCLEATION_LIMIT) return true;
			else return false;

		} else{
			return false;
		}
	}

	template<int DIM>
	int BubbleNucleation<DIM>::RandomBubbleNucleation() {
		// I feel the nucleation should be more natural to be after update fields.
		// and the nucleation is performed on the updated fields.
		const auto T = (this->time_step_ + 1) % CYCLE;
		int new_bubble_count = 0;

		//get Higgs field phase info from all the processes.
		auto global_size = this->lat_.get_visible_global_sites();
		std::vector<char> higgs_phase_info(global_size);
		
		GatherHiggsPhaseInfo(T, NUCLEATION_LIMIT, higgs_phase_info);
		
		// check if all the sites are in broken phase.
		int all_broken_phase = 1;
		if(this->parallel_.is_root()){
			for( auto site : higgs_phase_info){
				if(site == 0){
					all_broken_phase = 0;
					break;
				}
			}
		}
		//this->parallel_.Barrier();
		this->parallel_.Broadcast(all_broken_phase, this->parallel_.get_root()); 
		this->parallel_.Barrier();
		if(all_broken_phase) {
			this->new_bubbles_count_ = 0;
			return 0; //All sites in broken phase, return 0.
		}
		
		if(this->parallel_.is_root())
			this->ReorderHiggsPhaseInfo(higgs_phase_info);
		this->parallel_.Barrier();

		//select some number of sites by p_B.
		std::vector<IndexType> preselect_positions;
		std::vector<IndexType> picked_positions;
		if(this->parallel_.is_root()){
			// initialize an ordered site index vector.
			std::vector<IndexType> all_sites(global_size);
			IndexType n = 0;
			std::generate(all_sites.begin(), all_sites.end(), 
							[&n]() { return n++; });
			//random shuffle
			std::random_shuffle(all_sites.begin(), all_sites.end());

			//pick some positions with NUCLEATION_PROB
			auto seed = std::chrono::system_clock::now().time_since_epoch().count();
			std::default_random_engine generator(seed);
			std::uniform_real_distribution<Real> distrib(0.0, 1.0);
			preselect_positions.reserve(static_cast<int>(NUCLEATION_PROB*global_size*1.2));
			for (auto item : all_sites) {
				auto site_nucl_p = distrib(generator);
				if (site_nucl_p <= NUCLEATION_PROB) 
					preselect_positions.push_back(item);
			}
			//verify these locations with higgs_phase_info.
			this->VerifyPickedPositions(
				higgs_phase_info,
				preselect_positions,
				picked_positions);	
		}
		int picked_size = picked_positions.size();

		//this->parallel_.Barrier();
		this->parallel_.Broadcast(picked_size, this->parallel_.get_root());
		if(picked_size == 0) {
			this->new_bubbles_count_ = 0;
			return 0; //No new bubbles: return 0.
		}
		this->parallel_.Barrier();

		//send the picked_positions to all the processes.
		if(! this->parallel_.is_root() ){
			picked_positions.resize(picked_size);
		}
		this->parallel_.Broadcast(picked_positions.data(), 
								picked_size, 
								this->parallel_.get_root());

		//generate Higgs angular components on each process with the same seed.
		//The seed is broadcasted from the root process.
		//specify aligned allocator so that the optimization can work.
		std::vector<SU2vector, Eigen::aligned_allocator<SU2vector>> phi_hat(picked_size);
		for(auto i = 0; i < picked_size; ++i){
			auto ptr = phi_hat.data() + i;
			this->GenerateRandomHiggsComponents(*ptr);
		}

		// each process do nucleation, by information from
		// phi_hat and picked_positions.
		for(auto i = 0; i < picked_size; ++i){
			this->NucleateOneBubble_Exp(T, picked_positions[i], phi_hat[i]);
		}
		this->phi_.update_halo();
		this->new_bubbles_count_ = picked_size;
		return picked_size;
	}

	template<int DIM>
	void BubbleNucleation<DIM>::GatherHiggsPhaseInfo(
			const int nowTime,
			const double phi2_limit,
			std::vector<char>& higgs_phase_info) const{
		higgs_phase_info.resize(this->lat_.get_visible_global_sites());
		auto local_sites = this->lat_.get_visible_local_sites();
		std::vector<char> higgs_phase_local(local_sites);
		Site<DIM> x(this->lat_);
		IndexType i = 0;
		for(x.first(); x.test(); x.next()){
			if( this->phi_(x, nowTime).squaredNorm() < phi2_limit ){
				higgs_phase_local[i] = 0;
			} else {
				higgs_phase_local[i] = 1;
			}
			i++;
		}
		this->parallel_.Gather(higgs_phase_local.data(), local_sites,
								higgs_phase_info.data(), this->parallel_.get_root());
		return;
	}
	
	template<int DIM>
	void BubbleNucleation<DIM>::ReorderHiggsPhaseInfo(
		std::vector<char>& higgs_phase_info) const {
		auto local_sites_per_process = this->lat_.get_visible_local_sites();
		auto local_size = this->lat_.get_local_size();
		auto total_size = higgs_phase_info.size();
		std::vector<IndexType> gidx_list(total_size);
		for(auto i = 0; i < total_size; ++i){
			//get grid_loc
			GridIndexType grid_loc[2];
			auto world_rank = i / local_sites_per_process;
			this->parallel_.world_rank_to_rowcol(world_rank, grid_loc);
			
			//get global index
			auto local_vis_idx = i % local_sites_per_process;
			IndexType lv_coord[DIM], g_coord[DIM];
			this->lat_.local_vis_index2coord(local_vis_idx, lv_coord);
			local_vis_coord_to_global_coord<DIM>(lv_coord, g_coord, grid_loc, local_size);
			auto global_idx = this->lat_.global_coord2index(g_coord);

			gidx_list[i] = global_idx;
		}
		//sort higgs_phase_info by gidx_list
		std::vector<char> sorted_higgs_info(total_size);
		for(auto i = 0; i < total_size; ++i){
			sorted_higgs_info[gidx_list[i]] = higgs_phase_info[i];
		}
		higgs_phase_info = sorted_higgs_info;
		return;
	}

	template<int DIM>
	void BubbleNucleation<DIM>::VerifyPickedPositions(
			const std::vector<char>& sorted_higgs_phase_info,
			const std::vector<IndexType>& preselected_positions,
			std::vector<IndexType>& picked_positions) const {
		picked_positions.clear();
		std::vector<IndexType> region_list;
		//The modification of higgs_info will be done on a copy.
		auto cp_sorted_info = sorted_higgs_phase_info;
		for(auto gidx : preselected_positions){
			this->GetBubbleRegion(gidx, NUCLEATION_RADIUS_SITE, region_list);
			bool can_nucleate = true;
			for(auto s : region_list){
				if(cp_sorted_info[s] == 1){
					can_nucleate = false;
					break;
				}
			}
			if(can_nucleate){
				picked_positions.push_back(gidx);
				for(auto s : region_list){
					cp_sorted_info[s] = 1;
				}
			}
		}
		return;
	}

	template<int DIM>
	void BubbleNucleation<DIM>::NucleateOneBubble_Exp(
		const int nowTime,
		const IndexType global_index, 
		const SU2vector& phi_hat) const {

		std::vector<IndexType> region_list;
		this->GetBubbleRegion(global_index, NUCLEATION_RADIUS_SITE, region_list);

		Site<DIM> x(this->lat_);
		for(auto gidx : region_list){
			//check if gidx is a local visible site
			if(this->lat_.is_local(gidx)){
				//radial part
				auto r = this->lat_.global_lat_distance(global_index, gidx) * DX; //distance to bubble center
				auto mag = (1.0 + pow(sqrt(2.0) - 1.0, 2)) * exp(-mH * r / sqrt(2.0));
				mag /= 1 + pow(sqrt(2.0) - 1, 2)*exp(-sqrt(2.0)*mH*r);
				auto vis_idx = this->lat_.global_index_to_local_vis_index(gidx);
				auto mem_idx = this->lat_.local_vis_index_to_local_mem_index(vis_idx);
				x.set_index(mem_idx);
				this->phi_(x, nowTime) = mag * v * phi_hat;
			} else continue;
		}
		return;
	}

	template<int DIM>
	void BubbleNucleation<DIM>::GetBubbleRegion(
		const IndexType global_index, 
		const int bubble_lat_radius,
		std::vector< IndexType >& region_list) const {

		static_assert(DIM == 3, "Dimension not equal to 3.");
		const auto& r0 = bubble_lat_radius;
		const auto r2 = r0*r0;

		std::vector<IndexType> tmp_list;
		IndexType center_coord[DIM];
		this->lat_.global_index2coord(global_index, center_coord);
		IndexType admitted_coord[DIM];

		for(auto x = -r0; x <= r0; ++x){
			for(auto y = -r0; y <= r0; ++y){
				for(auto z = -r0; z <= r0; ++z){
					if(x*x + y*y + z*z <= r2){
						admitted_coord[0] = (center_coord[0] + x + this->lat_.get_global_size()[0]) 
												% this->lat_.get_global_size()[0];
						admitted_coord[1] = (center_coord[1] + y + this->lat_.get_global_size()[1]) 
												% this->lat_.get_global_size()[1];
						admitted_coord[2] = (center_coord[2] + z + this->lat_.get_global_size()[2]) 
												% this->lat_.get_global_size()[2];
						tmp_list.push_back(this->lat_.global_coord2index(admitted_coord));
					}
					//else continue.
				}
			}
		}
		region_list = tmp_list;
		return;
	}

	template<int DIM>
	void BubbleNucleation<DIM>::GenerateRandomHiggsComponents(SU2vector& phi_hat) const{
		auto seed = std::chrono::system_clock::now().time_since_epoch().count();
		this->parallel_.Broadcast(seed, this->parallel_.get_root());
		
		std::default_random_engine generator(seed);
		std::uniform_real_distribution<Real> distrib1(-1.0, 1.0);
		Real v1, v2, v3, v4;
		Real s1, s2, s3;
		do{
			v1 = distrib1(generator);
			v2 = distrib1(generator);
			s1 = v1*v1 + v2*v2;
		} while( s1 >= 1.0 );
		do{
			v3 = distrib1(generator);
			v4 = distrib1(generator);
			s2 = v3*v3 + v4*v4;
		} while(s2 >= 1.0);
		s3 = sqrt((1.0 - s1) / s2);
		phi_hat(0) = Cmplx(v1, v2);
		phi_hat(1) = s3 * Cmplx(v3, v4);
		return;
	}

	template<int DIM>
	void BubbleNucleation<DIM>::OneBubbleTest() const {
		Site<DIM> x(this->lat_);
		if (this->time_step_ == 0) {
			auto T = (this->time_step_ + 1) % CYCLE;
			SU2vector phi_hat;
			this->GenerateRandomHiggsComponents(phi_hat);
			IndexType global_coord[] = {nSize[0] / 2, nSize[1] / 2, nSize[2] / 2};
			this->NucleateOneBubble(T, global_coord, phi_hat);
			this->phi_.update_halo();
		}
		return;
	}

	template<int DIM>
	void BubbleNucleation<DIM>::TwoBubblesTest(const int half_sep) const {
		Site<DIM> x(this->lat_);
		if (this->time_step_ == 0) {
			auto T = (this->time_step_ + 1) % CYCLE;
			SU2vector phi_hat;
			this->GenerateRandomHiggsComponents(phi_hat);
			IndexType c1[DIM] = { nSize[0] / 2 , nSize[1] / 2 , nSize[2] / 2 - half_sep };
			this->NucleateOneBubble(T, c1, phi_hat);

			this->GenerateRandomHiggsComponents(phi_hat);
			IndexType c2[DIM] = { nSize[0] / 2 , nSize[1] / 2 , nSize[2] / 2 + half_sep };
			this->NucleateOneBubble(T, c2, phi_hat);
			this->phi_.update_halo();
		}
		return;
	}

	template<int DIM>
	void BubbleNucleation<DIM>::NonRandomTest(const int half_sep) {
		Site<DIM> x(this->lat_);
		if (this->time_step_ == 0) {
			auto T = (this->time_step_ + 1) % CYCLE;
			int num_in_line[DIM];
			std::transform(nSize, nSize + DIM,
							num_in_line,
							[half_sep](IndexType x){return x/(2*half_sep); });
			auto total_bubbles = std::accumulate(num_in_line, num_in_line + DIM, 1, 
												std::multiplies<int>());
			for(auto i = 0; i < num_in_line[0]; ++i){
				for(auto j = 0; j < num_in_line[1]; ++j){
					for(auto k = 0; k < num_in_line[2]; ++k){
						IndexType core_pos[3] = {(2*i+1)*half_sep, (2*j+1)*half_sep, (2*k+1)*half_sep};
						SU2vector phi_hat;
						this->GenerateRandomHiggsComponents(phi_hat);
						this->NucleateOneBubble(T, core_pos, phi_hat);
					}
				}
			}
			this->new_bubbles_count_ = total_bubbles;
			this->phi_.update_halo();
		} else{
			this->new_bubbles_count_ = 0;
		}
		return;
		/*
		if (_time_step == 0) {
			const auto num_in_line = nSize[0] / (2 * half_sep);
			vector<int> core_pos(num_in_line);
			for (int i = 0; i < num_in_line; ++i) {
				core_pos[i] = (2 * i + 1)*half_sep;
			}
			const auto len = core_pos.size();
			new_bubble_count = 0;
			for (int i = 0; i < len; ++i) {
				for (int j = 0; j < len; ++j) {
					for (int k = 0; k < len; ++k) {
						int c[DIM] = { core_pos[i], core_pos[j], core_pos[k] };
						x.setCoord(c);
						NucleateRegion_Exp(x, NUCLEATION_RADIUS_SITE, 1);
						new_bubble_count++;
					}
				}
			}
		}*/
		return;
	}

	template<int DIM>
	NucleationObserver<DIM>::NucleationObserver(const BubbleNucleation<DIM>& bubble)
		:
		nucl_(bubble),
		ElectroweakObserver<DIM>(bubble){

	}

	/*
	template<int DIM>
	void NucleationObserver<DIM>::SetObservables(
			const FlagType data_table_flags,
			const FlagType density_data_flags) {
		this->data_table_flags_ = data_table_flags;
		this->density_data_flags_ = density_data_flags;
		this->init_data_table();
		this->init_density_data();
	}
	*/

	template<int DIM>
	void NucleationObserver<DIM>::Measure() {
		this->basic_measure();
		auto time_step = this->evo_.get_time_step();
		if (this->data_table_flags_ & ObserverFlags::OBS_NewBubbleCount) {
			this->CalcNewBubbleCount(time_step);
		}
		return;
	}
/*
	template<int DIM>
	void NucleationObserver<DIM>::init_data_table(){
		std::vector<std::string> data_table_names;
		this->init_extend_name_vector(this->data_table_flags_, data_table_names);
		this->data_table_.initialize(data_table_names);
		this->data_table_names_ = data_table_names;
	}

	template<int DIM>
	void NucleationObserver<DIM>::init_density_data() {
		std::vector<std::string> density_data_names;
		this->init_extend_name_vector(this->density_data_flags_, density_data_names);
		this->density_names_ = density_data_names;
		this->density_data_.reinit(this->density_names_.size(), 1);
		return;
	}
*/
	template<int DIM>
	void NucleationObserver<DIM>::init_name_vector(
			const FlagType flags,
			std::vector<std::string>& names) const {
		names.clear();
		this->init_basic_name_vector(flags, names);
		if (flags & ObserverFlags::OBS_NewBubbleCount) {
				names.push_back("NewBubbleCount");
		}
		return;
	}

	template<int DIM>
	void NucleationObserver<DIM>::CalcNewBubbleCount(const int time_step){
		auto dt_idx = this->find_index(this->data_table_names_, "NewBubbleCount");
		auto den_idx = this->find_index(this->density_names_, "NewBubbleCount");
		if(den_idx != -1){
			std::cerr << "NewBubbleCount does not have density info." << std::endl;
		}
		if(dt_idx != -1)
			this->data_table_.append_value(dt_idx, this->nucl_.GetNewBubbleCount());
		return;
	}

}

#endif