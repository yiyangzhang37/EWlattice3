#ifndef EW_BUBBLE_NUCL_H
#define EW_BUBBLE_NUCL_H

#include "EW_Base.h"

namespace BubbleNucleation{
    using namespace Electroweak;

    /*The probability that one site can nucleate*/
	const Real NUCLEATION_PROB = std::stod(ReadConfigIni("NULEATION_PROBABILITY"));
	/*Higgs magnitude larger than this value cannot nucleate*/
	const Real NUCLEATION_LIMIT = 0.0001*v2;
	/*Hard cutoff for the nucleation area, should be > 3 * HIGGS_STD*/
	const unsigned int NUCLEATION_RADIUS_SITE = std::stoi(ReadConfigIni("BUBBLE_LAT_RADIUS"));
	const Real NUCLEATION_RADIUS_LEN = NUCLEATION_RADIUS_SITE*DX;
	
	const int TWO_BUBBLES_HALF_SEP = std::stoi(ReadConfigIni("BUBBLE_HALF_SEPARATION"));
	// time offset due to the starting point of bubble collision at different time 
	const int RECORD_TIME_OFFSET = std::stoi(ReadConfigIni("RECORD_TIME_OFFSET"));

	const int NUCL_REGION_LIMIT_SPPEDUP = 1000000; //11 if NUCLEATION_RADIUS_SITE <= 6

	const int STOP_MEAN_STEPS = 10;

	template<int DIM>
    class BubbleNucleation : public ElectroweakEvolution<DIM> {
    public:
		BubbleNucleation(
			const Lattice<DIM>& lat, 
			const Parallel2D& parallel, 
			const std::string& id);
		BubbleNucleation() = default;
        void run();
		void RecordCustomParameters() override;

    private:
        /*Initial condition*/
        void InitializeSymmetricPhase() const;

		/*nulceation*/
		/*
		The overall control of random bubble nucleation.
		Return: new bubble count.
		*/
        int RandomBubbleNucleation() const;

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
			const std::vector<char>& higgs_phase_info,
			const std::vector<Real>& preselected_positions,
			std::vector<Real>& picked_positions) const;
		
		/*
		Nucleate one bubble with the bubble center specified with global_coord,
		and the randomized Higgs angluar components specified by ha[4].
		*/
		void NucleateOneBubble_Exp(const IndexType global_index, const SU2vector& phi_hat) const;
		
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
		This function should only take place on the root process.
		*/
		void GenerateRandomHiggsDirection(SU2Vector& phi_hat) const;

        void OneBubbleTest();
        void TwoBubbleTest(const int sep);
        void NonRandomTest();

        int total_bubbles_count_ = 0;
	};

	template<int DIM>
	class BubbleObserver : public ElectroweakObserver<DIM>{
	private:
	public:
	};

	template<int DIM>
	BubbleNucleation<DIM>::BubbleNucleation(
		const Lattice<DIM>& lat,
		const Parallel2D& parallel,
		const std::string& id)
		:
		ElectroweakEvolution(lat, parallel, id) {}

	template<int DIM>
	void BubbleNucleation<DIM>::RecordCustomParameters() {
		this->param_.add("NucleationProbability", NUCLEATION_PROB);
		this->param_.add("NucleationLimit", NUCLEATION_LIMIT);
		this->param_.add("BubbleRadius(*DX)", NUCLEATION_RADIUS_SITE);
		this->param_.add("BubbleRadius", NUCLEATION_RADIUS_LEN);
		this->param_.add("BubbleFixedHalfSeparation(ifFixed)", TWO_BUBBLES_HALF_SEP);
		return;
	}

	template<int DIM>
	void BubbleNucleation<DIM>::InitializeSymmetricPhase() const {
		auto& x = this->x_;
		const auto nowTime = this->time_step_ % CYCLE;
		const auto nextTime = (this->time_step_ + 1) % CYCLE;
		for (x.first(); x.test(); x.next()) {
			this->phi_(x, nextTime) = this->phi_(x, nowTime) = SU2vector(0, 0);
			this->pi_(x, nextTime) = this->pi_(x, nowTime) = SU2vector(0, 0);
			for (unsigned int i = 0; i < DIM; ++i) {
				this->U_(x, i, nextTime) = this->U_(x, i, nowTime) = Ident;
				this->F_(x, i, nextTime) = this->F_(x, i, nowTime) = UNITY_F * Ident;
				this->V_(x, i, nextTime) = this->V_(x, i, nowTime) = Cmplx(1, 0);
				this->E_(x, i, nextTime) = this->E_(x, i, nowTime) = UNITY_E * Cmplx(1, 0);
			}
		}
		return;
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
		this->parallel_.Barrier();
		this->parallel_.Broadcast(all_broken_phase, this->parallel_.root());
		if(all_broken_phase) return 0; //All sites in broken phase, return 0.

		this->ReorderHiggsPhaseInfo(higgs_phase_info);

		//select some number of sites by p_B.
		std::vector<IndexType> preselect_positions;
		std::vector<IndexType> picked_positions;
		if(this->parallel_.is_root()){
			// initialize an ordered site index vector.
			std::vector<IndexType> all_sites(global_size);
			std::generate(all_sites.begin(), all_sites.end(), 
							[IndexType n = 0]() mutable { return n++; });
			//random shuffle
			std::random_shuffle(all_sites.begin(), all_sites.end());

			//pick some positions with NUCLEATION_PROB
			auto seed = std::chrono::system_clock::now().time_since_epoch().count();
			std::default_random_engine generator(seed);
			std::uniform_real_distribution<Real> distrib(0.0, 1.0);
			preselect_positions.reserve(static_cast<int>(NUCLEATION_PROB*global_size*1.2));
			for (auto item : all_sizes) {
				auto site_nucl_p = distrib(generator);
				if (site_nucl_p <= NUCLEATION_PROB) 
					preselect_positions.push_back(item);
			}
			//verify these locations with higgs_phase_info.
			this->VerifyPickedPositions(
				higgs_phase_info,
				preslect_positions,
				picked_positions);	
		}
		int picked_size = picked_positions.size();
		this->parallel_.Broadcast(picked_size, this->parallel_.root());
		if(picked_size == 0) return 0; //No new bubbles: return 0.
		this->parallel_.Barrier();

		//send the picked_positions to all the processes.
		this->parallel_.Broadcast(picked_size, this->parallel_.get_root());
		if(! this->parallel_.is_root() ){
			picked_positions.resize(picked_size);
		}
		this->parallel_.Broadcast(picked_positions.data(), 
								picked_size, 
								this->parallel_.get_root());

		//generate Higgs angular components and send to each process.
		//The size will be 4*picked_position_size
		std::vector<SU2vector> phi_hat(picked_position_size);
		if(this->parallel_.is_root()){
			for(auto i = 0; i < picked_position_size; ++i){
				this->GenerateRandomHiggsComponents(phi_hat.data() + i);
			}
		}
		this->parallel_.Broadcast(phi_hat.data(), 
								picked_position_size, 
								this->parallel_.get_root()); //transfer as MPI_BYTE
		
		// each process do nucleation, by information from
		// ha and picked_positions.
		for(auto i = 0; i < picked_size; ++i){
			this->NucleateOneBubble_Exp(picked_positions[i], phi_hat[i]);
		}
		this->phi_.update_halo();
		return;
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
								higgs_phase_info.data(), this->parallel_.root());
		return;
	}
	
	template<int DIM>
	void BubbleNucleation<DIM>::ReorderHiggsPhaseInfo(
		std::vector<char>& higgs_phase_info) const {
		;//todo
	}

	template<int DIM>
	void VerifyPickedPositions(
			const std::vector<char>& higgs_phase_info,
			const std::vector<Real>& preselected_positions,
			std::vector<Real>& picked_positions) const {
		//

	}

	template<int DIM>
	void NucleateOneBubble_Exp(
		const int nowTime,
		const IndexType global_index, 
		const SU2vector& phi_hat) const {
		IndexType center_coord[DIM];
		this->lat_.global_index2coord(global_index, center_coord);

		std::vector<IndexType> region_list;
		this->GetBubbleRegion(global_index, NUCLEATION_RADIUS_SITE, region_list);

		Site<DIM> x(this->lat_);
		for(auto gidx : region_list){
			//check if gidx is a local visible site
			if(this->lat_.is_local(gidx)){
				//radial part
				auto r = this->lat_.global_lat_distance(global_index, gidx); //distance to bubble center
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
	void GetBubbleRegion(
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
	void BubbleNucleation<DIM>::GenerateRandomHiggsComponents(SU2Vector& phi_hat) const{
		auto seed = std::chrono::system_clock::now().time_since_epoch().count();
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
	void BubbleNucleation<DIM>::OneBubbleTest() {
		if (_time_step == 0) {
			int c[DIM] = { nSize[0] / 2 , nSize[1] / 2 , nSize[1] / 2 };
			x.setCoord(c);
			NucleateRegion_Exp(x, NUCLEATION_RADIUS_SITE, 1);
		}
		return;
	}

	template<int DIM>
	void BubbleNucleation<DIM>::TwoBubblesTest(const int half_sep) {
		if (_time_step == 0) {
			int c[DIM] = { nSize[0] / 2 - half_sep , nSize[1] / 2 , nSize[2] / 2 };
			x.setCoord(c);
			NucleateRegion_Exp(x, NUCLEATION_RADIUS_SITE, 1);
			int c2[DIM] = { nSize[0] / 2 + half_sep , nSize[1] / 2 , nSize[2] / 2 };
			x.setCoord(c2);
			NucleateRegion_Exp(x, NUCLEATION_RADIUS_SITE, 1);
		}
		return;
	}

	template<int DIM>
	void BubbleNucleation<DIM>::NonRandomTest(const int half_sep, unsigned int& new_bubble_count) {
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
		}
		return;
	}


}

#endif