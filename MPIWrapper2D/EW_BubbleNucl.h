#ifndef EW_BUBBLE_NUCL_H
#define EW_BUBBLE_NUCL_H

#include "EW_Model.h"

namespace BubbleNucleation{
    using namespace Electroweak;

    /*The probability that one site can nucleate*/
	const Real NUCLEATION_PROB = std::stod(ReadConfigIni("NULEATION_PROBABILITY"));
	/*Higgs magnitude larger than this value cannot nucleate*/
	const Real NUCLEATION_LIMIT = 0.0001*v2;
	/*Hard cutoff for the nucleation area, should be > 3 * HIGGS_STD*/
	const unsigned int NUCLEATION_RADIUS_SITE = std::stoi(ReadConfigIni("BUBBLE_LAT_RADIUS"));
	const Real NUCLEATION_RADIUS_LEN = NUCLEATION_RADIUS_SITE*DX;
	//Higgs_std should be 1/3 of nucleation radius, since we assume the radius is 3*sigma.
	const Real HIGGS_STD_LEN = NUCLEATION_RADIUS_LEN / 3.0; //given \eta=6.0
	
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
		//The random procedure should be conducted only on one process.
        void RandomBubbleNucleation(unsigned int& new_bubble_count);
        bool NucleateRegionTest_Sphere(const Site<DIM>& x_,
			const int regionHalfSize,
			const int nowTime,
			const double phi2_limit);
        
		void NucleateRegion_Exp(const Site<DIM>& x_, const int regionHalfSize, const int nowTime);

        void OneBubbleTest();
        void TwoBubbleTest(const int sep);
        void NonRandomTest();

        int total_bubbles_count_ = 0;
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
	void BubbleNucleation<DIM>::RandomBubbleNucleation(unsigned int& new_bubble_count) {
		// I feel the nucleation should be more natural to be after update fields.
		// and the nucleation is performed on the updated fields.
		const auto T = (this->time_step_ + 1) % CYCLE;
		new_bubble_count = 0;
		auto seed = std::chrono::system_clock::now().time_since_epoch().count();
		// initialize an ordered site index vector.
		std::vector<IndexType> all_sites(nSize[0] * nSize[1] * nSize[2]);
		std::generate(all_sites.begin(), all_sites.end(), [n = 0]() mutable { return n++; });
		//random shuffle
		std::random_shuffle(all_sites.begin(), all_sites.end());

		//pick some positions with NUCLEATION_PROB
		std::default_random_engine generator(seed);
		std::uniform_real_distribution<Real> distrib(0.0, 1.0);
		std::vector<IndexType> picked_positions;
		picked_positions.reserve(static_cast<int>(NUCLEATION_PROB*nSize[0] * nSize[1] * nSize[2]));
		for (auto i = 0; i < nSize[0] * nSize[1] * nSize[2]; ++i) {
			const auto site_nucl_p = distrib(generator);
			if (site_nucl_p <= NUCLEATION_PROB) picked_positions.push_back(all_sites[i]);
		}

		for (auto i : picked_positions) {
			x.setIndex(i);
			if (NucleateRegionTest_Sphere(x, NUCLEATION_RADIUS_SITE, T, NUCLEATION_LIMIT)) {
				NucleateRegion_Exp(x, NUCLEATION_RADIUS_SITE, T);
				new_bubble_count++;
			}
		}
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