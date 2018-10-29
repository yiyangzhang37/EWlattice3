#ifndef EW_CSBUBBLE_H
#define EW_CSBUBBLE_H

#include <cmath>
#include "EW_BubbleNucl.h"
#include "EW_Transformation.h"

namespace EW_BubbleNucleation {
    using namespace Electroweak;

    /*Radius for the winding-changing profile function*/
	const Real NUCLEATION_CS_RADIUS = NUCLEATION_RADIUS_LEN * 0.25;

    template<int DIM>
    class CSBubble : public BubbleNucleation<DIM> {
    public:
        CSBubble(const Lattice<DIM>& lat, 
                const Parallel2D& parallel, 
                const std::string& id);
        virtual ~CSBubble() {}

        void OneBubbleTest_WithWinding(const int winding, const SU2vector& phi_hat) const;
        void TwoBubblesTest_WithWinding(const int winding1, const int winding2, 
            const SU2vector& phi1_hat, const SU2vector& phi2_hat) const;
        //Two bubble collision, similar with the above function.
        //But with some small purturbations inside one bubble.
        void TwoBubblesTest_WithWinding_Perturbed(const int winding1, const int winding2, 
            const SU2vector& phi1_hat, const SU2vector& phi2_hat) const;
        //Each pair of bubbles are supposed to be located within a 100*100*200 lattice
        //If one wants to make the whole configuration symmetric to each pair,
        //the size of the whole lattice should be a multiple of this size.
        void ArrayOfTwoBubblesTest_WithWinding(const int winding1, const int winding2, 
            const SU2vector& phi1_hat, const SU2vector& phi2_hat) const;
        void InitPureGauge(const int winding, const double r_scale) const;
    protected:
        /*
        Nucleate a bubble with the exponential profile,
        the Higgs field is uniform inside the bubble as in BubbleNucleation class.
        Then a winding-changing gauge transform is performed on both the Higgs
        and the SU(2) gauge field.
        */
        void NucleateOneBubble_Exp_WithWinding(
            const int nowTime,
            const IndexType global_index,
            const SU2vector& phi_hat,
            const HedgehogWinding& winding) const;
        
        //void SetLocalizedPureGaugeSU2Field(
        //    const int nowTime,
        //    const IndexType global_index,
        //    const double lat_radius,
        //    const HedgehogWinding& winding) const;

    };

    template<int DIM>
    CSBubble<DIM>::CSBubble(
        const Lattice<DIM>& lat, 
        const Parallel2D& parallel, 
        const std::string& id)
        :
        BubbleNucleation<DIM>(lat, parallel, id)
    {}

    //I have tested that: 
    //When the bubbls is either at center, surface, edge, or corner,
    //the evolutions (I did up to 5 steps), yield the same results.
    template<int DIM>
    void CSBubble<DIM>::OneBubbleTest_WithWinding(const int winding, const SU2vector& phi_hat) const {
        if (this->time_step_ == 0) {
            IndexType global_coord[] = {nSize[0] / 2, nSize[1] / 2, nSize[2] / 2};
            auto global_idx = this->get_lattice().global_coord2index(global_coord);
            double center_coord[DIM];
            std::transform(global_coord, global_coord + DIM, 
                        CENTER_POS, center_coord, 
                        [](IndexType x, Real c){return (x-c)*DX;});
            double periods[DIM] = {nSize[0]*DX, nSize[1]*DX, nSize[2]*DX};
            HedgehogWinding w(center_coord, winding, NUCLEATION_CS_RADIUS, periods);
            Site<DIM> x(this->lat_);
		
			auto T = (this->time_step_ + 1) % CYCLE;
			this->NucleateOneBubble_Exp_WithWinding(T, global_idx, phi_hat, w);
			this->phi_.update_halo();
            this->U_.update_halo();
		}
		return;
        
    }

    template<int DIM>
    void CSBubble<DIM>::TwoBubblesTest_WithWinding(const int winding1, const int winding2, const SU2vector& phi1_hat, const SU2vector& phi2_hat) const {
        if (this->time_step_ == 0) {
            IndexType c1[] = {nSize[0] / 2, nSize[1] / 2, nSize[2] / 2 - BUBBLES_HALF_SEP};
            IndexType c2[] = {nSize[0] / 2, nSize[1] / 2, nSize[2] / 2 + BUBBLES_HALF_SEP};
            double center1[DIM], center2[DIM];
            std::transform(c1, c1+DIM, CENTER_POS, 
                        center1, 
                        [](IndexType x, Real c){return (x-c)*DX;});
            std::transform(c2, c2+DIM, CENTER_POS, 
                        center2, 
                        [](IndexType x, Real c){return (x-c)*DX;});
            double periods[DIM] = {nSize[0]*DX, nSize[1]*DX, nSize[2]*DX};
            HedgehogWinding w1(center1, winding1, NUCLEATION_CS_RADIUS, periods);
            HedgehogWinding w2(center2, winding2, NUCLEATION_CS_RADIUS, periods);
            Site<DIM> x(this->lat_);
			auto T = (this->time_step_ + 1) % CYCLE;
            auto global_idx_1 = this->get_lattice().global_coord2index(c1);
            auto global_idx_2 = this->get_lattice().global_coord2index(c2);
			this->NucleateOneBubble_Exp_WithWinding(T, global_idx_1, phi1_hat, w1);
            this->NucleateOneBubble_Exp_WithWinding(T, global_idx_2, phi2_hat, w2);
			this->phi_.update_halo();
            this->U_.update_halo();
		}
        return;
    }

    template<int DIM>
    void CSBubble<DIM>::TwoBubblesTest_WithWinding_Perturbed(const int winding1, const int winding2, const SU2vector& phi1_hat, const SU2vector& phi2_hat) const {
        if (this->time_step_ == 0) {
            IndexType c1[] = {nSize[0] / 2, nSize[1] / 2, nSize[2] / 2 - BUBBLES_HALF_SEP};
            IndexType c2[] = {nSize[0] / 2, nSize[1] / 2, nSize[2] / 2 + BUBBLES_HALF_SEP};
            double center1[DIM], center2[DIM];
            std::transform(c1, c1+DIM, CENTER_POS, 
                        center1, 
                        [](IndexType x, Real c){return (x-c)*DX;});
            std::transform(c2, c2+DIM, CENTER_POS, 
                        center2, 
                        [](IndexType x, Real c){return (x-c)*DX;});
            double periods[DIM] = {nSize[0]*DX, nSize[1]*DX, nSize[2]*DX};
            HedgehogWinding w1(center1, winding1, NUCLEATION_CS_RADIUS, periods);
            HedgehogWinding w2(center2, winding2, NUCLEATION_CS_RADIUS, periods);
            
			auto T = (this->time_step_ + 1) % CYCLE;
            auto global_idx_1 = this->get_lattice().global_coord2index(c1);
            auto global_idx_2 = this->get_lattice().global_coord2index(c2);
			this->NucleateOneBubble_Exp_WithWinding(T, global_idx_1, phi1_hat, w1);
            this->NucleateOneBubble_Exp_WithWinding(T, global_idx_2, phi2_hat, w2);

            //set a purturbation inside one bubble
            IndexType c_pt[] = {nSize[0] / 2 - 2, nSize[1] / 2, nSize[2] / 2 - BUBBLES_HALF_SEP};
            IndexType c_local[] = {0,0,0};
            IndexType c_mem[] = {0,0,0};
            this->lat_.global_coord_to_local_vis_coord(c_pt, c_local);
            this->lat_.local_vis_coord_to_local_mem_coord(c_local, c_mem);
            Site<DIM> x(this->lat_);
            x.set_index(this->lat_.local_mem_coord2index(c_mem));
            auto phi_ori = this->phi_(x, T); //unpurturbed field.
            auto phi_mag = phi_ori.norm();
            SU2vector phi_pt;
            phi_pt(0) = Cmplx(0.141067, 0);
            phi_pt(1) = Cmplx(0.99, 0);
            this->phi_(x, T) = phi_mag *phi_pt;

			this->phi_.update_halo();
            this->U_.update_halo();
		}
        return;
    }

    template<int DIM>
    void CSBubble<DIM>::ArrayOfTwoBubblesTest_WithWinding(const int winding1, const int winding2, 
            const SU2vector& phi1_hat, const SU2vector& phi2_hat) const {
        if (this->time_step_ == 0) {
            Site<DIM> x(this->lat_);
			auto T = (this->time_step_ + 1) % CYCLE;
            for(int i = 0; i < nSize[0]/100; i++){
                for(int j = 0; j < nSize[1]/100; j++){
                    for(int k = 0; k < nSize[2]/200; k++){
                        IndexType c1[] = {50 + 100 * i, 50 + 100 * j, 100 + 200 * k - BUBBLES_HALF_SEP};
                        IndexType c2[] = {50 + 100 * i, 50 + 100 * j, 100 + 200 * k + BUBBLES_HALF_SEP};
                        double center1[DIM], center2[DIM];
                        std::transform(c1, c1+DIM, CENTER_POS, 
                                    center1, 
                                    [](IndexType x, Real c){return (x-c)*DX;});
                        std::transform(c2, c2+DIM, CENTER_POS, 
                                    center2, 
                                    [](IndexType x, Real c){return (x-c)*DX;});
                        double periods[DIM] = {nSize[0]*DX, nSize[1]*DX, nSize[2]*DX};
                        HedgehogWinding w1(center1, winding1, NUCLEATION_CS_RADIUS, periods);
                        HedgehogWinding w2(center2, winding2, NUCLEATION_CS_RADIUS, periods);
                        auto global_idx_1 = this->get_lattice().global_coord2index(c1);
                        auto global_idx_2 = this->get_lattice().global_coord2index(c2);
                        this->NucleateOneBubble_Exp_WithWinding(T, global_idx_1, phi1_hat, w1);
                        this->NucleateOneBubble_Exp_WithWinding(T, global_idx_2, phi2_hat, w2);
                    }
                }
            }
			this->phi_.update_halo();
            this->U_.update_halo();
		}
        return;
    }

    template<int DIM>
    void CSBubble<DIM>::NucleateOneBubble_Exp_WithWinding(
        const int nowTime,
        const IndexType global_index,
        const SU2vector& phi_hat,
        const HedgehogWinding& winding) const {
        auto w = winding;
        //First of all, the bubble is still limited within 
        //NUCLEATION_RADIUS_SITE lattice distance.

        std::vector<IndexType> region_list;
		this->GetBubbleRegion(global_index, NUCLEATION_RADIUS_SITE, region_list);
        // std::cout << region_list.size() << std::endl;
        // for(auto x : region_list) {
        //     IndexType coord[DIM];
        //     this->get_lattice().global_index2coord(x, coord);
        //     std::cout << "[" << coord[0] << ", " << coord[1] << ", " << coord[2] << "]" << std::endl;
        // }
        IndexType global_center_coord[DIM];
        this->get_lattice().global_index2coord(global_index, global_center_coord);

        // set the Higgs field
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
                
                //Higgs field is transformed.
                Real coord[] = {rx(x, 0), rx(x, 1), rx(x, 2)};
                w.set_location(coord);
				this->phi_(x, nowTime) = mag * v * w.gauge_transform(phi_hat);
			} else continue;
		}

        //Set the SU2 field.
        //x-component
        this->GetLinkFieldRegion(global_index, 0,
                                NUCLEATION_RADIUS_SITE, region_list);
        for(auto gidx : region_list){
			//check if gidx is a local visible site
			if(this->lat_.is_local(gidx)){
                auto vis_idx = this->lat_.global_index_to_local_vis_index(gidx);
				auto mem_idx = this->lat_.local_vis_index_to_local_mem_index(vis_idx);
                x.set_index(mem_idx);
                Real gc_x[] = {rgx(x, 0), rx(x, 1), rx(x, 2)};
                w.set_location(gc_x);
                this->U_(x, 0, nowTime) = w.pure_gauge(0, DX);
            } else continue;
        }
        //y-component
        this->GetLinkFieldRegion(global_index, 1,
                                NUCLEATION_RADIUS_SITE, region_list);
        for(auto gidx : region_list){
			//check if gidx is a local visible site
			if(this->lat_.is_local(gidx)){
                auto vis_idx = this->lat_.global_index_to_local_vis_index(gidx);
				auto mem_idx = this->lat_.local_vis_index_to_local_mem_index(vis_idx);
                x.set_index(mem_idx);
                Real gc_y[] = {rx(x, 0), rgx(x, 1), rx(x, 2)};
                w.set_location(gc_y);
                this->U_(x, 1, nowTime) = w.pure_gauge(1, DX);
            } else continue;
        }
        //z-component
        this->GetLinkFieldRegion(global_index, 2,
                                NUCLEATION_RADIUS_SITE, region_list);
        //std::cout << "z-component size: " << region_list.size() << std::endl;
        for(auto gidx : region_list){
			//check if gidx is a local visible site
			if(this->lat_.is_local(gidx)){
                auto vis_idx = this->lat_.global_index_to_local_vis_index(gidx);
				auto mem_idx = this->lat_.local_vis_index_to_local_mem_index(vis_idx);
                x.set_index(mem_idx);
                Real gc_z[] = {rx(x, 0), rx(x, 1), rgx(x, 2)};
                w.set_location(gc_z);
                this->U_(x, 2, nowTime) = w.pure_gauge(2, DX);
            } else continue;
        }
		return;
    }

    template<int DIM>
    void CSBubble<DIM>::InitPureGauge(const int winding, const double r_scale) const {
        IndexType global_coord[] = {nSize[0] / 2, nSize[1] / 2, nSize[2] / 2};
        double center_coord[DIM];
        std::transform(global_coord, global_coord + DIM, 
                    CENTER_POS, center_coord, 
                    [](IndexType x, Real c){return (x-c)*DX;});
        Site<DIM> x(this->lat_);
        double periods[DIM] = {nSize[0]*DX, nSize[1]*DX, nSize[2]*DX};
        HedgehogWinding w(center_coord, winding, r_scale, periods);
        //HedgehogWinding_Tanh2 w(winding, r_scale);
		for (auto t = 0; t < CYCLE; ++t) {
			for (x.first(); x.test(); x.next()) {
				this->phi_(x, t) = SU2vector(0, 0);
				this->pi_(x, t) = SU2vector(0, 0);
				for (auto i = 0; i < DIM; ++i) {
					this->F_(x, i, t) = UNITY_F * Ident;
					this->V_(x, i, t) = Cmplx(1, 0);
					this->E_(x, i, t) = UNITY_E * Cmplx(1, 0);
				}
                Real gc_x[] = {rgx(x, 0), rx(x, 1), rx(x, 2)};
                w.set_location(gc_x);
                //SU2matrix gWx = w.pure_gauge(0, dx);
                //this->U_(x, 0, t) = su2_W2U(gWx);
                this->U_(x, 0, t) = w.pure_gauge(0, DX);
                
                Real gc_y[] = {rx(x, 0), rgx(x, 1), rx(x, 2)};
                w.set_location(gc_y);
                this->U_(x, 1, t) = w.pure_gauge(1, DX);

                Real gc_z[] = {rx(x, 0), rx(x, 1), rgx(x, 2)};
                w.set_location(gc_z);
                this->U_(x, 2, t) = w.pure_gauge(2, DX);
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

}


#endif