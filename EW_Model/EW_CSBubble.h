#ifndef EW_CSBUBBLE_H
#define EW_CSBUBBLE_H

#include <cmath>
#include "EW_BubbleNucl.h"
#include "EW_Winding.h"

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
        ~CSBubble() = default;

        void OneBubbleTest_WithWinding(const int winding) const;
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

    };

    template<int DIM>
    CSBubble<DIM>::CSBubble(
        const Lattice<DIM>& lat, 
        const Parallel2D& parallel, 
        const std::string& id)
        :
        BubbleNucleation<DIM>(lat, parallel, id)
    {}

    template<int DIM>
    void CSBubble<DIM>::OneBubbleTest_WithWinding(const int winding) const {
        HedgehogWinding w(winding, NUCLEATION_CS_RADIUS);
        Site<DIM> x(this->lat_);
		if (this->time_step_ == 0) {
			auto T = (this->time_step_ + 1) % CYCLE;
			SU2vector phi_hat;
			phi_hat(0) = Cmplx(0,0);
            phi_hat(1) = Cmplx(1,0);
			IndexType global_coord[] = {nSize[0] / 2, nSize[1] / 2, nSize[2] / 2};
            auto global_idx = this->get_lattice().global_coord2index(global_coord);
			this->NucleateOneBubble_Exp_WithWinding(T, global_idx, phi_hat, w);
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
        /*
        This is quite sloppy, because the gauge field locations are on the links,
        and the region_list for gauge field is not exactly the same as that 
        for the Higgs field.
        */
        std::vector<IndexType> region_list;
		this->GetBubbleRegion(global_index, NUCLEATION_RADIUS_SITE, region_list);
        IndexType global_center_coord[DIM];
        this->get_lattice().global_index2coord(global_index, global_center_coord);
        Real gcc[DIM];
        std::transform(global_center_coord, global_center_coord + DIM,
            gcc, [](IndexType c){return static_cast<Real>(c); });
    
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
                Real coord[] = {rx(x, 0, DX, gcc), rx(x, 1, DX, gcc), rx(x, 2, DX, gcc)};
                w.set_location(coord);
				this->phi_(x, nowTime) = mag * v * w.gauge_transform(phi_hat);
                //Set pure-gauge SU(2) field.
                Real gc_x[] = {rgx(x, 0, DX, gcc), rx(x, 1, DX, gcc), rx(x, 2, DX, gcc)};
                w.set_location(gc_x);
                auto gWx = w.pure_gauge(0);
                Cmplx SU2_Ucomp_x[] = {(iPauli[0]*gWx).trace(), 
                                    (iPauli[1]*gWx).trace(), 
                                    (iPauli[2]*gWx).trace()};
                su2field2U(x, 0, this->U_, nowTime, SU2_Ucomp_x);

                Real gc_y[] = {rx(x, 0, DX, gcc), rgx(x, 1, DX, gcc), rx(x, 2, DX, gcc)};
                w.set_location(gc_y);
                auto gWy = w.pure_gauge(1);
                Cmplx SU2_Ucomp_y[] = {(iPauli[0]*gWy).trace(), 
                                    (iPauli[1]*gWy).trace(), 
                                    (iPauli[2]*gWy).trace()};
                su2field2U(x, 1, this->U_, nowTime, SU2_Ucomp_y);

                Real gc_z[] = {rx(x, 0, DX, gcc), rx(x, 1, DX, gcc), rgx(x, 2, DX, gcc)};
                w.set_location(gc_z);
                auto gWz = w.pure_gauge(2);
                Cmplx SU2_Ucomp_z[] = {(iPauli[0]*gWz).trace(), 
                                    (iPauli[1]*gWz).trace(), 
                                    (iPauli[2]*gWz).trace()};
                su2field2U(x, 2, this->U_, nowTime, SU2_Ucomp_z);
                
			} else continue;
		}
		return;
    }

}


#endif