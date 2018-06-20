#ifndef EW_WINDING_H
#define EW_WINDING_H

#include <array>
#include "EW_parameter.h"

namespace Electroweak{

    /*
    I am using the conventions listed in invertedBubble.lyx and SU2TopoDegree.lyx
    */

    inline int levi_civita_3d(const int i, const int j, const int k){
        return (i-j)*(j-k)*(k-i)/2;
    }

    class HedgehogWinding{
    public:
        HedgehogWinding(const int winding, const double r_scale);
        ~HedgehogWinding() = default;

        /*
        set the coordinate to compute the location-dependent transformations.
        Attention should be paid that the Higgs field is on the lattice site,
        while gauge fields are on the links. 
        The locations are different.
        */
        void set_location(const double r, const double* unit_dir);
        void set_location(const double* coord);

        virtual double profile_f() const;
        virtual double profile_df() const;


        SU2matrix winding_mat() const;
        SU2matrix inv_winding_mat() const;
        SU2matrix d_winding_mat(const int dir) const;
        SU2matrix d_inv_winding_mat(const int dir) const;
        
        //compute gW_i = U d_i(U^-1)
        SU2matrix pure_gauge(const int dir) const;
        //compute the component gW_i^a
        double pure_gauge(const int dir, const int wa) const;

        //compute the Higgs field gauge transform: U*phi.
        SU2vector gauge_transform(const SU2vector& phi) const;
        //compute the SU2 field gauge transform: U*gW_i*U^-1 + U*d_i(U^-1). 
        SU2matrix gauge_transform(const SU2matrix& gw, const int dir) const;

        constexpr int get_winding_number() const {return this->winding_;}
        const int r() const {return this->radius_;}
        const int coord(const int i) const {return this->coord_[i];}
    
    protected:
        double r_scale_;
        double radius_;
        std::array<double, 3> coord_; //DIM==3 is assumed.
        int winding_;

        void set_radius(){
            this->radius_ = std::sqrt(
                std::inner_product(this->coord_.begin(), this->coord_.end(), 
                this->coord_.begin(), 0));
        }

        double unit_dir(const int i) const {
            return this->coord_[i] / this->r();
        }
        double w(const int wa) const {
            return this->unit_dir(wa);
        }
        double dr(const int i) const {
            return this->unit_dir(i);
        }
        double dw(const int i, const int wa) const {
            return (static_cast<double>(i==wa) - this->unit_dir(wa)*this->unit_dir(i))/this->r();
        }

        SU2matrix d_mat_impl(const int dir, const int sign) const;
    
    };

}


#endif