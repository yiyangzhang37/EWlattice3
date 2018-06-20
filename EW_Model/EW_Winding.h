#ifndef EW_WINDING_H
#define EW_WINDING_H

#include <cmath>
#include <numeric>
#include <algorithm>
#include <array>
#include "EW_parameter.h"

namespace Electroweak{

    /*
    I am using the conventions listed in invertedBubble.lyx and SU2TopoDegree.lyx
    */

    inline int levi_civita_3d(const int i, const int j, const int k){
        return (i-j)*(j-k)*(k-i)/2;
    }

    template<int W>
    class HedgehogWinding{
    public:
        HedgehogWinding(const double r_scale);
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

        constexpr int get_winding_number() const {return W;}
        const int r() const {return this->radius_;}
        const int coord(const int i) const {return this->coord_[i];}
    
    protected:
        double r_scale_;
        double radius_;
        std::array<double, 3> coord_; //DIM==3 is assumed.

        void set_radius() const{
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

    template<int W>
    HedgehogWinding<W>::HedgehogWinding(const double r_scale)
        :
        r_scale_(r_scale) 
    {}

    template<int W>
    void HedgehogWinding<W>::set_location(const double r, const double* unit_dir) {
        std::array<double, 3> coord = {r*unit_dir[0], r*unit_dir[1], r*unit_dir[2]};
        this->set_location(coord.data());
    }

    template<int W>
    void HedgehogWinding<W>::set_location(const double* coord) {
        std::copy(coord, coord + 3, this->coord_.begin());
        this->set_radius();
    }
    
    template<int W>
    double HedgehogWinding<W>::profile_f() const {
        return 2.0 * PI * W * tanh(this->radius_/this->r_scale_);
    }

    template<int W>
    double HedgehogWinding<W>::profile_df() const {
        return 2.0 * PI * W / (this->r_scale_ * pow(cosh(this->radius_/this->r_scale_), 2));
    }

    template<int W>
    SU2matrix HedgehogWinding<W>::winding_mat() const {
        return cos(0.5*this->profile_f()) * Ident + sin(0.5*this->profile_f()) * 
        (iPauli[0]*this->unit_dir(0) + iPauli[1]*this->unit_dir(1) + iPauli[2]*this->unit_dir(2));
    }

    template<int W>
    SU2matrix HedgehogWinding<W>::inv_winding_mat() const {
        return cos(0.5*this->profile_f()) * Ident - sin(0.5*this->profile_f()) * 
        (iPauli[0]*this->unit_dir(0) + iPauli[1]*this->unit_dir(1) + iPauli[2]*this->unit_dir(2));
    }

    template<int W>
    SU2matrix HedgehogWinding<W>::d_mat_impl(
        const int dir, 
        const int sign) const {    
        auto sinf2 = sin(0.5*this->profile_f();
        auto cosf2 = cos(0.2*this->profile_f());
        auto dfi = this->profile_df() * this->dr(dir);
        return -0.5 * sinf2 * dfi * Ident
            + sign * 
            ( iPauli[0]*(this->dw(dir, 0)*sinf2 + 0.5*this->w(0)*cosf2*dfi)
            + iPauli[1]*(this->dw(dir, 1)*sinf2 + 0.5*this->w(1)*cosf2*dfi)
            + iPauli[2]*(this->dw(dir, 2)*sinf2 + 0.5*this->w(2)*cosf2*dfi) );
    }

    template<int W>
    SU2matrix HedgehogWinding<W>::d_winding_mat(
        const int dir) const {
        this->d_mat_impl(dir, 1);
    }

    template<int W>
    SU2matrix HedgehogWinding<W>::d_inv_winding_mat(
        const int dir) const {
        this->d_mat_impl(dir, -1);
    }

    template<int W>
    SU2matrix HedgehogWinding<W>::pure_gauge(const int dir) const {
        return -0.5*(
            iPauli[0]*this->pure_gauge(dir, 0) +
            iPauli[1]*this->pure_gauge(dir, 1) +
            iPauli[2]*this->pure_gauge(dir, 2) );
    }

    template<int W>
    SU2matrix HedgehogWinding<W>::pure_gauge(const int dir, const int wa) const {
        auto sinfr = sin(this->profile_f()) / this->radius_;
        auto dfr = this->profile_df();
        auto xaxi = this->unit_dir(dir) * this->unit_dir(wa);
        auto result = (dir==wa)*sinfr + (dfr-sinfr)*xaxi;
        if(dir != wa){
            auto cosr2 = (1-cos(this->profile_f())) / (this->radius_*this->radius_);
            auto b = 3-dir-wa; //this only works for DIM==3.
            result += cosr2 * levi_civita_3d(wa, dir, b) * this->coord_[b];
        }
        return result; 
    }

    template<int W>
    SU2vector HedgehogWinding<W>::gauge_transform(const SU2vector& phi) const{
        return this->winding_mat()*phi;
    }

    template<int W>
    SU2vector HedgehogWinding<W>::gauge_transform(const SU2matrix& gw, const int dir) const{
        return this->winding_mat() * gw * this->inv_winding_mat() 
            + this->winding_mat() * this->d_inv_winding_mat(dir);
    }




}


#endif