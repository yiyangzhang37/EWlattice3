#include <cmath>
#include <numeric>
#include <algorithm>
#include "EW_Winding.h"

namespace Electroweak {
    
    HedgehogWinding::HedgehogWinding(const int winding, const double r_scale)
        :
        r_scale_(r_scale),
        winding_(winding)
    {}

    void HedgehogWinding::set_location(const double r, const double* unit_dir) {
        std::array<double, 3> coord = {{r*unit_dir[0], r*unit_dir[1], r*unit_dir[2]}};
        this->set_location(coord.data());
    }

    void HedgehogWinding::set_location(const double* coord) {
        std::copy(coord, coord + 3, this->coord_.begin());
        this->set_radius();
    }
    
    double HedgehogWinding::profile_f() const {
        return 2.0 * PI * this->winding_ * tanh(this->radius_/this->r_scale_);
    }

    double HedgehogWinding::profile_df() const {
        return 2.0 * PI * this->winding_ / (this->r_scale_ * pow(cosh(this->radius_/this->r_scale_), 2));
    }
    
    SU2matrix HedgehogWinding::winding_mat_impl(const int sign) const {
        if(this->radius_ < 1e-6) return Ident;
        return cos(0.5*this->profile_f()) * Ident + sign * sin(0.5*this->profile_f())  / this->radius_ * 
        (iPauli[0]*this->coord_[0] + iPauli[1]*this->coord_[1] + iPauli[2]*this->coord_[2]);
    }
    SU2matrix HedgehogWinding::winding_mat() const {
        return this->winding_mat_impl(1);
    }

    SU2matrix HedgehogWinding::inv_winding_mat() const {
        return this->winding_mat_impl(-1);
    }

    SU2matrix HedgehogWinding::d_mat_impl(
        const int dir, 
        const int sign) const {
        if(this->radius_ < 1e-6) return sign*0.5*this->profile_df()*iPauli[dir];   
        auto sinf2 = sin(0.5*this->profile_f());
        auto cosf2 = cos(0.5*this->profile_f());
        auto dfi = this->profile_df() * this->dr(dir);
        return -0.5 * sinf2 * dfi * Ident
            + sign * 
            ( iPauli[0]*(this->dw(dir, 0)*sinf2 + 0.5*this->w(0)*cosf2*dfi)
            + iPauli[1]*(this->dw(dir, 1)*sinf2 + 0.5*this->w(1)*cosf2*dfi)
            + iPauli[2]*(this->dw(dir, 2)*sinf2 + 0.5*this->w(2)*cosf2*dfi) );
    }

    SU2matrix HedgehogWinding::d_winding_mat(
        const int dir) const {
        return this->d_mat_impl(dir, 1);
    }

    SU2matrix HedgehogWinding::d_inv_winding_mat(
        const int dir) const {
        return this->d_mat_impl(dir, -1);
    }

    SU2matrix HedgehogWinding::pure_gauge(const int dir) const {
        return -0.5*(
            iPauli[0]*this->pure_gauge(dir, 0) +
            iPauli[1]*this->pure_gauge(dir, 1) +
            iPauli[2]*this->pure_gauge(dir, 2) );
    }

    double HedgehogWinding::pure_gauge(const int dir, const int wa) const {
        if(this->radius_ < 1e-6){
            if(dir==wa) return this->profile_df();
            else return 0.0;
        } else{
            auto sinfr = sin(this->profile_f()) / this->radius_;
            auto dfr = this->profile_df();
            auto xaxi = this->unit_dir(dir) * this->unit_dir(wa);
            auto result = static_cast<double>(dir==wa)*sinfr + (dfr-sinfr)*xaxi;
            if(dir != wa){
                auto cosr2 = (1.0-cos(this->profile_f())) / (this->radius_*this->radius_);
                auto b = 3-dir-wa; //this only works for DIM==3.
                result += cosr2 * levi_civita_3d(wa, dir, b) * this->coord_[b];
            }
            return result; 
        }
    }

    SU2vector HedgehogWinding::gauge_transform(const SU2vector& phi) const{
        return this->winding_mat()*phi;
    }

    SU2matrix HedgehogWinding::gauge_transform(const SU2matrix& gw, const int dir) const{
        return this->winding_mat() * gw * this->inv_winding_mat() 
            + this->winding_mat() * this->d_inv_winding_mat(dir);
    }
}