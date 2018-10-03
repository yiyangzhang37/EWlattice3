#include <cmath>
#include <numeric>
#include <algorithm>
#include "EW_Transformation.h"

namespace Electroweak {

    SU2GaugeTransform::SU2GaugeTransform(const double* center, const double* periods){
        std::copy_n(center, 3, this->center_.begin());
        std::copy_n(periods, 3, this->periods_.begin());
    }

    void SU2GaugeTransform::set_location(const double r, const double* unit_dir) {
        std::array<double, 3> coord = {{r*unit_dir[0], r*unit_dir[1], r*unit_dir[2]}};
        this->set_location(coord.data());
    }

    void SU2GaugeTransform::set_location(const double* coord) {
        std::transform(coord, coord + 3, this->center_.begin(), 
                        this->coord_.begin(), 
                        [](double x, double c){return x-c;});
        for(int i = 0; i < 3; ++i){
            auto d1 = coord[i] + this->periods_[i] - this->center_[i];
            auto d2 = coord[i] - this->center_[i];
            auto d3 = coord[i] - this->periods_[i] - this->center_[i];
            auto absd1 = std::abs(d1);
            auto absd2 = std::abs(d2);
            auto absd3 = std::abs(d3);
            if(absd1 <= absd2 && absd1 <= absd3) this->coord_[i] = d1;
            else if(absd2 <= absd1 && absd2 <= absd3) this->coord_[i] = d2;
            else if(absd3 <= absd1 && absd3 <= absd2) this->coord_[i] = d3;
            else std::cout << "SU2GaugeTransform::set_location ERROR." << std::endl;
        }
    }

    SU2matrix SU2GaugeTransform::pure_gauge(const int dir, const double dx) {
        return this->gauge_transform(SU2matrix::Identity(), dir, dx);
    }

    SU2vector SU2GaugeTransform::gauge_transform(const SU2vector& phi) const {
        return this->transform_mat()*phi;
    }

    SU2matrix SU2GaugeTransform::gauge_transform(const SU2matrix& su2_link, const int dir, const double dx) {
        std::array<double, 3> c0 = this->get_coords();
        std::array<double, 3> center = this->get_center();
        std::transform(c0.begin(), c0.end(), center.begin(), c0.begin(), 
                        [](double x, double y){return x+y;});
        
        std::array<double, 3> c1 {{c0[0] - 0.5*dx*(dir==0), 
                                c0[1] - 0.5*dx*(dir==1), 
                                c0[2] - 0.5*dx*(dir==2)}};
        this->set_location(c1.data());
        auto mat1 = this->transform_mat();

        std::array<double, 3> c2 {{c0[0] + 0.5*dx*(dir==0), 
                                c0[1] + 0.5*dx*(dir==1), 
                                c0[2] + 0.5*dx*(dir==2)}};
        this->set_location(c2.data());
        auto mat2 = this->transform_mat();

        this->set_location(c0.data());
        return mat1 * su2_link * mat2.adjoint(); 
    }

    
    HedgehogWinding::HedgehogWinding(
        const double* center, 
        const int winding, 
        const double r_scale,
        const double* periods)
        :
        SU2GaugeTransform(center, periods),
        r_scale_(r_scale),
        winding_(winding)
    {}
 
    double HedgehogWinding::profile_f() const {
        return 2.0 * PI * this->winding_ * tanh(this->radius()/this->r_scale_);
    }

    double HedgehogWinding::profile_df() const {
        return 2.0 * PI * this->winding_ / (this->r_scale_ * pow(cosh(this->radius()/this->r_scale_), 2));
    }
    
    SU2matrix HedgehogWinding::winding_mat_impl(const int sign) const {
        if(this->radius() < 1e-6) return Ident;
        return cos(0.5*this->profile_f()) * Ident + sign * sin(0.5*this->profile_f())  / this->radius() * 
        (iPauli[0]*this->coord_[0] + iPauli[1]*this->coord_[1] + iPauli[2]*this->coord_[2]);
    }
    SU2matrix HedgehogWinding::transform_mat() const {
        return this->winding_mat_impl(1);
    }

    SU2matrix HedgehogWinding::inv_transform_mat() const {
        return this->winding_mat_impl(-1);
    }

    SU2matrix HedgehogWinding::d_mat_impl(
        const int dir, 
        const int sign) const {
        if(this->radius() < 1e-6) return sign*0.5*this->profile_df()*iPauli[dir];   
        auto sinf2 = sin(0.5*this->profile_f());
        auto cosf2 = cos(0.5*this->profile_f());
        auto dfi = this->profile_df() * this->dr(dir);
        return -0.5 * sinf2 * dfi * Ident
            + sign * 
            ( iPauli[0]*(this->dw(dir, 0)*sinf2 + 0.5*this->w(0)*cosf2*dfi)
            + iPauli[1]*(this->dw(dir, 1)*sinf2 + 0.5*this->w(1)*cosf2*dfi)
            + iPauli[2]*(this->dw(dir, 2)*sinf2 + 0.5*this->w(2)*cosf2*dfi) );
    }

    SU2matrix HedgehogWinding::d_transform_mat(
        const int dir) const {
        return this->d_mat_impl(dir, 1);
    }

    SU2matrix HedgehogWinding::d_inv_transform_mat(
        const int dir) const {
        return this->d_mat_impl(dir, -1);
    }

    SU2matrix HedgehogWinding::pure_gauge_direct(const int dir) const {
        return -0.5*(
            iPauli[0]*this->pure_gauge_direct(dir, 0) +
            iPauli[1]*this->pure_gauge_direct(dir, 1) +
            iPauli[2]*this->pure_gauge_direct(dir, 2) );
    }

    double HedgehogWinding::pure_gauge_direct(const int dir, const int wa) const {
        auto radius = this->radius();
        if(radius < 1e-6){
            if(dir==wa) return this->profile_df();
            else return 0.0;
        } else{
            auto sinfr = sin(this->profile_f()) / radius;
            auto dfr = this->profile_df();
            auto xaxi = this->unit_dir(dir) * this->unit_dir(wa);
            auto result = static_cast<double>(dir==wa)*sinfr + (dfr-sinfr)*xaxi;
            if(dir != wa){
                auto cosr2 = (1.0-cos(this->profile_f())) / (radius*radius);
                auto b = 3-dir-wa; //this only works for DIM==3.
                result += cosr2 * levi_civita_3d(wa, dir, b) * this->coord_[b];
            }
            return result; 
        }
    }
}