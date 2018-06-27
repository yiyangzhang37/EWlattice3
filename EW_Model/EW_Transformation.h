#ifndef EW_WINDING_H
#define EW_WINDING_H

#include <array>
#include <numeric>
#include "EW_parameter.h"

namespace Electroweak{

    /*
    I am using the conventions listed in invertedBubble.lyx and SU2TopoDegree.lyx
    */

    inline int levi_civita_3d(const int i, const int j, const int k){
        return (i-j)*(j-k)*(k-i)/2;
    }

    /*
    Base class of a position-dependent SU2 gauge transformation.
    The space dimension is assumed to be 3.
    */
    class SU2GaugeTransform {
    public:
        SU2GaugeTransform(const double* center);
        ~SU2GaugeTransform() = default;

        /*
        set the coordinate to compute the location-dependent transformations.
        Attention should be paid that the Higgs field is on the lattice site,
        while gauge fields are on the links. 
        The locations are different.
        */
        void set_location(const double r, const double* unit_dir);
        void set_location(const double* coord);

        virtual SU2matrix transform_mat() const {return Ident;}
        virtual SU2matrix inv_transform_mat() const {return Ident;}
        virtual SU2matrix d_transform_mat(const int dir) const {return SU2matrix::Zero();}
        virtual SU2matrix d_inv_transform_mat(const int dir) const {return SU2matrix::Zero();}

        //compute gW_i = U d_i(U^-1)
        SU2matrix pure_gauge(const int dir) const;
        //compute the component gW_i^a
        double pure_gauge(const int dir, const int wa) const;

        //compute the Higgs field gauge transform: U*phi.
        SU2vector gauge_transform(const SU2vector& phi) const;
        //compute the SU2 field gauge transform: U*gW_i*U^-1 + U*d_i(U^-1). 
        SU2matrix gauge_transform(const SU2matrix& gw, const int dir) const;

        const double r() const {return this->radius_;}
        const double coord(const int i) const {return this->coord_[i];}
    
    protected:
        // The referece center of the gauge transformation.
        std::array<double, 3> center_;
        // The current relative (to the reference center) position of the measured point.
        std::array<double, 3> coord_ = {};
        // The current relative radius.
        double radius_ = 0;

        void set_radius(){
            this->radius_ = std::sqrt(
                std::inner_product(this->coord_.begin(), this->coord_.end(), 
                this->coord_.begin(), 0.0));
        } 

        //This function does not check if this->radius_ is 0.
        double unit_dir(const int i) const {
            return this->coord_[i] / this->radius_;
        }

    };

    /*
    Gauge transformtion with winding number, with Hedgehog ansatz.
    */
    class HedgehogWinding : public SU2GaugeTransform{
    public:
        HedgehogWinding(const double* center, const int winding, const double r_scale);
        ~HedgehogWinding() = default;

        virtual double profile_f() const;
        virtual double profile_df() const;

        SU2matrix transform_mat() const override;
        SU2matrix inv_transform_mat() const override;
        SU2matrix d_transform_mat(const int dir) const override;
        SU2matrix d_inv_transform_mat(const int dir) const override;

        SU2matrix winding_mat() const {return this->transform_mat();}
        SU2matrix inv_winding_mat() const {return this->inv_transform_mat();}
        SU2matrix d_winding_mat(const int dir) const {return this->d_transform_mat(dir);}
        SU2matrix d_inv_winding_mat(const int dir) const {return this->d_inv_transform_mat(dir);}
        
        //compute gW_i = U d_i(U^-1)
        SU2matrix pure_gauge_direct(const int dir) const;
        //compute the component gW_i^a
        double pure_gauge_direct(const int dir, const int wa) const;

        const int get_winding_number() const {return this->winding_;}

    protected:
        double r_scale_;
        int winding_;
        
        double w(const int wa) const {
            return this->unit_dir(wa);
        }
        double dr(const int i) const {
            return this->unit_dir(i);
        }
        double dw(const int i, const int wa) const {
            return (static_cast<double>(i==wa) - this->unit_dir(wa)*this->unit_dir(i))/this->radius_;
        }

        SU2matrix d_mat_impl(const int dir, const int sign) const;
        SU2matrix winding_mat_impl(const int sign) const;
    
    };

    class HedgehogWinding_Tanh2 : public HedgehogWinding {
    public:
        HedgehogWinding_Tanh2(const double* center, const int winding, const double r_scale)
            :
            HedgehogWinding(center, winding, r_scale) 
        {}
        
        ~HedgehogWinding_Tanh2() = default;

        double profile_f() const override {
            return 2.0 * PI * this->winding_ * tanh( pow(this->radius_/this->r_scale_, 2) );
        }
        double profile_df() const override {
            auto& r = this->radius_;
            auto& rg = this->r_scale_;
            return 4.0*PI*this->winding_ * r / (rg*rg) / pow( cosh(pow(r/rg, 2)), 2 );
        }
    };

}


#endif