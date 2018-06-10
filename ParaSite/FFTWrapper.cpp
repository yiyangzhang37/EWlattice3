#include "FFTWrapper.h"
#include <functional>


namespace FFT_Wrapper{
/*
    FFT1D::FFT1D(
        const int size,
        const bool copy_input_data,
        const unsigned fft_flags)
        :
        FFTBase<1, complex_t, complex_t>(&size, copy_input_data, fft_flags)
        {}

    FFT1D::~FFT1D(){
    }

    void FFT1D::fft(const complex_t* data_in) {
        this->fft_impl(data_in, FFTW_FORWARD);
    }

    void FFT1D::ifft(const complex_t* data_in) {
        this->fft_impl(data_in, FFTW_BACKWARD);
        real_t norm = this->get_output_len();
        std::transform(this->data_out_, this->data_out_ + this->get_output_len(), 
                        this->data_out_,
                        [norm](const complex_t& x) { return x / norm; });
    }

    void FFT1D::fft_impl(const complex_t* data_in, const int sign){
        this->init_input_data(data_in);
        this->reallocate_output_data();
        auto plan = fftw_plan_dft_1d(this->get_input_size()[0], 
                                    reinterpret_cast<fftw_complex*>(this->data_in_),
                                    reinterpret_cast<fftw_complex*>(this->data_out_),
                                    sign,
                                    this->get_fft_flags());
        this->import_input_data(data_in);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        return;                            
    }


    FFT1Dr2c::FFT1Dr2c(
            const int input_size,
            const bool copy_input_data,
            const unsigned fft_flags)
        : FFTBase<1, real_t, complex_t>(&input_size, copy_input_data, fft_flags)
        {   
            auto output_size = input_size / 2 + 1;
            this->resize(nullptr, &output_size);
        }
    
    FFT1Dr2c::~FFT1Dr2c(){}
    
    void FFT1Dr2c::fft(const real_t* data_in) {
        this->fft_impl(data_in, FFTW_FORWARD);
        return;
    }

    void FFT1Dr2c::fft_impl(const real_t* data_in, const int sign) {
        this->init_input_data(data_in);
        this->reallocate_output_data();
        auto plan = fftw_plan_dft_r2c_1d(this->get_input_len(), 
                                        this->data_in_, 
                                        reinterpret_cast<fftw_complex*>(this->data_out_),
                                        this->get_fft_flags());
        this->import_input_data(data_in);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        return;
    }


    FFT1Dc2r::FFT1Dc2r(
        const int output_size,
        const bool copy_input_data,
        const unsigned fft_flags) 
        : FFTBase<1, complex_t, real_t>(&output_size, copy_input_data, fft_flags)
    {
        auto input_size = output_size / 2 + 1;
        this->resize(&input_size, nullptr);
    }

    FFT1Dc2r::~FFT1Dc2r(){}

    void FFT1Dc2r::ifft(const complex_t* data_in) {
        this->fft_impl(data_in, FFTW_BACKWARD);
        real_t norm = this->get_output_len();
        std::transform(this->data_out_, this->data_out_ + this->get_output_len(), 
                        this->data_out_,
                        [norm](const real_t x) { return x / norm; });
    }

    void FFT1Dc2r::fft_impl(const complex_t* data_in, const int sign) {
        this->init_input_data(data_in);
        this->reallocate_output_data();
        auto plan = fftw_plan_dft_c2r_1d(this->get_output_len(), 
                                        reinterpret_cast<fftw_complex*>(this->data_in_), 
                                        this->data_out_,
                                        this->get_fft_flags());
        this->import_input_data(data_in);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        return;
    }
*/
}




