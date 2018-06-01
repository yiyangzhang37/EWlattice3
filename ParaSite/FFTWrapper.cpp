#include "FFTWrapper.h"
#include <algorithm>
#include <functional>


namespace FFT_Wrapper{

    ComplexFFT1D::ComplexFFT1D(const int size)
        :
        FFTBase<1, complex_t, complex_t>(&size)
        {}

    ComplexFFT1D::~ComplexFFT1D(){
    }

    void ComplexFFT1D::fft(const complex_t* data_in) {
        this->fft_impl(data_in, FFTW_FORWARD);
    }

    void ComplexFFT1D::ifft(const complex_t* data_in) {
        this->fft_impl(data_in, FFTW_BACKWARD);
    }

    void ComplexFFT1D::fft_impl(const complex_t* data_in, const int sign){
        this->set_input_ptr(data_in);
        this->reallocate();
        auto plan = fftw_plan_dft_1d(this->size_[0], 
                                    reinterpret_cast<fftw_complex*>(const_cast<complex_t*>(this->data_in_)),
                                    reinterpret_cast<fftw_complex*>(this->data_out_),
                                    sign,
                                    FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        return;                            
    }

}




