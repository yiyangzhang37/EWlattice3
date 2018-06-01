#ifndef FFT_WRAPPER_H
#define FFT_WRAPPER_H

#include <complex>
#include <fftw3.h>

namespace FFT_Wrapper{

    typedef complex<double> complex_t;
    typedef double real_t;


    /*
    The multi-dimensional data will take the convention of row-major order,
    that is, the last index changes fastest.
    */
    template<int DIM, class InType, class OutType = InType>
    class FFTBase {
    public:
        FFTBase(const int* size);
        ~FFTBase();

        void resize(const int* size);

        /*
        run an FFT given the input data. The result is stored in the data_out_.
        The input data is assumed to be known, with size given by this->size_.
        */
        virtual void fft(const InType* data_in);

        /*
        run an inverse FFT given the input data. The result is stored in the data_out_.
        */
        virtual void ifft(const InType* data_in);

    protected:
        /*allocate memory to data_out_, according to len_*/
        void allocate();
        /*free the memory pointed to data_out_*/
        void deallocate();
        /*re-allocate the memory pointed by data_out_*/
        void reallocate();

        int[DIM] size_; 
        int len_;

        const InType* data_in_ = nullptr;
        OutType* data_out_ = nullptr;

    };



    class ComplexFFT1D : public FFTBase<1, complex_t, complex_t>{
    public:
        ComplexFFT1D(const int size);
        ~ComplexFFT1D() = default;

    };


}


#endif