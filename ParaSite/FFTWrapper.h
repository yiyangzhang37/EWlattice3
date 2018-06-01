#ifndef FFT_WRAPPER_H
#define FFT_WRAPPER_H

#include <complex>
#include <fftw3.h>
#include <iostream>
#include <numeric>
#include <cassert>


namespace FFT_Wrapper{

    typedef std::complex<double> complex_t;
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
        This class will neither allocate, nor deallocate the memory of data_in.
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
        /*
        re-allocate the memory pointed by data_out_, if data_out_ != nullptr.
        allocate the memory to data_out_, if data_out == nullptr.
        */
        void reallocate();

        void set_input_ptr(const InType* data_in);

        int size_[DIM]; 
        int len_;

        const InType* data_in_ = nullptr;
        OutType* data_out_ = nullptr;

    };



    class ComplexFFT1D : public FFTBase<1, complex_t, complex_t>{
    public:
        ComplexFFT1D(const int size);
        ~ComplexFFT1D();

        void fft(const complex_t* data_in) override;
        void ifft(const complex_t* data_in) override;
    
    private:
        void fft_impl(const complex_t* data_in, const int sign);

    };



    template<int DIM, class InType, class OutType>
    FFTBase<DIM, InType, OutType>::FFTBase(const int* size){
        this->resize(size);
    }

    template<int DIM, class InType, class OutType>
    FFTBase<DIM, InType, OutType>::~FFTBase() {
        this->deallocate();
    }


    template<int DIM, class InType, class OutType>
    void FFTBase<DIM, InType, OutType>::resize(const int* size){
        std::copy_n(size, DIM, this->size_);
        this->len_ = std::accumulate(size, size + DIM, 1, std::multiplies<int>());
        assert(this->len_ > 0);
        return;
    }


    template<int DIM, class InType, class OutType>
    void FFTBase<DIM, InType, OutType>::allocate() {
        if( this->data_out_ != nullptr ){
            std::cerr << "Memory already been allocated." << std::endl;
        }
        this->data_out_ = static_cast<OutType*>( fftw_malloc(sizeof(OutType)* this->len_) );
        return;
    }

    template<int DIM, class InType, class OutType>
    void FFTBase<DIM, InType, OutType>::deallocate() {
        if( this->data_out_ == nullptr ){
            std::cerr << "No memory was allocated." << std::endl;
        }
        fftw_free(this->data_out_);
        this->data_out_ = nullptr;
        return;
    }

    template<int DIM, class InType, class OutType>
    void FFTBase<DIM, InType, OutType>::reallocate() {
        if( this->data_out_ == nullptr ) this->allocate();
        else {
            this->deallocate();
            this->allocate();
        }
        return;
    }

    template<int DIM, class InType, class OutType>
    void FFTBase<DIM, InType, OutType>::set_input_ptr(const InType* data_in){
        if(data_in != nullptr){
        this->data_in_ = data_in;
        } else {
            std::cerr << "Input is a nullptr." << std::endl;
        }
        return;
    }


}


#endif