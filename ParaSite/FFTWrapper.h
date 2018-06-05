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

    /*
    For the output with size n,
    The k-th element corresponds to the frequency k/n;
    or in terms of positive and negative frequencies, the first-half stores the positive frequencies,
    and the second-half stores the negative frequencies 
    (the frequencies -k/n is the same as (n-k)/n).
    Specifically, the 0-th element stores the zero-frequency component.
    */

    /*
    The normalization factor 1/n, is put in the ifft() scheme.
    */
    template<int DIM, class InType, class OutType = InType>
    class FFTBase {
    public:
        FFTBase(
            const int* input_size,
            const bool copy_input_data = true,
            const unsigned fft_flags = FFTW_ESTIMATE);
        FFTBase(const int* input_size,
            const int* output_size,
            const bool copy_input_data = true,
            const unsigned fft_flags = FFTW_ESTIMATE);
        ~FFTBase();

        void resize(
            const int* input_size = nullptr, 
            const int* output_size = nullptr);

        /*
        run an FFT given the input data. The result is stored in the data_out_.
        The input data is assumed to be known, with size given by this->size_.
        This class will neither allocate, nor deallocate the memory of data_in.
        */
        virtual void fft(const InType* data_in) {;}

        /*
        run an inverse FFT given the input data. The result is stored in the data_out_.
        */
        virtual void ifft(const InType* data_in) {;}
        
        const int* get_input_size() const { return this->input_size_; }
        const int* get_output_size() const { return this->output_size_; }
        const int get_input_len() const { return this->input_len_; }
        const int get_output_len() const { return this->output_len_; }
        const InType* get_input_data() const { return this->data_in_; }
        const OutType* get_output_data() const { return this->data_out_; }
        const bool is_input_data_copied() const { return this->is_copy_input_data_; }
        const unsigned get_fft_flags() const { return this->fft_flags_; }

    protected:
        /*allocate memory to ptr, with size alloc_size*/
        template<class Type>
        void allocate(Type*& ptr, const int alloc_size);
        /*free the memory pointed from ptr*/
        template<class Type>
        void deallocate(Type*& ptr);
        /*
        re-allocate the memory pointed from ptr, if ptr != nullptr.
        allocate the memory to ptr, if ptr == nullptr.
        */
        template<class Type>
        void reallocate(Type*& ptr, const int alloc_size);

        void reallocate_input_data();
        void reallocate_output_data();
        
        /*
        if is_copy_input_data_ is true, then allocate the input data memory;
        if false, then simply direct the this->data_in_ ptr to the memory.
        */
        void init_input_data(const InType* data_in);
        /*
        copy the input data if this->is_copy_input_data_ == true.
        */
        void import_input_data(const InType* data_in);

        InType* data_in_ = nullptr;
        OutType* data_out_ = nullptr;

    private:
        /*size of the input data (locally)*/
        int input_size_[DIM]; 
        /*length of the (local) input memory*/
        int input_len_;
        /*size of the output data (locally)*/
        int output_size_[DIM];
        /*length of the (local) output memory */
        int output_len_;

        //status variable
        bool is_copy_input_data_ = false;
        unsigned fft_flags_ = FFTW_ESTIMATE;

    };

    template<int DIM>
    class FFT : public FFTBase<DIM, complex_t, complex_t>{
    public:
        FFT(const int* size,
            const bool copy_input_data = true,
            const unsigned fft_flags = FFTW_ESTIMATE);
        ~FFT();

        void fft(const complex_t* data_in) override;
        void ifft(const complex_t* data_in) override;
    private:
        void fft_impl(const complex_t* data_in, const int sign);
    };

    /*
    Input field is a real field.
    Then the output has the degeneracy: out[i] == out[n-i].
    The length of the output is n/2+1.
    */
    template<int DIM>
    class FFTr2c : public FFTBase<DIM, real_t, complex_t>{
    public:
        FFTr2c(const int* input_size,
            const bool copy_input_data = true,
            const unsigned fft_flags = FFTW_ESTIMATE);
        ~FFTr2c();

        void fft(const real_t* data_in) override;
        void ifft(const real_t* data_in) override {;} //undefined
    private:
        void fft_impl(const real_t* data_in, const int sign);
    };

    template<int DIM>
    class FFTc2r : public FFTBase<DIM, complex_t, real_t>{
    public:
        FFTc2r(const int* output_size,
            const bool copy_input_data = true,
            const unsigned fft_flags = FFTW_ESTIMATE);
        ~FFTc2r();

        void fft(const complex_t* data_in) override {;} //undefined
        void ifft(const complex_t* data_in) override;
    private:
        void fft_impl(const complex_t* data_in, const int sign);
    };

    class FFT1D : public FFTBase<1, complex_t, complex_t>{
    public:
        FFT1D(
            const int size,
            const bool copy_input_data = true,
            const unsigned fft_flags = FFTW_ESTIMATE);
        ~FFT1D();

        void fft(const complex_t* data_in) override;
        void ifft(const complex_t* data_in) override;
    
    private:
        void fft_impl(const complex_t* data_in, const int sign);

    };

    class FFT1Dr2c : public FFTBase<1, real_t, complex_t>{
    public:
        FFT1Dr2c(
            const int input_size,
            const bool copy_input_data = true,
            const unsigned fft_flags = FFTW_ESTIMATE);
        ~FFT1Dr2c();

        void fft(const real_t* data_in) override;
        void ifft(const real_t* data_in) override {;} //undefined
    private:
        void fft_impl(const real_t* data_in, const int sign);
    };

    class FFT1Dc2r : public FFTBase<1, complex_t, real_t>{
    public:
        FFT1Dc2r(
            const int output_size,
            const bool copy_input_data = true,
            const unsigned fft_flags = FFTW_ESTIMATE);
        ~FFT1Dc2r();

        void fft(const complex_t* data_in) override {;} //undefined
        void ifft(const complex_t* data_in) override;
    private:
        void fft_impl(const complex_t* data_in, const int sign);
    };


    template<int DIM, class InType, class OutType = InType>
    class FFTBaseMPI : public FFTBase<DIM, InType, OutType> {
    public:
        FFTBaseMPI();
        ~FFTBaseMPI();
    protected:
    };



    /***
     * Implementation of the FFTBase class
    ***/

    template<int DIM, class InType, class OutType>
    FFTBase<DIM, InType, OutType>::FFTBase(
        const int* input_size,
        const bool copy_input_data, 
        const unsigned fft_flags)
        :
        FFTBase(input_size, input_size, copy_input_data, fft_flags)
    {
    }

    template<int DIM, class InType, class OutType>
    FFTBase<DIM, InType, OutType>::FFTBase(
        const int* input_size,
        const int* output_size,
        const bool copy_input_data, 
        const unsigned fft_flags)
        :
        is_copy_input_data_(copy_input_data),
        fft_flags_(fft_flags) {
        this->resize(input_size, output_size);
    }


    template<int DIM, class InType, class OutType>
    FFTBase<DIM, InType, OutType>::~FFTBase() {
        if(this->is_copy_input_data_){
            this->deallocate<InType>(this->data_in_);
        }
        this->deallocate<OutType>(this->data_out_);
    }


    template<int DIM, class InType, class OutType>
    void FFTBase<DIM, InType, OutType>::resize(
        const int* input_size, const int* output_size){
        if(input_size != nullptr){
            std::copy_n(input_size, DIM, this->input_size_);
            this->input_len_ = std::accumulate(input_size, input_size + DIM, 1, std::multiplies<int>());
            assert(this->input_len_ > 0);
        }
        if(output_size != nullptr){
            std::copy_n(output_size, DIM, this->output_size_);
            this->output_len_ = std::accumulate(output_size, 
                                                output_size + DIM, 
                                                1, 
                                                std::multiplies<int>());
            assert(this->output_len_ > 0);
        }
        return;
    }


    template<int DIM, class InType, class OutType>
    template<class Type>
    void FFTBase<DIM, InType, OutType>::allocate(Type*& ptr, const int alloc_size) {
        assert( ptr == nullptr ); //ptr should be a nullptr to be allocated.
        ptr = static_cast<Type*>( fftw_malloc(sizeof(Type) * alloc_size) );
        return;
    }

    template<int DIM, class InType, class OutType>
    template<class Type>
    void FFTBase<DIM, InType, OutType>::deallocate(Type*& ptr) {
        if( ptr != nullptr ){
            fftw_free(ptr);
            ptr = nullptr;
        }
        return;
    }

    template<int DIM, class InType, class OutType>
    template<class Type>
    void FFTBase<DIM, InType, OutType>::reallocate(Type*& ptr, const int alloc_size) {
        if( ptr == nullptr ) this->allocate<Type>(ptr, alloc_size);
        else {
            this->deallocate<Type>(ptr);
            this->allocate<Type>(ptr, alloc_size);
        }
        return;
    }

    template<int DIM, class InType, class OutType>
    void FFTBase<DIM, InType, OutType>::reallocate_input_data() {
        this->reallocate<InType>(this->data_in_, this->input_len_);
        return;
    }

    template<int DIM, class InType, class OutType>
    void FFTBase<DIM, InType, OutType>::reallocate_output_data() {
        this->reallocate<OutType>(this->data_out_, this->output_len_);
        return;
    }

    template<int DIM, class InType, class OutType>
    void FFTBase<DIM, InType, OutType>::init_input_data(const InType* data_in){
        if(this->is_copy_input_data_){
            this->reallocate_input_data();
        } else{
            this->data_in_ = const_cast<InType*> (data_in);
            // cast away the const-ness of input_data ptr.
            // but check the flag is FFTW_ESTIMATE, which does not change the input data.
            assert(this->fft_flags_ == FFTW_ESTIMATE);
        }
        return;
    }

    template<int DIM, class InType, class OutType>
    void FFTBase<DIM, InType, OutType>::import_input_data(const InType* data_in){
        if(this->is_copy_input_data_){
            std::copy_n(data_in, this->input_len_, this->data_in_);
        }
        return;
    }

    /***
     * Implementations of the FFT derived classes
    ***/

   template<int DIM>
   FFT<DIM>::FFT(const int* size,
            const bool copy_input_data,
            const unsigned fft_flags) 
        :
        FFTBase<DIM, complex_t, complex_t>(size, size, copy_input_data, fft_flags)
    {}

    template<int DIM>
    FFT<DIM>::~FFT() {}

    template<int DIM>
    void FFT<DIM>::fft(const complex_t* data_in){
        this->fft_impl(data_in, FFTW_FORWARD);
    }

    template<int DIM>
    void FFT<DIM>::ifft(const complex_t* data_in){
        this->fft_impl(data_in, FFTW_BACKWARD);
        auto size = this->get_output_size();
        auto norm = static_cast<real_t>(std::accumulate(size, size + DIM, 1, std::multiplies<int>()));
        std::transform(this->data_out_, this->data_out_ + this->get_output_len(), 
                        this->data_out_,
                        [norm](const complex_t& x) { return x / norm; });
        return;
    }

    template<int DIM>
    void FFT<DIM>::fft_impl(const complex_t* data_in, const int sign){
        this->init_input_data(data_in);
        this->reallocate_output_data();
        auto plan = fftw_plan_dft(DIM, 
                                this->get_output_size(), 
                                reinterpret_cast<fftw_complex*>(this->data_in_), 
                                reinterpret_cast<fftw_complex*>(this->data_out_),
                                sign,
                                this->get_fft_flags());
        this->import_input_data(data_in);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        return;
    }

    template<int DIM>
    FFTr2c<DIM>::FFTr2c(
            const int* input_size,
            const bool copy_input_data,
            const unsigned fft_flags)
        :
        FFTBase<DIM, real_t, complex_t>(input_size, input_size, copy_input_data, fft_flags)
    {
        int output_size[DIM];
        std::copy_n(input_size, DIM, output_size);
        output_size[DIM-1] = input_size[DIM-1] / 2 + 1;
        this->resize(nullptr, output_size);
    }

    template<int DIM>
    FFTr2c<DIM>::~FFTr2c() {}

    template<int DIM>
    void FFTr2c<DIM>::fft(const real_t* data_in){
        this->fft_impl(data_in, FFTW_FORWARD);
        return;
    }

    template<int DIM>
    void FFTr2c<DIM>::fft_impl(const real_t* data_in, const int sign){
        this->init_input_data(data_in);
        this->reallocate_output_data();
        auto plan = fftw_plan_dft_r2c(DIM,
                                    this->get_input_size(),
                                    this->data_in_, 
                                    reinterpret_cast<fftw_complex*>(this->data_out_),
                                    this->get_fft_flags());
        this->import_input_data(data_in);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        return;
    }

    template<int DIM>
    FFTc2r<DIM>::FFTc2r(
            const int* output_size,
            const bool copy_input_data,
            const unsigned fft_flags)
        :
        FFTBase<DIM, complex_t, real_t>(output_size, output_size, copy_input_data, fft_flags)
    {
        int input_size[DIM];
        std::copy_n(output_size, DIM, input_size);
        input_size[DIM-1] = output_size[DIM-1] / 2 + 1;
        this->resize(input_size, nullptr);
    }

    template<int DIM>
    FFTc2r<DIM>::~FFTc2r() {}

    template<int DIM>
    void FFTc2r<DIM>::ifft(const complex_t* data_in) {
        this->fft_impl(data_in, FFTW_BACKWARD);
        auto size = this->get_output_size();
        auto norm = static_cast<real_t>(std::accumulate(size, size + DIM, 1, std::multiplies<int>()));
        std::transform(this->data_out_, this->data_out_ + this->get_output_len(), 
                        this->data_out_,
                        [norm](const real_t x) { return x / norm; });
        return;

    }

    template<int DIM>
    void FFTc2r<DIM>::fft_impl(const complex_t* data_in, const int sign){
        this->init_input_data(data_in);
        this->reallocate_output_data();
        auto plan = fftw_plan_dft_c2r(DIM,
                                    this->get_output_size(),
                                    reinterpret_cast<fftw_complex*>(this->data_in_), 
                                    this->data_out_,
                                    this->get_fft_flags());
        this->import_input_data(data_in);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        return;
    }





}


#endif