#ifndef PARASITE_FFT_H
#define PARASITE_FFT_H

#include <algorithm>
#include <type_traits>
#include "ParaSite_Field.h"
#include "FFTWrapper.h"
#include "HDF5Wrapper.h"

namespace ParaSite{

    using ParallelObject = MPI_Wrapper::Parallel2D;
    using FileProcessor = HDF5_Wrapper::HDF5Wrapper;
    using namespace FFT_Wrapper;

    /* Custom type traits: complex type */
    template<class T>
    struct is_complex_type : std::false_type {};
    
    template<>
    struct is_complex_type<complex_t> : std::true_type {};
    
    /*
    Perform FFT for a component of the ParaSite:Field object.
    The FFT is done on the root process, by first gathering data from all the processes.
    */
    template<int DIM, class TypeIn, class TypeOut = TypeIn>
    class FieldFFT {
    private:
        /*
        The global size of the field.
        ***!!!*** right now, Field<> class and this class have different conventions.
        So the size_ is in reverse order compared with Field<> class.
        */
        IndexType size_[DIM];
        /*
        The number of global sites of the field.
        */
        IndexType len_;

        std::vector<TypeIn> data_in_;
        std::vector<TypeOut> data_out_;

        const Field<TypeIn, DIM>& field_;
        const int row_;
        const int col_;

        const ParallelObject& parallel_object_;

        void gather_data();

        void fft_floating_type();
        void fft_complex_type();
        void ifft_floating_type();
        void ifft_complex_type();

    public:
        /*
        Take a Field object, and do the FFT on the (row, col)-component.
        */
        FieldFFT(const Field<TypeIn, DIM>& field, 
                const int row, 
                const int col);
        ~FieldFFT() = default;

        // The transforms require template specifications.
        // The following two have specifications:
        // 1. TypeIn == complex_t, TypeOut == complex_t
        // 2. TypeIn == some floating type, TypeOut == complex_t
        void fft();
        void ifft();
        
        // The following has the specifications:
        // 1. TypeIn == some floating type, TypeOut == complex_t
        void fft_r2c();

        // The following has the specifications:
        // 1. TypeIn == complex_t, TypeOut == real_t
        void ifft_c2r();
        
        // save the transformed data to an HDF5 file.
        // This function only performs on root process.
        void save(const std::string& file_name, const std::string& dataset_name);
    };

    template<int DIM, class TypeIn, class TypeOut>
    FieldFFT<DIM, TypeIn, TypeOut>::FieldFFT(
        const Field<TypeIn, DIM>& field, 
        const int row, 
        const int col)
        : 
        field_(field),
        row_(row),
        col_(col),
        parallel_object_(field.get_parallel_object())
    {   
        const auto& lat_size = field.get_lattice().get_global_size();
        std::reverse_copy(lat_size, lat_size + DIM, this->size_);
        this->len_ = field.get_lattice().get_visible_global_sites();
        this->gather_data();
    }

    template<int DIM, class TypeIn, class TypeOut>
    void FieldFFT<DIM, TypeIn, TypeOut>::gather_data(){
        this->data_in_ = this->field_.serialize(this->row_, this->col_);
        return;
    }

    template<int DIM, class TypeIn, class TypeOut>
    void FieldFFT<DIM, TypeIn, TypeOut>::fft_floating_type() {
        if(this->parallel_object_.is_root()){
            FFT<DIM> fourier(this->size_, true, FFTW_MEASURE);
            std::vector<complex_t> complex_in(this->data_in_.size());
            std::transform(this->data_in_.begin(), this->data_in_.end(), 
                            complex_in.begin(), 
                            [](TypeIn x){return complex_t(x, 0);});
            fourier.fft(complex_in.data());
            auto fft_ptr = fourier.get_output_data();
            this->data_out_ = std::vector<complex_t>(fft_ptr, fft_ptr + this->len_);
        }
        this->parallel_object_.Barrier();
        return;
    }

    template<int DIM, class TypeIn, class TypeOut>
    void FieldFFT<DIM, TypeIn, TypeOut>::fft_complex_type() {
        if(this->parallel_object_.is_root()){
            FFT<DIM> fourier(this->size_, true, FFTW_MEASURE);
            fourier.fft(this->data_in_.data());
            auto fft_ptr = fourier.get_output_data();
            this->data_out_ = std::vector<complex_t>(fft_ptr, fft_ptr + this->len_);
        }
        this->parallel_object_.Barrier();
        return;
    }

    template<int DIM, class TypeIn, class TypeOut>
    void FieldFFT<DIM, TypeIn, TypeOut>::ifft_floating_type() {
        if(this->parallel_object_.is_root()){
            FFT<DIM> fourier(this->size_, true, FFTW_MEASURE);
            std::vector<complex_t> complex_in(this->data_in_.size());
            std::transform(this->data_in_.begin(), this->data_in_.end(), 
                            complex_in.begin(), 
                            [](TypeIn x){return complex_t(x, 0);});
            fourier.ifft(complex_in.data());
            auto fft_ptr = fourier.get_output_data();
            this->data_out_ = std::vector<complex_t>(fft_ptr, fft_ptr + this->len_);
        }
        this->parallel_object_.Barrier();
        return;
    }
    
    template<int DIM, class TypeIn, class TypeOut>
    void FieldFFT<DIM, TypeIn, TypeOut>::ifft_complex_type() {
        if(this->parallel_object_.is_root()){
            FFT<DIM> fourier(this->size_, true, FFTW_MEASURE);
            fourier.ifft(this->data_in_.data());
            auto fft_ptr = fourier.get_output_data();
            this->data_out_ = std::vector<complex_t>(fft_ptr, fft_ptr + this->len_);
        }
        this->parallel_object_.Barrier();
        return;
    }


    template<int DIM, class TypeIn, class TypeOut>
    void FieldFFT<DIM, TypeIn, TypeOut>::fft() {
        if(std::is_floating_point<TypeIn>::value == true) {
            this->fft_floating_type();
            return;
        }
        if(is_complex_type<TypeIn>::value == true) {
            this->fft_complex_type();
            return;
        }
        std::cout << "FieldFFT::fft requires a floating type or complex_t type." << std::endl;
        return;
    }   

    template<int DIM, class TypeIn, class TypeOut>
    void FieldFFT<DIM, TypeIn, TypeOut>::ifft() {
        if(std::is_floating_point<TypeIn>::value == true) {
            this->ifft_floating_type();
            return;
        }
        if(is_complex_type<TypeIn>::value == true) {
            this->ifft_complex_type();
            return;
        }
        std::cout << "FieldFFT::ifft requires a floating type or complex_t type." << std::endl;
        return;
        
    }

    template<int DIM, class TypeIn, class TypeOut>
    void FieldFFT<DIM, TypeIn, TypeOut>::fft_r2c() {
        if(std::is_floating_point<TypeIn>::value == false) {
            std::cout << "FieldFFT::fft_r2c requires a floating type." << std::endl;
            return;
        }
        if(this->parallel_object_.is_root()){
            FFTr2c<DIM> fourier(this->size_, true, FFTW_MEASURE);
            std::vector<real_t> real_in(this->len_);
            std::transform(this->data_in_.begin(), this->data_out_.end(),
                            real_in.begin(), 
                            [](TypeIn x){return static_cast<real_t>(x); });
            fourier.fft(real_in);
            auto fft_ptr = fourier.get_output_data();
            this->data_out_ = std::vector<complex_t>(fft_ptr, fft_ptr + fourier.get_output_len());
        }
        this->parallel_object_.Barrier();
    }

    template<int DIM, class TypeIn, class TypeOut>
    void FieldFFT<DIM, TypeIn, TypeOut>::ifft_c2r() {
        std::cout << "Not implemented yet." << std::endl;
        return;
    }

    template<int DIM, class TypeIn, class TypeOut>
    void FieldFFT<DIM, TypeIn, TypeOut>::save(
        const std::string& file_name, 
        const std::string& dataset_name) {
        FileProcessor fp;

        std::string dst_name;
        if(dataset_name == "") dst_name = file_name;
        else dst_name = dataset_name;

        if(this->parallel_object_.is_root()){
            fp.SaveSingleDatasetFile<TypeOut, DIM>(
                                file_name,
                                dst_name,
                                this->size_, /*dataset_size*/
                                this->data_out_.data());
            fp.open_file(file_name);
            fp.open_dataset(dst_name);
            //int count = 0;
            for(auto i = 0; i < DIM; ++i){
                fp.attach_attribute_to_dataset<std::string>("col_" + std::to_string(i), 
                                                    "axis_" + std::to_string(DIM-1-i));
            } 
            fp.close_dataset();
            fp.close_file();
        }
        this->parallel_object_.Barrier();
        return;
    }



    using FieldFFT_float3d = FieldFFT<3, float, complex_t>;
    using FieldFFT_double3d = FieldFFT<3, double, complex_t>;
    using FieldFFT_complex3d = FieldFFT<3, complex_t, complex_t>;

}


#endif