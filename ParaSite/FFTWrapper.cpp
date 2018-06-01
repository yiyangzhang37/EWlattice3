#include "FFTWrapper.h"
#include <algorithm>
#include <functional>
#include <iostream>
#include <cassert>

template<int DIM, class InType, class OutType>
FFTBase::FFTBase(const int* size){
    this->resize(size);
}

template<int DIM, class InType, class OutType>
FFTBase::~FFTBase() {
    this->deallocate();
}


template<int DIM, class InType, class OutType>
void FFTBase::resize(const int* size){
    std::copy_n(size, DIM, this->size_);
    this->len_ = std::accumulate(size, size + DIM, 1, std::multiplies());
    assert(this->len_ > 0);
    return;
}


template<int DIM, class InType, class OutType>
void FFTBase::allocate() {
    if( this->data_out_ != nullptr ){
        std::cerr << "Memory already been allocated." << std::endl;
    }
    this->data_out_ = static_cast<OutType*>( fftw_malloc(sizeof(OutType)* this->len_) );
    return;
}

template<int DIM, class InType, class OutType>
void FFTBase::deallocate() {
    if( this->data_out_ == nullptr ){
        std::cerr << "No memory was allocated." << std::endl;
    }
    fftw_free(this->data_out_);
    this->data_out_ = nullptr;
    return;
}

template<int DIM, class InType, class OutType>
void FFTBase::reallocate() {
    this->deallocate();
    this->allocate();
    return;
}


