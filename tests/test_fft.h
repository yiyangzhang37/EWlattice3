#include "../ParaSite/FFTWrapper.h"
#include <array>
#include <complex>

using namespace FFT_Wrapper;

// I am just lazy. Don't want to write another cpp file.
template<int _ = 0>
void test_FFT1D(){

    const int N = 128;
    FFT1D fourier(N, true, FFTW_MEASURE);
    std::array<complex_t, N> in;
    for(auto i = 0; i < N; ++i){
        in[i] = std::complex<double>(sin(2.0*3.1415926/N*i), 0);
    }
    fourier.fft(in.data());
    auto out_ptr = fourier.get_output_data();
    for(auto i = 0; i < N; ++i){
        std::cout << out_ptr[i] << std::endl;
    }
    return;
}

template<int _ = 0>
void test_FFT1D_fft_ifft(){
    const int N = 128;
    FFT1D fourier(N, true, FFTW_MEASURE);
    std::array<complex_t, N> in;
    for(auto i = 0; i < N; ++i){
        in[i] = std::complex<double>(sin(2.0*3.1415926/N*i), 0);
    }
    fourier.fft(in.data());
    std::array<complex_t, N> out;
    auto out_ptr = fourier.get_output_data();
    std::copy_n(out_ptr, N, out.begin());
    fourier.ifft(out.data());
    out_ptr = fourier.get_output_data();
    for(auto i = 0; i < N; ++i){
        std::cout << in[i] << ":::" << out_ptr[i] << std::endl;
    }
    return;
}

template<int _ = 0>
void test_FFT1D_r2c_c2r() {
    const int N = 128;
    FFT1Dr2c fr2c(N, true, FFTW_MEASURE);
    FFT1Dc2r fc2r(N, true, FFTW_MEASURE);
    std::array<real_t, N> in;
    std::array<complex_t, N> out;
    for(auto i = 0; i < N; ++i){
        in[i] = cos(2.0*3.1415926/N*i);
    }
    fr2c.fft(in.data());
    auto out_ptr = fr2c.get_output_data();
    std::copy_n(out_ptr, N, out.begin());
    //for(auto i = 0; i < N; ++i){
    //    std::cout << out_ptr[i] << std::endl;
    //}
    fc2r.ifft(out.data());
    auto out_ptr2 = fc2r.get_output_data();
    for(auto i = 0; i < N; ++i){
        std::cout << in[i] << ":::" << out_ptr2[i] << std::endl;
    }
    return;
}

template<int _ = 0>
void test_FFT3D() {
    double PI = 3.1415927;
    int SIZE[3] = {16, 16, 16};
    int LEN = SIZE[0]*SIZE[1]*SIZE[2];
    int HALF_LEN = SIZE[0]*SIZE[1]*(SIZE[2]/2+1);
    std::vector<complex_t> complex_in(LEN);
    std::vector<complex_t> complex_out(LEN);
    std::vector<real_t> real_in(LEN);

    for(auto i = 0; i < LEN; ++i) {
        auto x = i / (SIZE[1]*SIZE[2]);
        auto y = (i % (SIZE[1]*SIZE[2])) / SIZE[2];
        auto z = (i % (SIZE[1]*SIZE[2])) % SIZE[2];
        real_in[i] = sin(4.0*PI*x/SIZE[0]) * cos(2.0*PI*y/SIZE[1]) * sin(10.0*PI*z/SIZE[2]);
        complex_in[i] = complex_t(real_in[i], 0);
    }
    FFT<3> fourier(SIZE, true, FFTW_MEASURE);
    fourier.fft(complex_in.data());
    auto out_ptr = fourier.get_output_data();
    std::copy(out_ptr, out_ptr + LEN, complex_out.begin());
    fourier.ifft(complex_out.data());
    auto out_ptr_back = fourier.get_output_data();

    for(auto i = 0; i < LEN; ++i) {
        std::cout << complex_in[i] << "::::" << out_ptr_back[i] << std::endl;
    }
    return;
}


template<int _ = 0>
void test_FFT3D_r2c() {
    double PI = 3.1415927;
    int SIZE[3] = {16, 16, 16};
    int LEN = SIZE[0]*SIZE[1]*SIZE[2];
    int HALF_LEN = SIZE[0]*SIZE[1]*(SIZE[2]/2+1);
    std::vector<complex_t> complex_in(LEN);
    std::vector<complex_t> complex_out(LEN);
    std::vector<real_t> real_in(LEN);
    std::vector<complex_t> r2c_out(HALF_LEN);

    for(auto i = 0; i < LEN; ++i) {
        auto x = i / (SIZE[1]*SIZE[2]);
        auto y = (i % (SIZE[1]*SIZE[2])) / SIZE[2];
        auto z = (i % (SIZE[1]*SIZE[2])) % SIZE[2];
        real_in[i] = sin(4.0*PI*x/SIZE[0]) * cos(2.0*PI*y/SIZE[1]) * sin(10.0*PI*z/SIZE[2]);
        complex_in[i] = complex_t(real_in[i], 0);
    }

    FFT<3> fourier(SIZE, true, FFTW_MEASURE);
    fourier.fft(complex_in.data());
    auto out_ptr = fourier.get_output_data();
    std::copy(out_ptr, out_ptr + LEN, complex_out.begin());

    FFTr2c<3> fr2c(SIZE, true, FFTW_MEASURE);
    fr2c.fft(real_in.data());
    auto out_ptr2 = fr2c.get_output_data();
    std::copy(out_ptr2, out_ptr2 + HALF_LEN, r2c_out.begin());

    int r2c_i = 0;
    for(auto i = 0; i < LEN; ++i) {
        auto x = i / (SIZE[1]*SIZE[2]);
        auto y = (i % (SIZE[1]*SIZE[2])) / SIZE[2];
        auto z = (i % (SIZE[1]*SIZE[2])) % SIZE[2];
        if(z < SIZE[2]/2 + 1){
            std::cout << complex_out[i] << ":::" << r2c_out[r2c_i++] << std::endl;
        }
    }

    return;
}

template<int _ = 0>
void test_FFT3D_r2cc2r() {
    double PI = 3.1415927;
    int SIZE[3] = {16, 16, 16};
    int LEN = SIZE[0]*SIZE[1]*SIZE[2];
    int HALF_LEN = SIZE[0]*SIZE[1]*(SIZE[2]/2+1);
    std::vector<real_t> real_in(LEN);
    std::vector<complex_t> r2c_out(HALF_LEN);

    for(auto i = 0; i < LEN; ++i) {
        auto x = i / (SIZE[1]*SIZE[2]);
        auto y = (i % (SIZE[1]*SIZE[2])) / SIZE[2];
        auto z = (i % (SIZE[1]*SIZE[2])) % SIZE[2];
        real_in[i] = sin(4.0*PI*x/SIZE[0]) * cos(2.0*PI*y/SIZE[1]) * sin(10.0*PI*z/SIZE[2]);
    }

    FFTr2c<3> fr2c(SIZE, true, FFTW_MEASURE);
    fr2c.fft(real_in.data());
    auto out_ptr2 = fr2c.get_output_data();
    std::copy(out_ptr2, out_ptr2 + HALF_LEN, r2c_out.begin());

    FFTc2r<3> fc2r(SIZE, true, FFTW_ESTIMATE);
    fc2r.ifft(r2c_out.data());
    auto out_ptr3 = fc2r.get_output_data();
    
    for(auto i = 0; i < LEN; ++i) {
        std::cout << real_in[i] << "::::" << out_ptr3[i] << std::endl;
    }

    return;
}



