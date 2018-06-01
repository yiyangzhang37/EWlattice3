#include "../ParaSite/FFTWrapper.h"
#include <array>

using namespace FFT_Wrapper;

// I am just lazy. Don't want to write another cpp file.
template<int _ = 0>
void test_ComplexFFT1D(){

    const int N = 128;
    ComplexFFT1D fourier(N);
    std::array<double, N> in;
    std::array<double, N> out;
    return;

}