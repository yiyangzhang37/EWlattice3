#ifdef FFT3D
#include "LATfield2.hpp"
using namespace LATfield2;

temporaryMemFFT::temporaryMemFFT()
{
	allocated_ = 0;
}
temporaryMemFFT::~temporaryMemFFT()
{
	if (allocated_ != 0) {
#ifdef SINGLE
		fftwf_free(temp1_);
		fftwf_free(temp2_);
#endif
#ifndef SINGLE
		fftw_free(temp1_);
		fftw_free(temp2_);
#endif
	}
}

temporaryMemFFT::temporaryMemFFT(long size)
{
#ifdef SINGLE
	temp1_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
	temp2_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
#endif

#ifndef SINGLE
	temp1_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
	temp2_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
#endif
	allocated_ = size;
}
int temporaryMemFFT::setTemp(long size)
{
	if (size>allocated_)
	{
#ifdef SINGLE
		if (allocated_ != 0) {
			fftwf_free(temp1_);
			fftwf_free(temp2_);
		}
		temp1_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));
		temp2_ = (fftwf_complex *)fftwf_malloc(size*sizeof(fftwf_complex));

		//debug cout<<"("<< parallel.grid_rank()[0]<< ","<<parallel.grid_rank()[1] <<"): called temporary resize. Old size: " <<allocated_<<" , new size: "<< size<<endl;

#endif

#ifndef SINGLE
		if (allocated_ != 0) {
			fftw_free(temp1_);
			fftw_free(temp2_);
		}
		temp1_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
		temp2_ = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
#endif
		allocated_ = size;


	}
	return 1;
}


//////////////////////Temp memory///////////////////////////

temporaryMemFFT tempMemory;

#endif