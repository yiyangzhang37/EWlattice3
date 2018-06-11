import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import h5py
import EWdata
from multiprocessing import Pool

root_path = '/data/yiyang.zhang/EWlattice3'
c1 = EWdata.constants(eta = 1.0, dx = 0.25, nSize = 256)
BINS = 120

files = ['nr_1_bfield_' + str(i) +'.h5' for i in range(16, 4000, 50)]

def fft(f):
    file_path = os.path.join(root_path, f)
    h5file = EWdata.H5Reader(file_path)
    dset = h5file.dataset_
    dshape = h5file.shape_

    fourier = EWdata.FourierTransform(dset)
    fourier.fft()
    fourier.energy_spectrum()
    seps, vals = fourier.configure_magnitude_bins(BINS)
    DK = seps[1] - seps[0]
    s = fourier.radial_spectrum(BINS)
    # norm_factor = 1/(2.0*(c1.nSize()/c1.dx())**3) * 1/(c1.LatV()*c1.mH()**4)
    # proper normalization: such that radial.sum() == total energy
    s = s /(2.0*(c1.nSize()/c1.dx())**3)
    # normalize to energy density in units of m_H^4
    s = s / (c1.LatV()*c1.mH()**4)
    ## normalize to area
    s = s / DK
    kmean = fourier.k_mean()
    # add the last item of both vals and s as the k_mean value.
    vals = np.append(vals, kmean)
    s = np.append(s, kmean)
    output = np.array([vals, s])
    output_name = os.path.join(root_path, f[:-3] + '.npy')
    np.save(output_name, output)

def fft_single_process():
    for f in files:
        fft(f)

def fft_multi_process():
    pool = Pool(8)
    pool.map(fft, files)

    
if __name__ == '__main__':
    fft_multi_process()

