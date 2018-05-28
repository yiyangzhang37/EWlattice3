import numpy as np
import matplotlib.pyplot as plt
import h5py

class H5Reader:
    '''
    This class reads h5 file with only one dataset.
    '''
    file_ = None
    dataset_ = None
    shape_ = None
    attrs_ = {}
    idx_attrs_ = {}
    def __init__(self, file_path):
        self.file_ = h5py.File(file_path, 'r')
        dataset_name = list(self.file_.keys())[0]
        self.dataset_ = self.file_[dataset_name]
        self.shape_ = self.dataset_.shape

        for (idx, item) in enumerate(self.dataset_.attrs):
            self.attrs_[item] = self.dataset_.attrs[item][0]
            self.idx_attrs_[idx] = self.dataset_.attrs[item][0]
    
    def close(self):
        self.file_.close()

class constants:
    eta_ = 6.0
    dx_ = 0.25
    nSize_ = 128
    lamb_ = 0.129
    
    def __init__(self, eta, dx, nSize):
        self.eta_ = eta
        self.dx_ = dx
        self.nSize_ = nSize
    
    def eta(self):
        return self.eta_
    
    def dx(self):
        return self.dx_
    
    def dt(self):
        return self.dx_*0.25
    
    def nSize(self):
        return self.nSize_
    
    def lamb(self):
        return self.lamb_
    
    def L(self):
        return self.nSize_ * self.dx_
    
    def LatV(self):
        return self.L()**3
    
    def mH(self):
        return 2 * self.lamb_**0.5 * self.eta_

class FourierTransform:
    '''
    xdata_: the x-space data of a vector field. shape = [size_x, size_y. size_z, 3]
    kdata_: the FFT-transformed k-space data of a vector field. shape = [size_x, size_y, size_z, 3]
    shape_: [size_x, size_y, size_z, 3]
    freq_: tuple. the FFT frequency for each dimension. (array of size_x, array of size_y, array of size_z)
    kenergy_: the enengy density in k-space. shape = [size_x, size_y. size_z]
    '''
    xdata_ = None
    kdata_ = None
    shape_ = None
    freq_ = None
    kenergy_ = None
    
    def __init__(self, dataset):
        self.xdata_ = np.array(dataset)
        self.shape_ = self.xdata_.shape
        assert(np.array(self.shape_).shape[0] == 4), "Dataset dimension should be 4."
        assert(self.shape_[-1] == 3), "Dataset should have 3 components."
    
    def plot_energy_spectrum(self, sep_number, 
                             norm_factor = None, 
                             ax = None, 
                             is_zero_mode_separated = True):
        self.fft()
        self.energy_spectrum()
        seps, vals = self.configure_magnitude_bins(sep_number)
        spec = self.radial_spectrum(sep_number)
        if norm_factor is not None:
            spec = spec * norm_factor
        if ax is not None:
            ax.bar(vals, spec, align='center', width = 0.01)
        else:
            plt.bar(vals, spec, align='center', width = 0.01)
        return seps, spec
        
    def fft(self):
        #kdata_0 = np.fft.fftshift( np.fft.fftn(self.xdata_[:,:,:,0]) )
        #kdata_1 = np.fft.fftshift( np.fft.fftn(self.xdata_[:,:,:,1]) )
        #kdata_2 = np.fft.fftshift( np.fft.fftn(self.xdata_[:,:,:,2]) )
        kdata_0 = np.fft.fftn(self.xdata_[:,:,:,0])
        kdata_1 = np.fft.fftn(self.xdata_[:,:,:,1])
        kdata_2 = np.fft.fftn(self.xdata_[:,:,:,2])
        self.kdata_ = np.stack([kdata_0, kdata_1, kdata_2], axis = 3)
        #freq_0 = np.fft.fftshift( np.fft.fftfreq(self.shape_[0]) )
        #freq_1 = np.fft.fftshift( np.fft.fftfreq(self.shape_[1]) )
        #freq_2 = np.fft.fftshift( np.fft.fftfreq(self.shape_[2]) )
        freq_0 = np.fft.fftfreq(self.shape_[0])
        freq_1 = np.fft.fftfreq(self.shape_[1])
        freq_2 = np.fft.fftfreq(self.shape_[2])
        self.freq_ = (freq_0, freq_1, freq_2)
    
    def ifft(self):
        xdata_0 = np.fft.ifftn(self.kdata_[:,:,:,0])
        xdata_1 = np.fft.ifftn(self.kdata_[:,:,:,1])
        xdata_2 = np.fft.ifftn(self.kdata_[:,:,:,2])
        return np.stack([xdata_0, xdata_1, xdata_2], axis = 3)
    
    def energy_spectrum(self):
        '''Compute the energy spectrum.
        kenergy_ = |B_k1|^2 + |B_k2|^2 + |B_k3|^2
        '''
        #self.kenergy_ = np.square(np.absolute(self.kdata_[:,:,:,0])) \
        #            + np.square(np.absolute(self.kdata_[:,:,:,1])) \
        #            + np.square(np.absolute(self.kdata_[:,:,:,2]))
        self.kenergy_ = np.sum(np.square(np.absolute(self.kdata_)), axis = 3)
    
    def configure_magnitude_bins(self, sep_number, is_zero_mode_separated = True):
        '''
        sep_number: the number of separation points.
            if is_zero_mode_separated == True, then sep_number == bin_number;
            otherwise, sep_number == bin_number + 1.
        generate the separation points of bins, and the mean value of the bins.
        return: (sep_points, mean_values)
        sep_points is not affected by the parameter is_zero_mode_separated.
        '''
        seps = np.linspace(0, self.max_k(), sep_number)
        if is_zero_mode_separated:
            vals = np.zeros(sep_number)
            vals[0] = 0
            for i in range(sep_number-1):
                vals[i+1] = 0.5*(seps[i] + seps[i+1])
        else:
            vals = np.zeros(sep_number - 1)
            for i in range(sep_number-1):
                vals[i] = 0.5*(seps[i] + seps[i+1])
        return seps, vals
    
    def radial_spectrum(self, sep_number, is_zero_mode_separated = True):
        '''
        compute the radial spectrum.
        '''
        seps, vals = self.configure_magnitude_bins(sep_number, is_zero_mode_separated)
        # The radial spectrum is calculated first supposing is_zero_mode_separated == True.
        radial_spectrum = np.zeros(sep_number, dtype = float)
        k_mag = self.k_magnitude()
        flatk = k_mag.flatten()
        flate = self.kenergy_.flatten()
        for k, e in zip(flatk, flate):
            bn = self.query_bin_number(k, seps)
            radial_spectrum[bn] += e
        #for xi in range(self.shape_[0]):
        #    for yi in range(self.shape_[1]):
        #        for zi in range(self.shape_[2]):
        #            k_mag = np.sqrt(self.freq_[0][xi]**2 \
        #                              + self.freq_[1][yi]**2 \
        #                              + self.freq_[2][zi]**2)
        #            bn = self.query_bin_number(k_mag, seps)
        #            radial_spectrum[bn] += self.kenergy_[xi, yi, zi]
        if(is_zero_mode_separated == False):
            radial_spectrum[1] = radial_spectrum[0] + radial_spectrum[1]
            radial_spectrum[0] = 0
        return radial_spectrum
        
    def query_bin_number(self, val, seps):
        '''
        bin_seps: ndarray of bin separation points.
        Return: i, if bin_seps[i-1] < val <= bin_seps[i].
        it returns 0 only if val == 0.
        '''
        _seps = np.array(seps)
        if(val == _seps[0]):
            return 0
        else:
            lp = 0
            rp = _seps.shape[0] - 1
            mp = (lp+rp)//2
            while rp - lp > 1:
                if val<= _seps[mp]:
                    rp = mp
                else:
                    lp = mp
                mp = (lp+rp)//2
            return rp
    
    
    def k_mean(self):
        '''
        calculate the mean value of k, weighted by self.kenergy_.
        '''
        gridx, gridy, gridz = self.build_kspace_mesh()
        kspec = np.sum( self.kenergy_ * np.sqrt( np.square(gridx) + np.square(gridy) + np.square(gridz) ) )
        return kspec / self.kenergy_.sum()
        
    def build_kspace_mesh(self):
        '''
        build a k-space from self.freq_.
        return: tuple, k-space mesh, (gridx, gridy, gridz), each with shape [size_x, size_y, size_z].
        '''
        return np.meshgrid(self.freq_[0], self.freq_[1], self.freq_[2], indexing = 'ij')
   
    def k_magnitude(self):
        gridx, gridy, gridz = self.build_kspace_mesh()
        return np.sqrt( np.square(gridx) + np.square(gridy) + np.square(gridz) )
        
    def max_k(self):
        '''
        return the maximum magnitude of k in the current FFT.
        '''
        max0 = np.absolute(self.freq_[0]).max()
        max1 = np.absolute(self.freq_[1]).max()
        max2 = np.absolute(self.freq_[2]).max()
        return (max0**2 + max1**2 + max2**2)**0.5