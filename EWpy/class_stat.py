# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 11:29:45 2017

@author: YiyangZhang
"""

import numpy as np
import scipy.spatial as spatial
import pandas as pd
import codecs
import os.path
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
import re

import time

'''class of FFT'''

class FourierTransform:
    # The dict is used to translate the data_key into proper file name extensions.
    data_dict = {'Bx': '.Bx.datfast', \
                 'By': '.By.datfast', \
                 'Bz': '.Bz.datfast', \
                 'Ex': '.Ex.datfast', \
                 'Ey': '.Ey.datfast', \
                 'Ez': '.Ez.datfast', \
                 'EME': '.EME.datfast', \
                 'EMB': '.EMB.datfast'};
    DIM = 3;
                 
    def __init__(self, file_path, run_id, lat_size, stop_point = ''):
        '''
        file_path: str. file path of the field data files.
        run_id: str. 
        lat_size: array-like, int. size = [DIM(3)]
        stop_point: str. specify when the data is taken. either can specify time or minHiggs2.
        '''
        self.file_path = file_path;
        self.run_id = run_id;
        self.stop_point = stop_point;
        self.nDim = lat_size;
        self.is_ffted = False;

    def fft(self, vec_id):
        '''FFT for one field, or one component of a vector field.
        It does not apply to fft of energy.
        ---
        vec_id: str. used to retrieve the proper file extension from data_dict. 
        '''
        if self.stop_point == 0:
            self.x_data = np.fromfile(self.file_path + self.run_id + \
                                      FourierTransform.data_dict[vec_id], \
                                      dtype=np.dtype('float32'));
        else:
            self.x_data = np.fromfile(self.file_path + self.run_id + \
                                      '.' + str(self.stop_point) + \
                                      FourierTransform.data_dict[vec_id], \
                                      dtype=np.dtype('float32')); 
        #This is the correct way to reshape the data
         #the element access is given by data[xi,yi,zi]
        self.x_data = self.x_data.reshape((self.nDim,self.nDim,self.nDim)).transpose();
        self.k_data = np.fft.fftn(self.x_data)
        self.freq = np.fft.fftfreq(self.nDim)
        #The k space data are shifted to put 0 at center
        self.k_data = np.fft.fftshift(self.k_data)
        self.freq = np.fft.fftshift(self.freq)
        self.is_ffted = True
        
    def fft_ensemble(self, vec_ids, ensemble_size):
        ''' FFT Transform of an ensemle of (3D) vector field
        The function will directly compute the ensemble average of energy spectrum (self.k_energy)
        ---
        vec_ids: array-like, [3]. used to retrieve the proper file extension from data_dict.
        ensemble_size: int. The size of an ensemble.
        '''
        run_id_ = self.run_id
        ensemble_spec = np.zeros((self.nDim,self.nDim,self.nDim))
        for s in range(ensemble_size):
            self.run_id = run_id_ + str(s+1).zfill(3)
            self.fft_vector(vec_ids)
            self.fft_energy_spec()
            ensemble_spec += self.k_energy
        # average over the ensemble
        self.k_energy = ensemble_spec / ensemble_size
        # rewrite the run_id back
        self.run_id = run_id_
  
    def fft_1d(self, axis):
        sum_axis1 = (axis+1) % FourierTransform.DIM
        sum_axis2 = (axis+2) % FourierTransform.DIM
        return np.sum(self.k_data, axis=(sum_axis1, sum_axis2))
    
    def fft_vector(self, vec_ids):
        '''FFT Transform of a vector field (3D)
        vec_ids: array-like, [3]. use to retrieve the proper file extension from data_dict.
        '''
        if self.stop_point == 0:
            self.x_data_1 = np.fromfile(self.file_path + self.run_id + \
                                         FourierTransform.data_dict[vec_ids[0]], \
                                         dtype=np.dtype('float32'))
            self.x_data_2 = np.fromfile(self.file_path + self.run_id + \
                                         FourierTransform.data_dict[vec_ids[1]], \
                                         dtype=np.dtype('float32'))
            self.x_data_3 = np.fromfile(self.file_path + self.run_id + \
                                         FourierTransform.data_dict[vec_ids[2]], \
                                         dtype=np.dtype('float32'))
        else:
            self.x_data_1 = np.fromfile(self.file_path + self.run_id + \
                                        '.' + str(self.stop_point) + \
                                         FourierTransform.data_dict[vec_ids[0]], \
                                         dtype=np.dtype('float32'))
            self.x_data_2 = np.fromfile(self.file_path + self.run_id + \
                                        '.' + str(self.stop_point) + \
                                         FourierTransform.data_dict[vec_ids[1]], \
                                         dtype=np.dtype('float32'))
            self.x_data_3 = np.fromfile(self.file_path + self.run_id + \
                                        '.' + str(self.stop_point) + \
                                         FourierTransform.data_dict[vec_ids[2]], \
                                         dtype=np.dtype('float32'))
        self.x_data_1 = self.x_data_1.reshape((self.nDim,self.nDim,self.nDim)).transpose()
        self.x_data_2 = self.x_data_2.reshape((self.nDim,self.nDim,self.nDim)).transpose()
        self.x_data_3 = self.x_data_3.reshape((self.nDim,self.nDim,self.nDim)).transpose()
        self.k_data_1 = np.fft.fftn(self.x_data_1)
        self.k_data_2 = np.fft.fftn(self.x_data_2)
        self.k_data_3 = np.fft.fftn(self.x_data_3)
        self.freq = np.fft.fftfreq(self.nDim)
        #The k space are shifted to put 0 at center
        self.k_data_1 = np.fft.fftshift(self.k_data_1)
        self.k_data_2 = np.fft.fftshift(self.k_data_2)
        self.k_data_3 = np.fft.fftshift(self.k_data_3)
        self.freq = np.fft.fftshift(self.freq)
        self.is_ffted = True
        
    def fft_energy_spec(self):
        '''Compute the energy spectrum.
        k_energy = |E_k1|^2 + |E_k2|^2 + |E_k3|^2
        '''
        self.k_energy = np.square(np.absolute(self.k_data_1)) \
            + np.square(np.absolute(self.k_data_2)) \
            + np.square(np.absolute(self.k_data_3))
    
    def fft_radial_kdtree(self, xyz_data, bin_number):
        bins = self.bin_seps(bin_number);
        bin_size = bin_number;
        radial_spectrum = np.zeros(bin_size, dtype = complex);
        '''first consider zero modes separately'''
        coord = [ np.where(self.freq==0)[0][0],  \
                         np.where(self.freq==0)[0][0], np.where(self.freq==0)[0][0] ];
        radial_spectrum[0] = xyz_data[coord[0], coord[1], coord[2]];
        '''then consider other modes.'''
        for i in np.arange(bin_size-1):
            shell_points_index = self.query_shell_points([0,0,0], bins[i+1], bins[i]);
            shell_points = self.mesh[shell_points_index];
            sum_coeff = 0;
            for p in shell_points:
                coord = [ np.where(self.freq==p[0])[0][0],  \
                         np.where(self.freq==p[1])[0][0], \
                         np.where(self.freq==p[2])[0][0] ];
                sum_coeff += xyz_data[coord[0], coord[1], coord[2]];
            radial_spectrum[i+1] = sum_coeff;
        return np.absolute(radial_spectrum);
    
    def fft_radial(self, xyz_data, bin_number, zero_mode_sep = True):
        '''
        xyz_data: 
        bin_number: int. Number of bins.
        zero_mode_sep: bool. Whether separate the first bin as to count zero mode only.
        if zero_mode_sep is False, then the first bin is always zero.
        ---
        Return: ndarray, float. radial spectrum.
        '''
        bin_sep = self.bin_seps(bin_number);
        radial_spectrum = np.zeros(bin_number, dtype = float);
        for xi in range(self.nDim):
            for yi in range(self.nDim):
                for zi in range(self.nDim):
                    k_mag = math.sqrt(self.freq[xi]**2 \
                                      + self.freq[yi]**2 \
                                      + self.freq[zi]**2)
                    bn = self.query_bin_number(k_mag, bin_sep)
                    radial_spectrum[bn] += xyz_data[xi, yi, zi]
        if(zero_mode_sep == False):
            radial_spectrum[1] = radial_spectrum[0] + radial_spectrum[1]
            radial_spectrum[0] = 0
        return radial_spectrum
                    
        
    def query_bin_number(self, val, bin_seps):
        '''
        bin_seps: ndarray of bin separation points.
        Return: i, if bin_seps[i-1] < val <= bin_seps[i].
        it returns 0 only if val == 0.
        '''
        if(val == bin_seps[0]):
            return 0
        else:
            lp = 0
            rp = len(bin_seps) - 1
            mp = (lp+rp)//2
            while rp - lp > 1:
                if val<= bin_seps[mp]:
                    rp = mp
                else:
                    lp = mp
                mp = (lp+rp)//2
            return rp
    
    def bin_seps(self, bin_number):
        '''
        Return: ndarray, float; len = bin_number. The separation points of the bins
        '''
        return np.linspace(0, self.max_freq_magnitude(), bin_number);
    
    def fft_bins(self, bin_number):
        ''' Get the bin labels (averave value) for each bin.
        The first bin is used exclusively for zero modes.
        So that the len(bin_label) ==  len(bin_seps)
        the bin_sep are the starting and end points for each bin.
        ---
        Return: ndarray, float; len = bin_number. bin_labels, which are the average value for each bin.
        '''
        bin_sep = self.bin_seps(bin_number);
        bin_label = np.zeros(bin_number);
        for i in np.arange(bin_number - 1):
            bin_label[i+1] = (bin_sep[i+1] + bin_sep[i])/2.0;
        return bin_label;
    
    def fft_max_bin(self, spectrum):
        '''
        ---
        spectrum: 1darray, float. usually the radial spectrum.
        ---
        Return: tuple, (max_val, index): spectral value of the max bin and the bin index: 
        '''
        max_val = spectrum.max()
        index = np.where(spectrum == max_val)[0][0]
        return max_val, index
    
    def mean_bin(self, bins, vals):
        '''Compute the mean bin value for the spectrum.
        This is just a weighted average.
        ---
        bins: array-like. len(bins)=bin_number
        vals: array-like. len(vals)=bin_number. magnitude of each bin.
        '''
        return np.dot(vals,bins) / vals.sum()
    
    def k_mean(self, energy_spec):
        '''Compute the k_mean in the original way. usually takes more time.
        ---
        energy_spec: 3d array, float. The 3D energy spectrum.
        ---
        Return: k_mean, float.
        '''
        '''
        kspec = 0.0
        for xi in range(self.nDim):
            for yi in range(self.nDim):
                for zi in range(self.nDim):
                    k_mag = math.sqrt(self.freq[xi]**2 \
                                      + self.freq[yi]**2 \
                                      + self.freq[zi]**2)
                    kspec += k_mag * energy_spec[xi,yi,zi]
        '''
        gridx, gridy, gridz = np.meshgrid(self.freq, self.freq, self.freq)
        kspec = np.sum( energy_spec * np.sqrt( np.square(gridx) + np.square(gridy) + np.square(gridz) ) )
        return kspec / energy_spec.sum()
    
    def fit_powerlaw_two_parts(self, bin_labels, spectrum, 
                               small_k_samples = None,
                               large_k_samples = None):
        '''
        Fit two parts separated by the maximum in the spectrum.
        The zero mode will be neglected
        '''
        size = bin_labels.size
        max_val, max_index  = self.fft_max_bin(spectrum)
        if small_k_samples == None:
            first_bins = bin_labels[1 : max_index]
            first_spec = spectrum[1 : max_index]
        else:
            first_bins = bin_labels[1 : small_k_samples + 1]
            first_spec = spectrum[1 : small_k_samples + 1]
        
        if large_k_samples == None:
            second_bins = bin_labels[max_index + 1 : size]
            second_spec = spectrum[max_index + 1 : size]
        else:
            second_bins = bin_labels[size - large_k_samples : size]
            second_spec = spectrum[size - large_k_samples : size]
            
        log_first_bins = np.log(first_bins)
        log_second_bins = np.log(second_bins)
        log_first_spec = np.log(first_spec)
        log_second_spec = np.log(second_spec)
        
        linear = lambda x, a, b: a*x+b
        popt1, pcov1 = curve_fit(linear, log_first_bins, log_first_spec)
        popt2, pcov2 = curve_fit(linear, log_second_bins, log_second_spec)
        return popt1, popt2
    
    def fit_powerlaw(self, bin_labels, spectrum, 
                     sample_number = None,
                     fig = None):
        '''
        fit the ascending part with power law
        The zero mode will be neglected
        '''
        size = bin_labels.size
        max_val, max_index  = self.fft_max_bin(spectrum)
        if sample_number == None:
            first_bins = bin_labels[1 : max_index]
            first_spec = spectrum[1 : max_index]
        else:
            first_bins = bin_labels[1 : sample_number + 1]
            first_spec = spectrum[1 : sample_number + 1]
        log_first_bins = np.log(first_bins)
        log_first_spec = np.log(first_spec)
        linear = lambda x, a, b: a*x+b
        popt1, pcov1 = curve_fit(linear, log_first_bins, log_first_spec)
        fitx1 = np.linspace(bin_labels[0], bin_labels[max_index], 50)
        powE = lambda x, a, b: np.power(x, a)*np.exp(b)
        fity1 = powE(fitx1, *popt1)
        if(fig == None):
            fig = plt.figure(100,figsize=(6,4))
        plt.scatter(first_bins, first_spec, marker = 's')
        plt.plot(fitx1, fity1)
        return popt1
    
    def fit_powerexp(self, bin_labels, spectrum, fig = None):
        '''
        Fit the descending part with y=k*x^a*exp(-b*x)
        that is, lny = a*lnx-bx +lnk
        '''
        size = bin_labels.size
        max_val, max_index  = self.fft_max_bin(spectrum)
        bins = bin_labels[max_index + 1 : size]
        norm_spec = spectrum[max_index + 1 : size] / max_val
        powerexp = lambda x, a, b, k: k*np.power(x, a)*np.exp(-b*x)
        
        popt2, pcov2 = curve_fit(powerexp, bins, norm_spec)
        fitx2 = np.linspace(bin_labels[max_index + 1], bin_labels[-1], 50)
        fity2 = max_val*powerexp(fitx2, *popt2)
        if(fig == None):
            fig = plt.figure(100,figsize=(6,4))
        plt.scatter(bins, max_val*norm_spec, marker = 'D')
        plt.plot(fitx2, fity2)
        return popt2
    
        
    def k_space_meshgrid(self):
        mx, my, mz = np.meshgrid(self.freq, self.freq, self.freq, indexing = 'ij');
        self.mesh = np.array(list(zip(mx.ravel(), my.ravel(), mz.ravel())));
        del mx, my, mz;
    
    def build_kdtree(self):
        self.k_space_meshgrid();
        self.tree = spatial.cKDTree(self.mesh);
    
    '''auxillary functions'''
    #query all points larger than min_r (not included) but no larger than max_r (included).
    #that is (min_r, max_r]
    def query_shell_points(self, p, max_r, min_r):
        max_r_list = self.tree.query_ball_point(p, max_r);
        min_r_list = self.tree.query_ball_point(p, min_r);
        shell_points = [x for x in max_r_list if x not in min_r_list];
        return shell_points;    
    
    def max_freq_magnitude(self):
        '''
        Return: max magnitude of the frequency.
        assume a cubic 3D lattice
        '''
        max_freq = np.max(np.absolute(self.freq));
        return (self.DIM * max_freq**2)**0.5;
    
    def distance(self, p1, p2):
        return spatial.distance.euclidean(p1,p2);

    def distance_list(self, p1, p_list):
        dist = [];
        for p in p_list:
            dist.append(spatial.distance.euclidean(p1,p));
        return dist;
    

'''class of 2-point correlation function'''
'''not efficient enough'''
class Correlation:
    DIM = 3;
    data_dict = {'Bx': ['.Bx.datfast'],\
                 'By': ['.By.datfast'],\
                 'Bz': ['.Bz.datfast'],\
                 'EME': ['.EME.datfast'], \
                 'EMB': ['.EMB.datfast']}
    
    def __init__(self, file_path, run_id, lat_size, data_key, is_periodic = True):
        self.file_path = file_path;
        self.run_id = run_id;
        self.nDim = lat_size;
        self.periodic = is_periodic;
        data_param = Correlation.data_dict[data_key];
        self.x_data = np.fromfile(self.file_path + self.run_id + \
                                  data_param[0], \
                                  dtype=np.dtype('float32'));
        self.x_mean = self.x_data.mean();
        self.x_stdev = self.x_data.std();
        self.x_data = self.x_data.reshape((self.nDim,self.nDim,self.nDim)).transpose();
        #by default, y_data will be the same as x_data
        self.y_data = self.x_data;
        self.y_mean = self.x_mean;
        self.y_stdev = self.x_stdev;
        #real_range indicates the point range that is not ghost points
        #in a tuple (a,b), where a is included and b is not.
        self.real_range = (0, 0);
        if(is_periodic == False):
            full_nDim = self.nDim;
            self.real_range = (0, self.nDim);
        if(is_periodic == True):
            full_nDim = self.nDim * 3; #3 is more than needed. Here use 3 to do test  
            self.real_range = (self.nDim, 2*self.nDim);
        coord1d = np.arange(full_nDim);
        mx, my, mz = np.meshgrid(coord1d, coord1d, coord1d, indexing = 'ij');
        #construct a mesh that includes both real points and ghost points.
        #a periodic boundary will require ghost points.
        #For a given index n, the transform to the 3D coordinates is 
        # TODO
        self.mesh = np.array(list(zip(mx.ravel(), my.ravel(), mz.ravel())));
        return;
    
    #if one needs to compute correlation between two data sets
    #then a second data needs to be added.
    def add_second_data(self, second_data_key):
        data_param = Correlation.data_dict[second_data_key];
        self.y_data = np.fromfile(self.file_path + self.run_id + \
                                  data_param[0], \
                                  dtype=np.dtype('float32'));
        self.y_mean = self.y_data.mean();
        self.y_stdev = self.y_data.std();
        return;
    
    def build_kdtree(self):
        self.tree = spatial.cKDTree(self.mesh);
        return;
    
    def two_point_corr(self, bin_number):
        bin_sep = np.linspace(0, self.max_corr_distance(), bin_number + 1);
        corr_val = np.zeros(bin_number);
        for i in np.arange(bin_number):
            t1 = time.time();
            corr_val[i] = self.corr_in_one_bin(bin_sep[i+1], bin_sep[i]);
            t2 = time.time();
            print("bin # ", i, "; corr_val = ", corr_val[i]);
            print('time cost: ', t2-t1);
        '''add self-point contribution'''
        zero_corr = 0;
        for i in range(self.nDim):
            for j in range(self.nDim):
                for k in range(self.nDim):
                    p1 = self.x_data[i,j,k];
                    p2 = self.y_data[i,j,k];
                    zero_corr += (p1 - self.x_mean)*(p2 - self.y_mean);
        zero_corr /= self.x_stdev * self.y_stdev;
        corr_val[0] += zero_corr;
        return corr_val;
    
    def corr_in_one_bin(self, max_r, min_r):
        relevant_pairs = self.query_relevant_pairs(max_r, min_r);
        print('relevant pair size = ', len(relevant_pairs));
        corr = 0;
        for p in relevant_pairs:
            coords = self.point_pair_to_coords(p);
            rc1 = self.to_real_coordinates(coords[0]);
            rc2 = self.to_real_coordinates(coords[1]);
            dp1 = self.x_data[rc1[0], rc1[1], rc1[2]];
            dp2 = self.y_data[rc2[0], rc2[1], rc2[2]];
            corr += (dp1 - self.x_mean)*(dp2 - self.y_mean);
        corr = corr / ( self.x_stdev * self.y_stdev );
        return corr;
    
    def corr_bins(self, bin_number):
        '''The bin_sep are the starting and end points for each bin'''
        bin_sep = np.linspace(0, self.max_corr_distance(), bin_number + 1);
        bin_val = np.zeros(bin_number);
        for i in np.arange(bin_number):
            bin_val[i] = 0.5*(bin_sep[i]+bin_sep[i+1]);
        return bin_val;
        
    
    def query_relevant_pairs(self, max_r, min_r):
        max_pair_list = self.tree.query_pairs(max_r);
        min_pair_list = self.tree.query_pairs(min_r);
        shell_pairs = [x for x in max_pair_list if x not in min_pair_list];
        del max_pair_list, min_pair_list;
        relevant_pairs = shell_pairs;
        for p in shell_pairs:
            if(not self.is_relevant_pair(self.point_pair_to_coords(p))):
                relevant_pairs.remove(p);
        del shell_pairs;
        return relevant_pairs; #returns a set of tuples
    
    #sort all the pairs into different bin sections
    #For memory efficiency, the result will be stored in a 1d array.
    def sort_all_pairs(self, bin_number):
        bin_sep = np.linspace(0, self.max_corr_distance(), bin_number + 1);
    
    
    #Translate pair index to coordinates
    #return a tuple of coordintes pairs.
    def point_pair_to_coords(self, pair):
        first = self.mesh[pair[0]];
        second = self.mesh[pair[1]];
        pair_coords = (first, second);
        return pair_coords;
    
    #test if a pair of coordinates is relevant
    #a relevant pair is one in which at least one point is in the real range
    #if both points are in ghost region, then the pair is irrelevant.
    def is_relevant_pair(self, pair_coords):
        if(self.periodic == False):
            return True;
        first = pair_coords[0];
        second = pair_coords[1];
        ghost_first = False;
        ghost_second = False;
        for i in range(Correlation.DIM):
            if( first[i] <= self.real_range[0] or first[i] > self.real_range[1] ):
                ghost_first = True;
                break;
            if( second[i] <= self.real_range[0] or second[i] > self.real_range[1] ):
                ghost_second = True;
                break;
        if(ghost_first == True and ghost_second == True):
            return False;
        else:
            return True;
    
    def to_real_coordinates(self, coord):
        real_coords = np.ndarray(Correlation.DIM, dtype = np.uint64);
        for i in range(Correlation.DIM):
            real_coords[i] = coord[i] % self.nDim;
        return real_coords;
    
    def max_corr_distance(self):
        if(self.periodic == True):
            return self.nDim//2 - 0.5;
        else:
            return self.nDim - 1;
    
    def distance(self, p1, p2):
        return spatial.distance.euclidean(p1,p2);