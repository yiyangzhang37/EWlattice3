# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 09:07:50 2017

@author: YiyangZhang
"""

import numpy as np
import pandas as pd
import codecs
import os.path
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
from scipy.optimize import curve_fit


class ThermalData:
    eta = 6.0
    eta2 = eta*eta 
    
    def read_path(self, path_):
        self.path = path_
    
    def read_infos(self, path_):
        self.path = path_
        all_files = glob.glob( path_+ '/*.info')
        self.df_info = pd.DataFrame()
        #loop over all files
        for file_ in all_files:
            param_ = {}
            param_['file'] = file_.strip().split('\\')[-1]
            param_['ID'] = param_['file'][:-5]
            param_file_ = open(file_, 'r')
            for line_ in param_file_:
                key, val = line_.strip().split(':')
                param_[key.strip()] = val.strip()
            param_file_.close()
            df_ = pd.DataFrame([param_])
            self.df_info = pd.concat([self.df_info, df_])
        self.df_info.reset_index(inplace=True, drop=True)
        #convert to numerics
        for c in self.df_info.columns:
            try:
                self.df_info[c] = pd.to_numeric(self.df_info[c])
            except:
                pass
    def read_multi_infos(self, paths_):
        self.path = paths_[0]
        all_files = []
        for p in paths_:    
            all_files += glob.glob( p + '/*.info')
        self.df_info = pd.DataFrame()
        #loop over all files
        for file_ in all_files:
            param_ = {}
            param_['file'] = file_.strip().split('\\')[-1]
            param_['ID'] = param_['file'][:-5]
            param_file_ = open(file_, 'r')
            for line_ in param_file_:
                key, val = line_.strip().split(':')
                param_[key.strip()] = val.strip()
            param_file_.close()
            df_ = pd.DataFrame([param_])
            self.df_info = pd.concat([self.df_info, df_])
        self.df_info.reset_index(inplace=True, drop=True)
        #convert to numerics
        for c in self.df_info.columns:
            try:
                self.df_info[c] = pd.to_numeric(self.df_info[c])
            except:
                pass
            
    def add_info_source_xls(self, file_name, sheet_name_ = 0):
        source = pd.read_excel(os.path.join(self.path, file_name), \
                               sheet_name = sheet_name_).transpose()
        source.reset_index(inplace=True)
        source.rename(columns = source.iloc[0], inplace = True)
        source.drop(index = 0, inplace = True)
        source.reset_index(inplace=True, drop = True)
        for s in source.columns:
            try:
                source[s] = pd.to_numeric(source[s])
            except:
                pass
        self.df_info = pd.merge(self.df_info, source, on = 'ID', suffixes = ('','_dup'))
        return source
                
    
    def output_info_to_html(self, file_name = 'list_info.html', is_output_short = False):
        self.df_info.transpose().to_html(file_name)
        if(is_output_short):
            short = self.df_info.copy()
            for c in short.columns:
                if(len(short[c].unique()) == 1):
                    del short[c]
            short.transpose().to_html('s_' + file_name)
    
    def output_info_to_csv(self, file_name = 'list_info.csv'):
        self.df_info.to_csv(file_name)
    
    
    def assign_pt_end_crit(self, Higgs2_, duration_):
        '''
        The parameter Higgs2_ is taken as the normalized value
        The phase transition end criteria is such that the mean value of min(|Phi|^2) 
        meets Higgs2_ in a duration_ time.
        '''
        self.end_Higgs2_crit = Higgs2_ * ThermalData.eta2
        self.end_dura_crit = duration_
    def scan_pt_end_points(self):
        kwargs = {'EndPoint(TimeStep)' : np.nan}
        self.df_info = self.df_info.assign(**kwargs)
        for i in self.df_info['ID']:
            self.read_single_run_data(i)
            end_p = self.find_pt_end_point(i)
            info_index = self.find_curr_index()
            if np.isnan(end_p):
                self.df_info.at[info_index,'EndPoint(TimeStep)'] = len(self.one_data.index) - 1
            else:
                self.df_info.at[info_index,'EndPoint(TimeStep)'] = end_p

    
    def scan_run_data_quantities(self):
        '''
        Scan the whole folder (all runs) and compute quantities that can be obtained 
        from the data_table data
        Total bubble
        Actual TimeSteps
        '''
        kwargs = {'TotalBubbles': np.nan, 
                  'ActualTimeSteps': np.nan,
                  'PercolationSize': np.nan,
                  'PercolationTime': np.nan, #end of bubble nucleation
                  'BubbleNumberDensity': np.nan,
                  'GaugeEnergyRatio': np.nan,# (E_U1+E_SU2)/E_TOTAL
                  'EMEnergyRatio': np.nan, # (E_emb+E_eme) / E_TOTAL
                  'MagneticEnergyRatio': np.nan # (E_emb) / E_TOTAL
                  } 
        self.df_info = self.df_info.assign(**kwargs)
        for i in self.df_info['ID']:
            self.read_single_run_data(i)
            info_index = self.find_curr_index()
            end_p = self.find_pt_end_point(i)
            self.df_info.at[info_index, 'TotalBubbles'] = self.compute_total_bubbles()
            self.df_info.at[info_index, 'ActualTimeSteps'] = len(self.one_data.index)    
            try:
                self.df_info.at[info_index, 'PercolationTime'] = \
                    self.one_data.index[self.one_data['new_bubble']>0][-1]
            except IndexError:
                self.df_info.at[info_index, 'PercolationTime'] = 0
        
        for i in self.df_info['ID']:
            #compute the percolation size
            info_index = self.df_info['ID'][self.df_info['ID']==i].index[0]
            dx = self.df_info['DX'].iloc[info_index]
            vol = self.df_info['LatticeSize(X)'].iloc[info_index] * \
                self.df_info['LatticeSize(Y)'].iloc[info_index] * \
                self.df_info['LatticeSize(Z)'].iloc[info_index] * dx**3
            self.df_info.loc[info_index, 'PercolationSize'] = \
            (vol / self.df_info.loc[info_index, 'TotalBubbles'])**(1./3.)
            #compute bubble number density
            self.df_info.loc[info_index, 'BubbleNumberDensity'] = \
                self.df_info.loc[info_index, 'TotalBubbles'] / vol
            #compute GaugeEnergyRatio
            self.df_info.loc[info_index, 'GaugeEnergyRatio'] = \
                (self.df_info.loc[info_index, 'U1Energy'] + self.df_info.loc[info_index, 'SU2Energy'] )\
                / self.df_info.loc[info_index, 'TotalEnergy']
            #compute EMEnergyRatio
            self.df_info.loc[info_index, 'EMEnergyRatio'] = \
                (self.df_info.loc[info_index, 'MagneticEnergy'] + self.df_info.loc[info_index, 'ElectricEnergy'] )\
                / self.df_info.loc[info_index, 'TotalEnergy']
            self.df_info.loc[info_index, 'MagneticEnergyRatio'] = \
                (self.df_info.loc[info_index, 'MagneticEnergy'])\
                / self.df_info.loc[info_index, 'TotalEnergy']
            
   
    def scan_pt_end_quantities(self, quantities):
        '''
        find the values of the quantities at the pt_end_point step,
        and store them in the df_info DataFrame
        '''
        for q in quantities:
            kwargs = {q : np.nan}
            self.df_info = self.df_info.assign(**kwargs)
        
        for i in self.df_info['ID']:
            self.read_single_run_data(i)
            info_index = self.find_curr_index()
            end_index = self.df_info['EndPoint(TimeStep)'][info_index]
            for q in quantities:
                if(np.isnan(end_index)):
                    val = np.nan
                else:
                    val = self.one_data[quantities[q]][end_index]
                #self.df_info.set_value(info_index, q, val)
                #deprecated
                self.df_info.loc[info_index, q] = val
        del kwargs, info_index, end_index, val
     
    def fit(self, x_col, y_col, fixed_param, fixed_val, \
                 logx_ = False, logy_ = False, scatter_only_ = False):
        plot_indices = np.where(self.df_info[fixed_param] == fixed_val)
        if(logx_):
            x_vals = np.log10(self.df_info.iloc[plot_indices][x_col]).values
        else:
            x_vals = self.df_info.iloc[plot_indices][x_col].values
        if(logy_):
            y_vals = np.log10(self.df_info.iloc[plot_indices][y_col]).values
        else:    
            y_vals = self.df_info.iloc[plot_indices][y_col].values
        plt.scatter(x_vals, y_vals)
        
        popt, pcov = curve_fit(self.linearfit, x_vals, y_vals)
        if(not scatter_only_):
            fit_x = np.linspace(x_vals.min(), x_vals.max(), 100)
            fit_y = self.linearfit(fit_x, *popt)
            plt.plot(fit_x, fit_y, label = str(fixed_val))
            
        plt.title( fixed_param + ' fixed' )
        if(logx_):
            plt.xlabel( 'Log('+x_col+')' )
        else:
            plt.xlabel( x_col )
        if(logy_):
            plt.ylabel( 'Log('+y_col+')' )
        else:
            plt.ylabel( y_col )
        return popt
    
    def fit_general(self, x_col, y_col, \
                 logx_ = False, logy_ = False, scatter_only_ = False):
        if(logx_):
            x_vals = np.log10(self.df_info[x_col]).values
        else:
            x_vals = self.df_info[x_col].values
        if(logy_):
            y_vals = np.log10(self.df_info[y_col]).values
        else:    
            y_vals = self.df_info[y_col].values
        plt.scatter(x_vals, y_vals)
        
        popt, pcov = curve_fit(self.linearfit, x_vals, y_vals)
        if(not scatter_only_):
            fit_x = np.linspace(x_vals.min(), x_vals.max(), 100)
            fit_y = self.linearfit(fit_x, *popt)
            plt.plot(fit_x, fit_y)
            
        plt.title( y_col + ' vs. ' + x_col)
        if(logx_):
            plt.xlabel( 'Log('+x_col+')' )
        else:
            plt.xlabel( x_col )
        if(logy_):
            plt.ylabel( 'Log('+y_col+')' )
        else:
            plt.ylabel( y_col )
        return popt
    
                        
#one run data table functions    
    def read_single_run_data(self, run_id_):
        self.curr_id = run_id_
        data_table_name = self.path + run_id_ + '.txt';
        self.one_data = pd.read_csv(data_table_name) \
            .dropna(axis = 'columns', how ='any')
    
    def find_curr_index(self):
        '''
        find the curr_id in the df_info DataFrame.
        This only works with read_single_run_data()
        Return: the index.
        '''
        return self.df_info['file'][self.df_info['file'] == self.curr_id + '.info'].index[0]
    
    def add_new_value(self, new_col_, val_, index_):
        '''
        add new value to the df_info DataFrame.
        a new column will be constructed if it does not exist.
        '''
        if(new_col_ not in self.df_info.columns):
            kwargs = {new_col_ : np.nan}
            self.df_info = self.df_info.assign(**kwargs)
            del kwargs
        self.df_info.set_value(index_, new_col_, val_)          
        
    def compute_total_bubbles(self):
        return self.one_data['new_bubble'].sum()
        
    def find_pt_end_point(self, run_id_):
        if(self.curr_id != run_id_):
            print('Different runs.')
            return
        data_size = self.one_data['minHiggs2'].size
        mean = np.zeros(data_size)
        for i in range(self.end_dura_crit, data_size):
            mean[i] = self.one_data['minHiggs2'][i - self.end_dura_crit : i].mean();
        end_p = np.nonzero(mean > self.end_Higgs2_crit)
        if(end_p[0].size != 0):
            return end_p[0][0]
        else:
            return np.nan
        
    def plot_one_run_data(self, run_id_, is_save_ = False):
        self.read_single_run_data(run_id_)
        #self.one_data = self.one_data.iloc[0:100]
        df_normalize = self.one_data / self.one_data.max();
        
        fig = plt.figure(1,figsize=(15,12))
        gs = gridspec.GridSpec(3,3)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])
        ax3 = plt.subplot(gs[0,2])
        ax4 = plt.subplot(gs[1,0])
        ax5 = plt.subplot(gs[1,1])
        ax6 = plt.subplot(gs[1,2])      
        ax7 = plt.subplot(gs[2,0])
        ax8 = plt.subplot(gs[2,1])
        ax9 = plt.subplot(gs[2,2])
        
        self.one_data[['energy','E_kinetic','E_grad','E_pot','E_U1','E_SU2']].plot(ax=ax1, logy=True)
        (self.one_data['EM_energy']/self.one_data['EM_electricEnergy']).plot(ax=ax2, title='EMB/EME')
        self.one_data[['CS_number','W_Higgs']].plot(title='windings', ax=ax3)
        self.one_data[['Higgs_mag2']].plot(title='|PHI|^2', ax=ax4)
        #modify EM energy
        self.one_data[['EM_energy','EM_electricEnergy']].plot(ax=ax6,logy=True,ylim = (10**-5,10**6))
        df_normalize[['new_bubble','E_pot']].plot(ax=ax9)
        self.one_data[['Hamiltonian_GC']].plot(ax=ax5)
        self.one_data[['helicity']].plot(title='helicity', ax=ax7)
        self.one_data['minHiggs2'].plot(ax=ax8, title='minHiggs2')

        try:
            self.one_data[['EM_electricEnergy','EM_electricEnergy_noScalar','EM_electricEnergy_pureScalar']] \
            .plot(ax=ax7, logy=True, title='electric energy',ylim=(10**-5,10**6))
            self.one_data[['EM_energy','EM_energy_noScalar','EM_energy_pureScalar']] \
            .plot(ax=ax8, logy=True, title='magnetic energy',ylim=(10**-5,10**6))
        except:
            pass
        plt.show()
        if(is_save_):
            fig.savefig(run_id_ + '.png')
    
    def plot_all_runs(self, is_save_ = True):
        for i in self.df_info['ID']:
            self.plot_one_run_data(i, is_save_)
    
    #Fitting functions
    def linearfit(self, x, a, b):
        return a*x+b  

             
        
    