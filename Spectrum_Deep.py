# -*- coding: utf-8 -*-
"""
Created on Sun Jun  1 10:47:26 2025

@author: VINICIUSSALES
"""

import tkinter as tk
from tkinter import filedialog
import pandas as pd 
import pysptools.spectro as spectro
import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sns
from scipy.signal import savgol_filter


class Retruturation_Spectrum():
  def __init__(self,names,X,wl,Afloramento):
    self.names = names
    self.X = X
    self.data_ = {}
    self.lista = []
    self.df =[]

    for cont in range(len(names)):

        bands = FeaturesConvexHullQuotient(X[cont,:], wl, baseline = 1).get_information()
        # self.data_[names[cont]]={dict_['abs_wvl']:dict_ for dict_ in bands if dict_['abs_wvl'] is not None and dict_['abs_depth']}
        
        # # print(f" Espectro {cont}º ->{names[cont]}------------- ABS ->{list(self.data_[names[cont]].keys())}")
        # self.lista.append(np.array(list(self.data_[names[cont]].keys())))
        
        dict_={str(dict_['abs_wvl']):dict_['abs_depth'] for dict_ in bands if dict_['abs_wvl'] is not None and dict_['abs_depth']<0.95}
        dict_['Afloramento'] = Afloramento[cont]
        self.data_[names[cont]]=dict_
        
        
        
  def return_(self):
      
      return self.data_

  def ocorrencia_total(self, p= 0):
      p = p/100
      #p -> porcentagem de ocorrencia
      flattened_array = np.concatenate(self.lista)
      unique_values = np.unique(flattened_array,return_counts = True)
      if p >= 1:
        return -1
      position = np.sort(np.where(unique_values[1]>len(self.lista)-len(self.lista)* p)[0])

      for cont in list(position):
        print(f'WVL ----> {unique_values[0][cont]}   Frequência de Ocorrência ----> {unique_values[1][cont]}')


  def filtragem_por_band(self,band):

    list_cont = []
    index =[]
    try:
        for id in self.data_.keys():
          try:
             list_cont.append(self.data_[id][band])
             index.append(id)
          except:
            continue
        self.df =pd.DataFrame(list_cont, index = index)

        return   self.df.drop(columns = ['seq', 'id', 'state'])
    except:
      return -1

  def expand_column(self,df):

    columns_to_expand = ['hx', 'hy', 'FWHM_x', 'FWHM_y']
    for col_name in columns_to_expand:
      if df[col_name].apply(lambda x: isinstance(x, (list, tuple))).all():
          expanded_cols = pd.DataFrame(df[col_name].to_list(), index=df.index)
          expanded_cols.columns = [f"{col_name}_{i+1}" for i in range(expanded_cols.shape[1])]
          df = df.drop(columns=[col_name]).join(expanded_cols)

    return df

  def plot_graphics(self,df):

    for name in df.index:
      plt.figure(figsize=(10, 6))
      plt.plot(df.loc[name]['wvl'], df.loc[name]['crs'], label='Espectro')
      x = (df.loc[name]['abs_wvl'],df.loc[name]['abs_wvl'])
      y = (df.loc[name]['abs_depth'],1)
      plt.plot(x, y, 'r',color='r', linestyle='--', label='Pico de Absorção  '+ str(round(df.loc[name]['abs_wvl'],4)))
      x=(df.loc[name]['FWHM_x_1'],df.loc[name]['FWHM_x_2'])
      y=(df.loc[name]['FWHM_y_1'],df.loc[name]['FWHM_y_2'])
      plt.plot(x,y,color='g', linestyle='--',label='FWHM  '+str(round(df.loc[name]['FWHM_delta'],4)))
      plt.fill_between(np.asarray(df.loc[name]['wvl']),1, df.loc[name]['crs'], alpha=0.2, label='Área sob a curva  '+str(round(df.loc[name]['area'],4)))
      plt.plot([], [], ' ', label= 'Slope'+str(round(df.loc[name]['cslope'],4)))
      plt.title('Espectro com Identificação de Parâmetros '+ str(name))
      plt.xlabel('Comprimento de Onda (nm)')
      plt.grid()
      plt.ylabel('Intensidade')
      plt.legend()
      plt.show()
      
      
class FeaturesConvexHullQuotient(spectro.SpectrumConvexHullQuotient):
    """
    Remove the convex-hull of the signal by hull quotient and auto-extract features.
    A baseline can be applied to avoid non-significant features.

    Reference:
        Kokaly F. Raymond, PRISM: Processing Routines in IDL for Spectroscopic
        Measurements (Installation Manual and User's Guide, Version 1.0),
        U.S. Geological Survey, Reston, Virginia: 2011.
    """

    def __init__(self, spectrum, wvl, startContinuum=None, stopContinuum=None, baseline=0.98, normalize=False):
        if startContinuum is not None and stopContinuum is not None:
            start = np.searchsorted(wvl, startContinuum, side='left')
            stop = np.searchsorted(wvl, stopContinuum, side='right')
            spectrum = spectrum[start:stop]
            wvl = wvl[start:stop]

        super().__init__(list(spectrum), list(wvl), normalize)
        self.baseline = baseline
        self.features = []
        self.features_all = []
        self._extract_features()
        self._clean_baseline()

    def _extract_features(self):
        for feat_no in range(len(self.hx) - 1):
            start_idx = np.where(self.wvl == self.hx[feat_no])[0][0]
            stop_idx = np.where(self.wvl == self.hx[feat_no + 1])[0][0] + 1
            feature = {
                'seq': feat_no,
                'id': None,
                'state': None,
                'spectrum': self.spectrum[start_idx:stop_idx],
                'wvl': self.wvl[start_idx:stop_idx],
                'crs': self.crs[start_idx:stop_idx],
                'hx': self.hx[feat_no:feat_no+2],
                'hy': self.hy[feat_no:feat_no+2],
                'cstart_wvl': None,
                'cstop_wvl': None,
                'abs_wvl': None,
                'abs_depth': None,
                'area': None,
                'cslope': None,
                'FWHM_x': None,
                'FWHM_y': None,
                'FWHM_delta': None,
            }
            self.features_all.append(feature)

    def _clean_baseline(self):
        id_counter = 1
        for feat in self.features_all:
            if np.min(feat['crs']) < self.baseline:
                feat['state'] = 'keep'
                feat['id'] = id_counter
                id_counter += 1
                self._add_stats(feat)
                self.features.append(feat)
            else:
                feat['state'] = 'reject'

    def _add_stats(self, feat):
        feat['area'] = self._calculate_area(feat['crs'])
        feat['cstart_wvl'] = feat['wvl'][0]
        feat['cstop_wvl'] = feat['wvl'][-1]
        feat['abs_wvl'] = feat['wvl'][np.argmin(feat['crs'])]
        feat['abs_depth'] = np.min(feat['crs'])
        feat['cslope'] = (feat['hy'][1] - feat['hy'][0]) / (feat['hx'][1] - feat['hx'][0])
        feat['FWHM_x'], feat['FWHM_y'], feat['FWHM_delta'] = self._calculate_FWHM(feat)

    def _calculate_area(self, y):
        yy = np.abs(np.array(y) - 1)
        deltax = self.wvl[1] - self.wvl[0]
        return np.trapz(yy, dx=deltax)

    def _calculate_FWHM(self, feat):
        crs = np.array(feat['crs'])
        wvl = np.array(feat['wvl'])
        depth = np.min(crs)
        half_max = depth + (1 - depth) / 2
        left_indices = np.where(crs <= half_max)[0]
        right_indices = np.where(crs <= half_max)[0]
        if len(left_indices) == 0 or len(right_indices) == 0:
            return (None, None), (None, None), None
        left_wvl = wvl[left_indices[0]]
        right_wvl = wvl[right_indices[-1]]
        return (left_wvl, right_wvl), (half_max, half_max), right_wvl - left_wvl

    def get_information(self):
      return self.features_all
  
    

class Picos_Abs():
    
    def __init__(self):
        
        print('Import Espectros')
        X_normal = pd.read_excel(self.select_file(),sheet_name = 1,header = 0, index_col = 0 )
        X_Afloramento = list(X_normal.iloc[:,0])
        X_normal=X_normal.iloc[:,1:]
        names = list(X_normal.index.values)
        wl = list(X_normal.columns.values.astype('float'))
        X_cont_remove = self.continue_remove(X_normal.values,wl,names)
        
        
        self.plot_array(np.asarray(wl), X_normal)
        self.plot_array(np.asarray(wl), X_cont_remove)
        self.depth = Retruturation_Spectrum(names,X_normal.values,wl,X_Afloramento)
        
        
    def return_p(self):
        dados =  self.depth.return_()
        
        list_=[]
        for chave_amostra, dict_ in dados.items():
            
            afloramento = dict_.get("Afloramento")
            bandas = [float(k) for k, v in dict_.items() if k != "Afloramento"]
            
            list_.append({'nome_amostra':"Amostra_"+str(chave_amostra),
                          'bandas_detectadas':bandas,
                          'nome_afloramento':afloramento})
        return list_
        
    def select_file(self):
                    root = tk.Tk()
                    root.withdraw()
                                       
                    return filedialog.askopenfilename(parent=None, 
                                                              title="Select file", 
                                                              filetypes=[("excel",'*.xlsx')])
                
    def continue_remove(self,data,wvl,names):
            list_ =[]
            dataframe ={}
        
            for count in range(data.shape[0]):
        
              schq = spectro.SpectrumConvexHullQuotient(list(data[count]), list(wvl))
              dataframe[names[count]] = schq.get_continuum_removed_spectrum()
        
            return pd.DataFrame(dataframe, index = list(wvl)).T
    
    def plot_spectral_dataframe(self,df):
          plt.figure(figsize=(20, 16))
          for index in range(len(df.index)):
            plt.plot(df.columns.values.astype('float'),df.iloc[index,:],label = df.index[index])
        
          plt.title('Assinaturas Espectrais de Múltiplas Amostras')
          plt.xlabel('Comprimento de Onda (nm)')
          plt.ylabel('Reflectância')
          plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.05))
          plt.grid(True,linewidth=0.5)
          plt.show(block=False)
          plt.pause(0.1)
    
    def plot_array(self,wvl, X):
        
        plt.figure(figsize=(12,8))
        plt.plot(wvl, X.T, linewidth=1)
        plt.xlabel('Wavelength (nm)', fontsize=20)
        plt.ylabel('Reflectance', fontsize=20)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.grid(True,linewidth=0.5)
        plt.show(block=False)
        plt.pause(0.1)
    
