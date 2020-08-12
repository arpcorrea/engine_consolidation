# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:51:32 2020

@author: andre
"""

import numpy as np

class BoundaryCondition():
    
    def __init__(self, deck,geometry,meshes):
        self.deck = deck
        self.geometry = geometry
        self.meshes=meshes
        self.set_temperatures()
        self.set_conductivity()
        self.set_density()
        self.set_specific_heat()
        self.set_diffusivity()
        
    def set_temperatures(self):    
        
        T = np.zeros((self.meshes.ny, self.meshes.nx))
        T[0:,0:] = self.deck.doc["Boundary Condition"]["Initial Temperature"]
        self.T=T.copy()
        self.T0=T.copy()
        
        self.Troom=self.deck.doc["Boundary Condition"]["Room Temperature"]
        
    def set_conductivity(self):
        
        Kx = np.zeros((self.meshes.ny, self.meshes.nx))
        Kx[0:,0:] = self.deck.doc["Materials"]["Aluminium"]["Thermal Conductivity X"]
        self.Kx=Kx.copy()
      
        
        Ky = np.zeros((self.meshes.ny, self.meshes.nx))
        Ky[0:,0:] = self.deck.doc["Materials"]["Aluminium"]["Thermal Conductivity Y"]
        self.Ky=Ky.copy()

        
        
    def set_density(self):
        
        Rho = np.zeros((self.meshes.ny, self.meshes.nx))
        Rho[0:,0:] = self.deck.doc["Materials"]["Aluminium"]["Density"]
        self.Rho=Rho.copy()
        
    def set_specific_heat(self):
        
        Cp = np.zeros((self.meshes.ny, self.meshes.nx))
        Cp[0:,0:] = self.deck.doc["Materials"]["Aluminium"]["Specific Heat Capacity"]
        self.Cp=Cp.copy()
        
        

    def set_diffusivity(self):
        Dx=np.zeros((np.shape(self.T)))
        Dx[0:,0:]=self.Kx[0:,0:]/(self.Rho[0:,0:]*self.Cp[0:,0:])
        self.Dx=Dx.copy()
        
        Dy=np.zeros((np.shape(self.T)))
        Dy[0:,0:]=self.Ky[0:,0:]/(self.Rho[0:,0:]*self.Cp[0:,0:])
        self.Dy=Dy.copy()
    
    


    

     
        
        
        