# -*- coding: utf-8 -*-
import numpy as np

class Plate():

    def __init__(self, deck,geometry):
        self.deck = deck
        self.geometry=geometry
        self.set_mesh_grid() 
        self.set_temperatures()
        self.set_conductivity()
        self.set_density()
        self.set_specific_heat()
        self.set_diffusivity()
        self.set_heat_density()
        
        
        
    def set_mesh_grid(self):
        self.nx = int(self.deck.doc["Simulation"]["Number of Elements X"])
        self.ny = int(self.deck.doc["Simulation"]["Number of Elements Y"])
        self.dx=self.geometry.Lx/self.nx
        self.dy=self.geometry.Ly/self.ny
        X, Y = np.meshgrid(np.arange(0, self.ny), np.arange(0, self.nx))
        X=X[1,:].copy()
        Y=Y[:,1].copy()
        self.X = X
        self.Y = Y

        

    
    def set_temperatures(self):
        T = np.zeros((self.ny, self.nx))        
        T[0:, 0:] = self.deck.doc["Materials"]["Aluminium"]["Initial Temperature"] # Set array size and set the interior value with Tini
        self.T = T.copy()
        self.T0=T.copy()
        
      
    def set_conductivity(self):
        KtotalX= np.zeros((self.ny, self.nx)) 
        KtotalX[0:, 0:] = self.deck.doc["Materials"]["Aluminium"]["Thermal Conductivity X"]
        self.KtotalX=KtotalX
                                                                                         
        KtotalY= np.zeros((self.ny, self.nx)) 
        KtotalY[0:, 0:] = self.deck.doc["Materials"]["Aluminium"]["Thermal Conductivity Y"]                                                                                    
        self.KtotalY=KtotalY         

    def set_density(self):
        RhoTotal= np.zeros((self.ny, self.nx)) 
        RhoTotal[0:, 0:] = self.deck.doc["Materials"]["Aluminium"]["Density"]                                                                                   
        self.RhoTotal=RhoTotal  

    def set_specific_heat(self):
        CpTotal= np.zeros((self.ny, self.nx)) 
        CpTotal[0:, 0:] = self.deck.doc["Materials"]["Aluminium"]["Cp"]                                                                                
        self.CpTotal=CpTotal  

    def set_diffusivity(self):                                                                                  
        DiffTotalX = np.zeros((self.ny, self.nx)) 
        DiffTotalX[0:,0:]=self.KtotalX[0:,0:]/(self.RhoTotal[0:,0:]*self.CpTotal[0:,0:])
        self.DiffTotalX = DiffTotalX.copy()
        
        DiffTotalY = np.zeros((self.ny, self.nx)) 
        DiffTotalY[0:,0:]=self.KtotalY[0:,0:]/(self.RhoTotal[0:,0:]*self.CpTotal[0:,0:])
        self.DiffTotalY = DiffTotalY.copy()


    def set_heat_density(self):       
        self.q=float(self.deck.doc["Boundary Conditions"]["Power"])
        Q=np.zeros((self.ny, self.nx))
        Q[0:, 0] = self.q
        self.Q=Q.copy()
        
        
        
        
        
