# -*- coding: utf-8 -*-
import numpy as np

        
        
        
        
        
        
class MeshOnePlate():

    def __init__(self, deck,geometry):
        self.deck = deck
        self.geometry=geometry
        self.set_mesh_grid() 
        self.set_dx()
        self.set_dy()
        self.set_Element_Position()
        
        
    def set_mesh_grid(self):
        self.nx = int(self.deck.doc["Simulation"]["Number of Elements X"])        
        self.ny = int(self.deck.doc["Simulation"]["Number of Elements Y"]) 
  
        self.dx = self.geometry.Lx/self.nx
        self.dy = self.geometry.Ly/self.ny
        
    
    def set_dx(self):
        Mdx = np.zeros((self.ny, self.nx))
        Mdx[0:,0:]=self.dx
        self.Mdx=Mdx.copy()
        
        
    def set_dy(self):
        Mdy = np.zeros((self.ny, self.nx))
        Mdy[0:,0:]=self.dy
        self.Mdy=Mdy.copy()
        
    def set_Element_Position(self):
        
        ElementXPosition=  np.arange(self.nx)*self.dx+self.dx  
        self.ElementXPosition=ElementXPosition.copy()
        
        ElementYPosition= np.arange(self.ny)*self.dy+self.dy      
        self.ElementYPosition=ElementYPosition.copy()
        
        self.Xposition=self.ElementXPosition
        self.Yposition=self.ElementYPosition
