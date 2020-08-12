# -*- coding: utf-8 -*-
import numpy as np


class SolvesOnePlates:
    
    def __init__(self, deck,meshes,BC,model_HT,plots):
        self.deck = deck
        self.meshes=meshes
        self.BC=BC
        self.model_HT=model_HT
        self.plots = plots
        self.do_solver()
        
            
        
        
        
    def do_solver(self):
        for m in range(int(self.deck.doc["Simulation"]["Number Time Steps"])):
# -------------- CALCULATE TEMPERATURE FOR EACH STEP INCREMENT----------             
            # self.A=self.model_HT.convection(self.meshes.T0, self.meshes.T, self.meshes.DiffTotalX, self.meshes.DiffTotalY)
            self.BC.T0, self.BC.T = self.model_HT.do_convection(self.BC.T0, self.BC.T, self.BC.Dx, self.BC.Dy)            
            self.BC.T0, self.BC.T = self.model_HT.do_timestep(self.BC.T0, self.BC.T, self.BC.Dx, self.BC.Dy)
            
            # -------------- FORCE TEMPERATURE AT THE INTERFACE: ISOTHERMAL CONDITION----------             
# -------------- UPDATE T0----------             
            self.BC.T0=self.BC.T.copy()

# -------------- DO PLOT ACCORDING TO THE SELECTED INTERVAL----------             
            if m in self.plots.mfig:
                
                self.plots.update_T(self.BC.T)       
                self.plots.do_plots(m)    
        self.plots.do_animation()

            
     