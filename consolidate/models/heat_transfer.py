import numpy as np

class HeatTransfer:

    def __init__(self, deck,meshes):
        self.dt = float(deck.doc["Simulation"]["Time Step"])
        self.dx2 = meshes.dx*meshes.dx
        self.dy2 = meshes.dy*meshes.dy
        self.ny=meshes.ny
        self.nx=meshes.nx
        self.dx=meshes.dx
        self.dy=meshes.dy
        self.rho=meshes.RhoTotal
        self.cp=meshes.CpTotal
        self.kx=meshes.KtotalX
        self.ky=meshes.KtotalY
        self.h=float(deck.doc["Boundary Conditions"]["convection coefficient"])
        self.Troom=float(deck.doc["Boundary Conditions"]["Room Temperature"])
        self.q=float(deck.doc["Boundary Conditions"]["Power"])
# -------------- BEGIN HEAT TRANSFER CALCULATION---------- 
    def do_timestep(self, u0, u, Diffx, Diffy,Q):
        # Propagate with forward-difference in time, central-difference in space
        

        Q=np.zeros((self.ny, self.nx))
        Q[0,:]=self.h*(self.Troom - u[0,:])
        Q[-1,:]=self.h*(self.Troom - u[-1,:])
        Q[:,-1]=self.h*(self.Troom - u[:,-1])     
        Q[:,0]=self.q
        self.Q=Q
        
        
        
        value=np.zeros((self.ny, self.nx))
        value[0,:] = -Q[0,:]/self.ky[0,:]
        value[-1,:]= -Q[-1,:]/self.ky[-1,:]
        value[:,-1]= -Q[:,-1]/self.kx[:,-1]
        value[:,0]= -Q[:,0]/self.kx[:,0]
        self.value=value
        
        
        Uout=np.zeros((self.ny,self.nx))
        Uout[0,:] =u[0,:]-2*self.dy*value[0,:]
        Uout[:,-1]=u[:,-1]-2*self.dx*value[:,-1]
        Uout[-1,:]=u[-1,:]-2*self.dy*value[-1,:]
        Uout[:,0]= u[:,0]-2*self.dx*value[:,0]
        self.Uout=Uout
        
        
        
       
        
    
       
        
        
        u[0,1:-1]  = u0[0,1:-1]+ Diffy[0,1:-1]*self.dt*((u0[1,1:-1]-2*u[0,1:-1]+ Uout[0,1:-1])/self.dy2) + Diffx[0,1:-1]*self.dt*((u0[0, 0:-2]-2*u0[0,1:-1]+u0[0,2:])/self.dx2) 
        u[-1,1:-1] = u0[-1,1:-1]+ Diffy[-1,1:-1]*self.dt*((u0[-2,1:-1]-2*u[-1,1:-1]+ Uout[-1,1:-1])/self.dy2) + Diffx[-1,1:-1]*self.dt*((u0[-1, 0:-2]-2*u0[-1,1:-1]+u0[-1,2:])/self.dx2)       
        u[1:-1,-1] = u0[1:-1,-1]+ Diffy[1:-1,-1]*self.dt* ((u0[2:,-1]-2*u[1:-1,-1]+ Uout[:-2,-1])/self.dy2) + Diffx[1:-1,-1]   * self.dt * ((u0[1:-1, -2 ]   -2*u0[1:-1,-1]    +u0[1:-1,-1])/self.dx2) 
        u[1:-1,0]  = u0[1:-1,0] + Diffy[1:-1,0] *self.dt* ((u0[2:,0 ]-2*u[1:-1,0]+  Uout[:-2, 0])/self.dy2) + Diffx[1:-1, 0]    *self.dt*((u0[1:-1, 1 ]       -2*u0[1:-1,0]     +u0[1:-1,0])/self.dx2) 


        u[0,0]   = u0[0,0]   +Diffy[0,0]   *self.dt*((u0[1,0]   -2*u[0,0]   + Uout[0,0])/self.dy2)   +Diffx[0,0]  *self.dt*((u0[0, 1]   -2*u0[0,0]   +Uout[0,0])/self.dx2)
        u[-1,0]  = u0[-1,0]  +Diffy[-1,0]  *self.dt*((u0[-2,0]  -2*u[-1,0]  + Uout[-1,0])/self.dy2)  +Diffx[-1,0] *self.dt*((u0[-1, 1]  -2*u0[-1,0]  +Uout[-1,0])/self.dx2)  
        u[0,-1]  = u0[0,-1]  +Diffy[0,-1]  *self.dt*((u0[1,-1]  -2*u[0,-1]  + Uout[0,-1])/self.dy2)  +Diffx[0,-1] *self.dt*((u0[0, -2]  -2*u0[0,-1]  +Uout[0,-1])/self.dx2)
        u[-1,-1] = u0[-1,-1] +Diffy[-1,-1] *self.dt*((u0[-2,-1] -2*u[-1,-1] + Uout[-1,-1])/self.dy2) +Diffx[-1,-1]*self.dt*((u0[-1,-2]  -2*u0[-1,-1] +Uout[-1,-1])/self.dx2)

      
        
        
        
        u[1:-1, 1:-1] = u0[1:-1, 1:-1] + Diffy[1:-1, 1:-1]* self.dt * ((u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/self.dy2 ) + Diffx[1:-1, 1:-1]* self.dt * ( (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/self.dx2 ) + self.dt*Q[1:-1,1:-1]/(self.cp[1:-1,1:-1]*self.rho[1:-1,1:-1])

        u0 = u.copy()
        
        return u0, u
    
# -------------- END HEAT TRANSFER CALCULATION---------- 
        
    
    # def convection(self, u0, u, Diffx, Diffy):
    #     Q=np.zeros((self.ny, self.nx))
    #     Q[0,:]=self.h*(self.Troom - u[0,:])
    #     Q[-1,:]=self.h*(self.Troom - u[-1,:])
    #     Q[:,-1]=self.h*(self.Troom - u[:,-1])
        
    #     value=np.zeros((self.ny, self.nx))
    #     value [0,:]=Q[0,:]/self.kx[0,:]
    #     value [-1,:]=Q[-1,:]/self.kx[-1,:]
    #     value [:,-1]=Q[:,-1]/self.ky[:,-1]
        
        
    #     Uout[0,:]=u[0,:]-2*self.dx*value[0,:]
    #     Uout[:,-1]=u[:,-1]-2*self.dy*value[:,-1]
    #     Uout[-1,:]=u[-1,:]-2*self*dx*value[-1,:]
        
        
    #     u[0,1:-1] = u0[0,1:-1]+ Diffy*self.dt*((u0[1,1:-1]-2*u[0,1:-1]+ Uout_X[0,1:-1])/self.dy2) + Diffx*self.dt*((u0[0, 0:-2]-2*u0[0:1:-1]+u0[0,2:])/self.dx2) 
    #     # u[1:-1,0] = u0[0,1:-1]+ Diffy*self.dt*((u0[1,1:-1]-2*u[0,1:-1]+ U0[0,1:-1])/self.dy2) + Diffx*self.dt*((u0[0, 0:-2]-2*u0[0:1:-1]+Uout_Y[0,2:])/self.dx2) 
        

        
    #     return value