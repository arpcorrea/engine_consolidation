import numpy as np

class HeatTransfer:

    def __init__(self, deck,meshes,BC):
        self.dt = float(deck.doc["Simulation"]["Time Step"])
        self.dx2 = meshes.Mdx*meshes.Mdx
        self.dy2 = meshes.Mdy*meshes.Mdy
        self.BC=BC
        self.Mdy=meshes.Mdy
        self.Mdx=meshes.Mdx
        self.h= float(deck.doc["Boundary Condition"]["Convective Coefficient"])
# -------------- BEGIN HEAT TRANSFER CALCULATION---------- 

    def do_convection(self, u0, u, Diffx, Diffy,Q):
        # Propagate with forward-difference in time, central-difference in space        

        Q=np.zeros(np.shape(self.BC.T))
        Q[0,:] =self.h*(float(self.BC.Troom) - u[0,:])
        Q[-1,:]=self.h*(float(self.BC.Troom) - u[-1,:])
        Q[:,-1]=self.h*(float(self.BC.Troom) - u[:,-1])     
        Q[:,0] =self.BC.q  
        Q[0,0]=self.h*(float(self.BC.Troom)  - u[0,0])
        Q[-1,0]=self.h*(float(self.BC.Troom) - u[-1,0])
        Q[0,-1]=self.h*(float(self.BC.Troom) - u[0,-1])
        Q[-1,-1]=self.h*(float(self.BC.Troom) - u[-1,-1])
        self.Q=Q
        
        
        
        
        value=np.zeros(np.shape(self.BC.T))
        value[0,:] = -Q[0,:]/ self.BC.Ky[0,:]
        value[-1,:]= -Q[-1,:]/self.BC.Ky[-1,:]
        value[:,-1]= -Q[:,-1]/self.BC.Kx[:,-1]
        value[:,0 ]= -Q[:,0]/ self.BC.Kx[:,0]
        value[0,0]   =  -Q[0,0]/self.BC.Kx[0,0]     - Q[0,0]/self.BC.Ky[0,0]
        value[-1,0]  =  -Q[-1,0]/self.BC.Kx[-1,0]   - Q[-1,0]/self.BC.Ky[-1,0]
        value[0,-1]  =  -Q[0,-1]/self.BC.Kx[0,-1]   - Q[0,-1]/self.BC.Ky[0,-1]
        value[-1,-1] =  -Q[-1,-1]/self.BC.Kx[-1,-1] - Q[-1,-1]/self.BC.Ky[-1,-1]       
        self.value=value
        
        
        Uout=np.zeros(np.shape(self.BC.T))
        Uout[0,:] =u[0,:] -2*self.Mdy[0,:] *value[0,:]        
        Uout[-1,:]=u[-1,:]-2*self.Mdy[-1,:]*value[-1,:]        
        Uout[:,-1]=u[:,-1]-2*self.Mdx[:,-1]*value[:,-1]
        Uout[:,0]= u[:,0] -2*self.Mdx[:,0] *value[:,0]
        Uout[0,0] = u[0,0]- self.Mdx[0,0]*value[0,0] - self.Mdy[0,0]*value[0,0]
        Uout[0,-1] = u[0,-1]- self.Mdx[0,-1]*value[0,-1] - self.Mdy[0,-1]*value[0,-1]
        Uout[-1,0] = u[-1,0]- self.Mdx[-1,0]*value[-1,0] - self.Mdy[-1,0]*value[-1,0]
        Uout[-1,-1] = u[-1,-1]- self.Mdx[-1,-1]*value[-1,-1] - self.Mdy[-1,-1]*value[-1,-1]    
        self.Uout=Uout
        
    
        
        
        u[0,1:-1]     = u0[0,1:-1]     + Diffy[0,1:-1] *self.dt*      ((u0[1,1:-1]    -2*u[0,1:-1]      + Uout[0,1:-1]) /self.dy2[0,1:-1])      + Diffx[0,1:-1] *self.dt*       ((u0[0, 2:]    -2*u0[0,1:-1]     +u0[0,:-2])/self.dx2[0,1:-1]) 
        u[-1,1:-1]    = u0[-1,1:-1]    + Diffy[-1,1:-1]*self.dt*      ((u0[-2,1:-1]   -2*u[-1,1:-1]     + Uout[-1,1:-1])/self.dy2[-1,1:-1])     + Diffx[-1,1:-1]*self.dt*       ((u0[-1, 2:]   -2*u0[-1,1:-1]    +u0[-1,:-2])/self.dx2[-1,1:-1]) 
        # u[1:-1, 1:-1] = u0[1:-1, 1:-1] + Diffy[1:-1, 1:-1]* self.dt * ((u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/self.dy2[1:-1, 1:-1] ) + Diffx[1:-1, 1:-1]* self.dt *  ((u0[1:-1, 2:] -2*u0[1:-1, 1:-1] +u0[1:-1, :-2])/self.dx2[1:-1, 1:-1] )
        
        u[1:-1,0]     = u0[1:-1,0]     + Diffy[1:-1,0] *self.dt*      ((u0[2:,0 ]     -2*u[1:-1,0]      + u0[0:-2,0]) /self.dy2[1:-1,0])       + Diffx[1:-1, 0]*self.dt*       ((u0[1:-1, 1 ] -2*u0[1:-1,0]     +Uout[1:-1,0])/self.dx2[1:-1,0])        
        u[1:-1,-1]    = u0[1:-1,-1]    + Diffy[1:-1,-1]*self.dt*      ((u0[2:,-1]     -2*u[1:-1,-1]     + u0[0:-2,-1])/self.dy2[1:-1,-1])      + Diffx[1:-1,-1]*self.dt*       ((u0[1:-1, -2 ]-2*u0[1:-1,-1]    +Uout[1:-1,-1])/self.dx2[1:-1,-1]) 
        # u[1:-1, 1:-1] = u0[1:-1, 1:-1] + Diffy[1:-1, 1:-1]* self.dt * ((u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/self.dy2[1:-1, 1:-1] ) + Diffx[1:-1, 1:-1]* self.dt *  ((u0[1:-1, 2:] -2*u0[1:-1, 1:-1] +u0[1:-1, :-2])/self.dx2[1:-1, 1:-1] )

        u[0,0]   = u0[0,0]      +Diffy[0,0]      *self.dt*((u0[1,0]   -2*u[0,0]    + Uout[0,0])/self.dy2[0,0])         +Diffx[0,0]     *self.dt*((u0[0, 1]     -2*u0[0,0]    +Uout[0,0])/self.dx2[0,0])
        u[-1,0]  = u0[-1,0]     +Diffy[-1,0]     *self.dt*((u0[-2,0]  -2*u[-1,0]   + Uout[-1,0])/self.dy2[-1,0])       +Diffx[-1,0]    *self.dt*((u0[-1, 1]    -2*u0[-1,0]   +Uout[-1,0])/self.dx2[-1,0])  
        u[0,-1]  = u0[0,-1]     +Diffy[0,-1]     *self.dt*((u0[1,-1]  -2*u[0,-1]   + Uout[0,-1])/self.dy2[0,-1])       +Diffx[0,-1]    *self.dt*((u0[0, -2]    -2*u0[0,-1]   +Uout[0,-1])/self.dx2[0,-1])
        u[-1,-1] = u0[-1,-1]    +Diffy[-1,-1]    *self.dt*((u0[-2,-1] -2*u[-1,-1]  + Uout[-1,-1])/self.dy2[-1,-1])     +Diffx[-1,-1]   *self.dt*((u0[-1,-2]    -2*u0[-1,-1] +Uout[-1,-1])/self.dx2[-1,-1])

        u0=u.copy()
        
        return u0,u      
        
    def do_timestep(self, u0, u, Diffx, Diffy,Q):
        # Propagate with forward-difference in time, central-difference in space       
        
        u[1:-1, 1:-1] = u0[1:-1, 1:-1] + Diffy[1:-1, 1:-1]* self.dt * ((u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/self.dy2[1:-1, 1:-1] ) + Diffx[1:-1, 1:-1]* self.dt * ( (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/self.dx2[1:-1, 1:-1] ) + self.dt*Q[1:-1,1:-1]/(self.BC.Cp[1:-1,1:-1]*self.BC.Rho[1:-1,1:-1])

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