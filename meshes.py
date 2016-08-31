import numpy as np
def f_quad(x,Lx, Ly):
    #-------------------------
    # quadratic mesh function
    #-------------------------
    return 2*(x-0.5*Lx)**2/Lx/np.sqrt(3) + 0.5*Ly-0.5*0.25*Lx/np.sqrt(3)
def f_V(x,Lx, Ly):
    #-------------------------
    # V shape mesh function
    #-------------------------
    return np.where(x<=0.5*Lx,
                    (-1/np.sqrt(3))*x+Lx/4./np.sqrt(3)+0.5*Lx,
                    np.where(x>0.5*Lx,
                    1/np.sqrt(3)*(x-0.5*Lx)-Lx/4./np.sqrt(3)+0.5*Lx,0))

def compt_to_phys_SB(Y, Ly, fx):
    #-------------------------------------
    # non-orthogonal computational domain 
    # for solid body rotation test case
    #-------------------------------------
    return np.where(Y>0.5*Ly, fx+(Y-0.5*Ly)*(Ly-fx)/(0.5*Ly),Y*fx/(0.5*Ly))

def phys_to_compt_SB(Y, Ly, fx):
    #-------------------------------------
    # non-orthogonal computational domain 
    # for solid body rotation test case
    #-------------------------------------
    return np.where(Y>fx, (Y-fx)*0.5*Ly/(Ly-fx) + (0.5*Ly),Y*0.5*Ly/fx)