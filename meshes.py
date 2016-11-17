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
def f_W(x,Lx):
    #-------------------------
    # W shape mesh function
    #-------------------------
    return np.where(x<=0.25*Lx,
                    -(1./np.sqrt(3))*(x)+Lx/8./np.sqrt(3),
                    np.where((x<0.5*Lx)&(x>=0.25*Lx),
                    1./np.sqrt(3)*(x-0.25*Lx)-Lx/8./np.sqrt(3),
                    np.where((x<0.75*Lx)&(x>=0.5*Lx),
                    -(1./np.sqrt(3))*(x-0.5*Lx)+Lx/8./np.sqrt(3),1./np.sqrt(3)*(x-0.75*Lx)-Lx/8./np.sqrt(3))))
def h(x, h0):
    #-------------------------
    # function for basic terrain
    # following coordinate (BTF)
    #-------------------------
    a = 25e3
    lam = 8e3
    return np.where(np.abs(x) <= a,
                    h0*np.cos(np.pi*x/lam)**2 * np.cos(np.pi*x*0.5/a)**2,
                    0)

def compt_to_phys_SB(Y, ymax, ymin, fx):
    #-------------------------------------
    # non-orthogonal computational domain 
    # for solid body rotation test case
    #-------------------------------------
    Ly = ymax + ymin
    # return np.where((Y>=0.5*Ly) & (Y <= ymax), ymax + (Y-ymax)*(fx-ymax)/(0.5*Ly - ymax), \
    #     np.where((Y<0.5*Ly) & (Y >= ymin),  ymin + (Y-ymin)*(fx-ymin)/(0.5*Ly-ymin), Y))#
    return np.where(Y>0.5*(ymax + ymin), fx+(Y-0.5*Ly)*(Ly-fx)/(0.5*Ly),Y*fx/(0.5*Ly))

def phys_to_compt_SB(Y, Ly, fx):
    #-------------------------------------
    # non-orthogonal computational domain 
    # for solid body rotation test case
    #-------------------------------------
    return np.where(Y>fx, (Y-fx)*0.5*Ly/(Ly-fx) + (0.5*Ly),Y*0.5*Ly/fx)


def compt_to_phys_oro(x,Z,zmax, h0):
    #------------------------------
    # computational domain for BTF
    #-----------------------------
    hx = h(x, h0)
    return hx + Z*(zmax - hx)/zmax

def phys_to_compt_oro(x,Z,zmax, h0):
    #------------------------------
    # computational domain for BTF
    #-----------------------------
    hx = h(x, h0)
    return (Z - hx)*zmax/(zmax - hx)

def compt_to_phys_deform(X,Y, Lx, ymin, ymax):
    #-------------------------
    # computational domain for
    # W shape mesh
    #-------------------------
    fx = f_W(X,Lx)
    return np.where(Y>=0, fx+Y*(1-fx/ymax),fx+Y*(1-fx/ymin))

def phys_to_compt_deform(X,Y, Lx, ymin, ymax):
    #-------------------------
    # computational domain for
    # W shape mesh
    #-------------------------
    fx = f_W(X,Lx)
    return np.where(Y>=fx, (Y - fx)*ymax/(ymax - fx),(Y - fx)*ymin/(ymin - fx))