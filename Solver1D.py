import numpy as np
import sys 
from scipy import sparse
from scipy.sparse.linalg import spsolve




#-----------------------------------------------------------------

# Reaction-advection-diffusion matrix solver for 1D geometry
# no flux condition on both boundaries 

#------------------------------------------------------------------


def rad_solver_1D(c0, n, w, vw, param_dict): 


    """
    # 1D reaction-advection-diffusion solver for solute c

    Variables:
    c0 - solute distribution at previous timestep
    w - water (solvent) distribution
    n - cell volume fraction that upregulates/uptakes solute
    vw - water (solvent) velocity

    Parameters:
    D - diffusion coefficient of solute
    alpha - production rate of solute by cells
    kappa - solution uptake rate by cells (Michaelis-Menten kinetics)
    K - concentration of solute at which uptake is half-maximal 
    delta - degradation rate of solute in solvent
    ny - number of grid points
    dy - grid spacing
    dy2 - dy*dy square of grid spacing 
    dt - timestep 

    """
    L, T, ny, dy, dy2, phi, m, dt, cint, Dc, alpha, kappa, K, delta = param_dict.values()

    #-- Matrix A (to construct Ac = B matrix equation)

    A = np.zeros((ny, ny))

    #no flux bottom boundary 

    A[0, 0] = (w[0] / dt) + (w[0] / (2 * dy)) * (- 3 * vw[0] + 4 * vw[1] - vw[2]) \
        + (vw[0] / (2 * dy)) * (- 3 * w[0] + 4 * w[1] - w[2]) + ((2 * Dc * w[0]) / dy2) + delta
        
    A[0, 1] = - ((2 * Dc * w[0]) / dy2) 

    for i in range(1, ny-1):
        A[i, i] = (w[i] / dt) + (w[i] / (2 * dy)) * (vw[i+1] - vw[i-1]) \
            + (vw[i] / (2 * dy)) * (w[i+1] - w[i-1]) + ((2 * Dc * w[i]) / (dy2)) + delta

    for i in range(1, ny-1):
        A[i, i-1] = - ((w[i] * vw[i]) / (2 * dy)) + (Dc / (4 * dy2)) * (w[i+1] - w[i-1]) - ((Dc * w[i]) / dy2)

    for i in range(1, ny-1):
        A[i, i+1] =  ((w[i] * vw[i]) / (2 * dy)) - (Dc / (4 * dy2)) * (w[i+1] - w[i-1]) - ((Dc * w[i]) / dy2)
                                      
    #no flux on upper boundary

    A[ny - 1, ny - 1] = (w[ny-1]/ dt) + (w[ny-1] / (2 * dy)) * (3 * vw[ny-1] - 4 * vw[ny-2] + vw[ny-3]) \
         + (vw[ny-1] / (2*dy)) * (3 * w[ny-1] - 4 * w[ny-2] + w[ny-3]) + ((2 * Dc * w[ny-1])/ dy2) + delta

    A[ny - 1, ny - 2] = - ((2 * Dc * w[ny - 1]) / dy2) 

    A = sparse.csr_matrix(A)

    #-- source terms B

    B = np.zeros((ny))

    for i in range(ny):     
        B[i] =  ((w[i] * c0[i]) / dt) 
        # additional (optional) reaction terms (production and uptake)
        B[i] += (alpha * n[i] * w[i]) - ((kappa * n[i] * c0[i] * w[i]) / (K + c0[i]))             

    #-- solve matrix equation for solute c

    c = spsolve(A, B)
    
    return c


