import numpy as np
import sys 
from scipy import sparse
from scipy.sparse.linalg import spsolve


#-----------------------------------------------------------------

# Reaction-advection-diffusion matrix solver for 2D axisymmetric geomtery
# (i, j) corresponds to (z, r) for global variables 
# r = 0 is line of symmetry, r = L no flux condition 
# z = 0, L no flux conditions 

#------------------------------------------------------------------


#--- INNER MATRICES

#-- diagonal matrices
def B_matrix(j, w, uw, vw, D, delta, ny, dy, dy2, dt):
    #valid for j = 1, ..., ny-2 

    B = sparse.lil_matrix((ny, ny))

    B[0, 0] = (4 * w[0, j]) + delta + (dy2  * w[0, j] / (D * dt)) \
        + (dy * w[0, j] / (2 * D)) * (vw[0, j+1] - vw[0, j-1]) + (dy * vw[0, j] / (2 * D)) * (w[0, j+1] - w[0, j-1])  + ((dy * vw[0, j] * w[0, j]) / (j * D)) \
        + (dy * w[0, j] / (2 * D)) * (-3*uw[0, j] + 4*uw[1, j] - uw[2, j]) + (dy * uw[0, j] / (2 * D)) * (-3 * w[0, j] + 4 * w[1, j] - w[2, j])
    B[0, 1] = - 2 * w[0, j]

    for i in range(1, ny-1):
        B[i, i] = (4 * w[i, j]) + delta + (dy2  * w[i, j] / (D * dt)) \
        + (dy * w[i, j] / (2 * D)) * (vw[i, j+1] - vw[i, j-1]) + (dy * vw[i, j] / (2 * D)) * (w[i, j+1] - w[i, j-1])  + ((dy * vw[i, j] * w[i, j]) / (j * D)) \
        + (dy * w[i, j]/ (2 * D)) * (uw[i+1, j] - uw[i-1, j]) + (dy * uw[i, j] / (2 * D)) * (w[i+1, j] - w[i-1, j])
        B[i, i-1] = - ((dy * w[i, j] * uw[i, j]) / (2 * D)) - w[i, j] + (w[i, j] / (2*j)) + (w[i+1, j] - w[i-1, j])/4
        B[i, i+1] = ((dy * w[i, j] * uw[i, j]) / (2 * D)) - w[i, j] - (w[i, j] / (2*j)) - (w[i+1, j] - w[i-1, j])/4

    B[ny-1, ny-1] = (4 * w[ny-1, j]) + delta + (dy2  * w[ny-1, j] / (D * dt)) \
        + (dy * w[ny-1, j] / (2 * D)) * (vw[ny-1, j+1] - vw[ny-1, j-1]) + (dy * vw[ny-1, j] / (2 * D)) * (w[ny-1, j+1] - w[ny-1, j-1])  + ((dy * vw[ny-1, j] * w[ny-1, j]) / (j * D)) \
        + (dy * w[ny-1, j]/ (2 * D)) * (3*uw[ny-1, j] - 4*uw[ny-2, j] + uw[ny-3, j]) + (dy * uw[ny-1, j] / (2 * D)) * (3 * w[ny-1, j] - 4 * w[ny-2, j] + w[ny-3, j])
    B[ny-1, ny-2] = - 2 * w[ny-1, j]

    
    return B

#-- off-diagonal matrices 
def D1_matrix(j, w, vw, D, ny, dy):
    # off-diagonal matrix for A valid for j = 1, ... n-2

    D1 = sparse.lil_matrix((ny, ny))

    for i in range(0, ny):
        D1[i, i] = - ((dy * w[i, j] * vw[i, j]) / (2 * D)) - w[i, j] + (w[i, j] / (2*j)) + (w[i, j+1] - w[i, j-1])/4

    return D1

def D2_matrix(j, w, vw, D, ny, dy):
    # off-diagonal matrix for A valid for j = 1, ... n-2

    D2 = sparse.lil_matrix((ny, ny))

    for i in range(0, ny):
        D2[i, i] = ((dy * w[i, j] * vw[i, j]) / (2 * D)) - w[i, j] - (w[i, j] / (2*j)) - (w[i, j+1] - w[i, j-1])/4

        
    return D2

#----- MAIN MATRIX SOLVER 

def rad_solver_2Daxi(c0, w, n, uw, vw, param_dict):

    """
    # 2D reaction-advection diffusion solver for solute 

    Variables:
    c0 - solute distribution at previous timestep
    w - current water (solvent) distribution
    w0 - water (solvent) distribution at previous timestep
    n - cell volume fraction that upregulates/uptakes solute
    uw, vw - z and r components of water (solvent) velocity

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

    L, T, nx, ny, dy, dy2, phi, m, dt, cint, Pe, alpha, kappa, K, delta = param_dict.values()
    D = 1 / Pe #FIX #TODO 

    #-- First diagonal matrix C (arises from line of symmetry boundary, valid for j=0)
    # alternative form of B_matrix at boundary 
    C = sparse.lil_matrix((ny, ny))

    C[0, 0] = (6 * w[0, 0]) + delta + (dy2  * w[0, 0] / (D * dt)) \
        + (dy * w[0, 0] / D) * (-3*vw[0, 0] + 4 * vw[0, 1] - vw[0, 2])  \
        + (dy * w[0, 0]/ (2 * D)) * (-3*uw[0, 0] + 4*uw[1, 0] - uw[2, 0]) + (dy * uw[0, 0] / (2 * D)) * (-3 * w[0, 0] + 4 * w[1, 0] - w[2, 0])
    
    C[0, 1] = - 2 * w[0, 0] 

    for i in range(1, ny-1):
        C[i, i] = (6 * w[i, 0]) + delta + (dy2  * w[i, 0] / (D * dt)) \
        + (dy * w[i, 0] / D) * (-3*vw[i, 0] + 4 * vw[i, 1] - vw[i, 2])  \
        + (dy * w[i, 0]/ (2 * D)) * (uw[i+1, 0] - uw[i-1, 0]) + (dy * uw[i, 0] / (2 * D)) * (w[i+1, 0] - w[i-1, 0])

        C[i, i-1] =  - ((dy * w[i, 0] * uw[i, 0]) / (2 * D)) - w[i, 0] + (w[i+1, 0] - w[i-1, 0])/4

        C[i, i+1] = ((dy * w[i, 0] * uw[i, 0]) / (2 * D)) - w[i, 0] - (w[i+1, 0] - w[i-1, 0])/4
        
    C[ny-1, ny-2] = - 2 * w[ny-1, 0]

    C[ny-1, ny-1] =  (6 * w[ny-1, 0]) + delta + (dy2  * w[ny-1, 0] / (D * dt)) \
        + (dy * w[ny-1, 0] / D) * (-3*vw[ny-1, 0] + 4 * vw[ny-1, 1] - vw[ny-1, 2]) \
        + (dy * w[ny-1, 0]/ (2 * D)) * (3*uw[ny-1, 0] - 4*uw[ny-2, 0] + uw[ny-3, 0]) + (dy * uw[ny-1, 0] / (2 * D)) * (3 * w[ny-1, 0] - 4 * w[ny-2, 0] + w[ny-3, 0])

    #-- Final matrix B (boundary at j = ny-1 (r = L))
    B = sparse.lil_matrix((ny, ny))

    B[0, 0] = (4 * w[0, ny-1]) + delta + (dy2  * w[0, ny-1] / (D * dt)) \
        + (dy * w[0, ny-1] / (2 * D)) * (3*vw[0, ny-1] - 4*vw[0, ny-2] + vw[0, ny-3]) + (dy * vw[0, ny-1] / (2 * D)) * (3*w[0, ny-1] - 4*w[0, ny-2] + w[0, ny-3]) + ((dy * vw[0, ny-1] * w[0, ny-1]) / ((ny-1) * D)) \
        + (dy * w[0, ny-1]/ (2 * D)) * (-3*uw[0, ny-1] + 4*uw[1, ny-1] - uw[2, ny-1]) + (dy * uw[0, ny-1] / (2 * D)) * (-3 * w[0, ny-1] + 4 * w[1, ny-1] - w[2, ny-1])
    B[0, 1] = - 2 * w[0, ny-1]

    for i in range(1, ny-1):
        B[i, i] = (4 * w[i, ny-1]) + delta + (dy2  * w[i, ny-1] / (D * dt)) \
        + (dy * w[i, ny-1] / (2 * D)) * (3*vw[i, ny-1] - 4*vw[i, ny-2] + vw[0, ny-3]) + (dy * vw[i, ny-1] / (2 * D)) * (3*w[i, ny-1] - 4*w[i, ny-2] + w[0, ny-3]) + ((dy * vw[i, ny-1] * w[i, ny-1]) / ((ny-1) * D)) \
        + (dy * w[i, ny-1]/ (2 * D)) * (uw[i+1, ny-1] - uw[i-1, ny-1]) + (dy * uw[i, ny-1] / (2 * D)) * (w[i+1, ny-1] - w[i-1, ny-1])
        B[i, i-1] = - ((dy * w[i, ny-1] * uw[i, ny-1]) / (2 * D)) - w[i, ny-1] + (w[i, ny-1] / (2*(ny-1))) + (w[i+1, ny-1] - w[i-1, ny-1])/4
        B[i, i+1] = ((dy * w[i, ny-1] * uw[i, ny-1]) / (2 * D)) - w[i, ny-1] - (w[i, ny-1] / (2*(ny-1))) - (w[i+1, ny-1] - w[i-1, ny-1])/4

    B[ny-1, ny-1] = (4 * w[ny-1, ny-1]) + delta + (dy2  * w[ny-1, ny-1] / (D * dt)) \
        + (dy * w[ny-1, ny-1] / (2 * D)) * (3*vw[ny-1, ny-1] - 4*vw[ny-1, ny-2] + vw[ny-1, ny-3]) + (dy * vw[ny-1, ny-1] / (2 * D)) * (3*w[ny-1, ny-1] - 4*w[ny-1, ny-2] + w[ny-1, ny-3]) + ((dy * vw[ny-1, ny-1] * w[ny-1, ny-1]) / ((ny-1) * D)) \
        + (dy * w[ny-1, ny-1]/ (2 * D)) * (3*uw[ny-1, ny-1] - 4*uw[ny-2, ny-1] + uw[ny-3, ny-1]) + (dy * uw[ny-1, ny-1] / (2 * D)) * (3 * w[ny-1, ny-1] - 4 * w[ny-2, ny-1] + w[ny-3, ny-1])
    B[ny-1, ny-2] = - 2 * w[ny-1, ny-1]

    #-- First off-diagonal E1 (for boundary j=0 (r=0))
    E1 = sparse.lil_matrix((ny, ny))

    for i in range(1, ny-1):
        E1[i, i] = - 4 * w[i, 0]

    #-- Last off-diagonal E2 (for boundary j=ny-1 (R=L))
    E2 = sparse.lil_matrix((ny, ny))

    for i in range(1, ny-1):
        E2[i, i] = - 2 * w[i, ny-1]


    #-- Main matrix A
    ny2 = ny * ny
    A = sparse.lil_matrix((ny2, ny2))

    A[0:ny, 0:ny] = C
    A[0:ny, ny:2*ny] = E1

    for i in range(1, ny-1):
        for j in range(0, ny):
            if i==j:
                A[i*ny:(i+1)*ny, j*ny:(j+1)*ny] = B_matrix(i, w, uw, vw, D, delta, ny, dy, dy2, dt)
            if j == i+1:
                A[i*ny:(i+1)*ny, j*ny:(j+1)*ny] = D2_matrix(i, w, vw, D, ny, dy)
            if j == i-1:
                A[i*ny:(i+1)*ny, j*ny:(j+1)*ny] = D1_matrix(i, w, vw, D, ny, dy)


    A[ny2-ny:ny2, ny2-2*ny:ny2-ny] = E2
    A[ny2-ny:ny2, ny2-ny:ny2] = B

    A = sparse.csr_matrix(A)

    #-- source terms F
    F = (c0 * w * dy2) / (D * dt)
    # additional (optional) reaction term (production and uptake)
    F += (alpha * n * w)  - ((kappa * n * c0 * w) / (K + c0))  

    # convert to 1D vector for use in matrix equation 
    F1d = np.zeros((ny2))
    for i in range(ny):
        for j in range(ny):
            F1d[j*ny + i] = F[i, j]

    #-- solve matrix equation for solute c 

    c = np.zeros((ny, ny))

    c = spsolve(A, F1d)

    return n 

