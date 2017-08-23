import numpy as np
import itertools as it
from numpy import linalg as LA
from scipy.integrate import simps

#---------------------------------------------------------------------------------
def kinetic_op_1d(m,deltax,n):
    """
    4th order centered difference scheme for the kinetic energy operator in 1D
    f"(x0) = (-f(-2) +16f(-1) -30f(0) +16f(+1) -f(+2))/12Dx^2
    """
    T = np.zeros((n,n),dtype=float)
    
    #--- upper edge
    i=0
    T[i,i]   = -30.0e+0 
    T[i,i+1] = +16.0e+0
    T[i,i+2] = -1.00e+0

    i=1
    T[i,1-i] = +16.0e+0
    T[i,i]   = -30.0e+0 
    T[i,i+1] = +16.0e+0
    T[i,i+2] = -1.00e+0

    #--- matrix bulk
    for i in range(2,(n-2)):
        T[i,i-2] = -1.00e+0
        T[i,i-1] = +16.0e+0
        T[i,i]   = -30.0e+0 
        T[i,i+1] = +16.0e+0
        T[i,i+2] = -1.00e+0
    
    #--- lower edge
    i=n-2
    T[i,i-2] = -1.00e+0
    T[i,i-1] = +16.0e+0
    T[i,i]   = -30.0e+0 
    T[i,i+1] = +16.0e+0

    i=n-1
    T[i,i-2] = -1.00e+0
    T[i,i-1] = +16.0e+0
    T[i,i]   = -30.0e+0 
    
    #---------------------------
    T=T/(-12.0e+0 * 2.0e+0 * m * (deltax**2))

    #----------------------------
    return T

#------------------------------------------------------------

def hamiltonian_1d(T,V):
    """
    Generates 1D Hamiltonian given the kinetic energy operator T and the diagonal potential operator V
    """
    H=np.empty_like(T)
    H[:,:] = T
    for i in range(V.size):
        H[i,i] = H[i,i] + V[i]
    return H

#------------------------------------------------------------


def solve_eigenstates(H):
    """
    performs diagonalization of Hamiltonian
    """
    eigen_val,eigen_vec=np.linalg.eigh(H)
    return eigen_val,eigen_vec


#------------------------------------------------------------

def hamiltonian_diag_1d(mu,V,deltax):
    """
     Diagonalizes a 1D hamiltonian with 4th order centered differenece kinetic energy operator
    
     mu - mass in a.u.
     V - potential in a.u.
     deltax - discretization step to be used in the kinetic energy operator
    """


    #constructs kinetic energy operator
    T=kinetic_op_1d(mu,deltax,V.size)
    #constructs hamiltonian
    H=hamiltonian_1d(T,V)
    #performs diagonalization
    eigen_val,eigen_vec=solve_eigenstates(H)

    return eigen_val,eigen_vec

#---------------------------------------------------------------------------------

def franck_condon(x,u,w,dipole=None):
    """
    computes the franck-condon amplitude between real vectors u and w FC = < u | w >
    """
    if(dipole is None):
        dipole = np.ones(x.shape,dtype=float)

    deltax=x[1]-x[0]
    fc=simps((u * dipole * w),x)/deltax
    return fc

#---------------------------------------------------------------------------------

def all_franck_condon(x,u,w,dipole=None):
    """
    computes the franck condon amplitude matrix between two sets of vectors u and w
    """
    n=u.shape[1]
    m=w.shape[1]

    if(dipole == None):
        dipole = np.ones(x.shape,dtype=float)

    fc=np.zeros((n,m),dtype=float)
    for i in range(n):
        for j in range(m):
            fc[i,j]=franck_condon(x,u[:,i],w[:,j],dipole)
    return fc

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

def gen_indexes(nmodes,fc_nv):
    nt=1
    for i in range(nmodes):
        nt = nt * fc_nv[i]
    all_indx=np.zeros((nt,nmodes),dtype=int)

    l=[]
    for i in range(nmodes):
        l.append([x for x in range(fc_nv[i])])
    indexes=[list(m) for m in list(it.product(*l))]
    if(len(indexes) != nt):
        print('error')
        exit()
    return indexes

def multd_fc(e_M,fc_M,fc_nv,fc_indx):
    fc=1e0
    e=0e0
    for i in range(fc_nv.shape[0]):
        fc = fc * fc_M[i,fc_indx[i]] 
        e = e + e_M[i,fc_indx[i]] 
    return fc,e

def do_all_fc(e_M,fc_M,fc_nv,all_indx):
    nt=len(all_indx)
    FC=np.zeros(nt,dtype=float)
    E=np.zeros(nt,dtype=float)
    for i in range(nt):
         FC[i],E[i]=multd_fc(e_M,fc_M,fc_nv,all_indx[i])
    return FC,E
