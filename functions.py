import numpy as np
import itertools as it
from numpy import linalg as LA
from scipy.integrate import simps
import matplotlib.pyplot as plt

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

#-------------------------------------------------------------------------------
def compute_1d_fc(mu,x,V_g,V_c,V_f,nc,nf,dipole=None):
    """
    computes the 1d FC factors between the gs-ce and ce-gs

    mu       : mass in a.u.
    x        : coordinate space vector 
    V_g      : PES for state |g>
    V_c      : PES for state |c>
    V_f      : PES for state |f>
    omega_gc : energy difference between the minimum of the |g> and |c> PES (in a.u.)
    omega_gf : energy difference between the minimum of the |g> and |f> PES (in a.u.)
    Gamma    : lifetime broadening of state |c> (in a.u.)
optional:
    omega    : desired incoming photon energies (in a.u.)
    nc       : desired number of vibrational levels considereg for |c>
    """

    omega_gc=V_c.min() - V_g.min()
    omega_gf=V_f.min() - V_g.min()

    deltax=(x[1] - x[0])
    eg,psi_g=hamiltonian_diag_1d(mu,V_g - V_g.min(),deltax)
    ec,psi_c=hamiltonian_diag_1d(mu,V_c - V_c.min(),deltax)
    ef,psi_f=hamiltonian_diag_1d(mu,V_f - V_f.min(),deltax)
    fc_gc=all_franck_condon(x,psi_g[:,0:3],psi_c[:,0:nc])
    fc_fc=all_franck_condon(x,psi_c[:,0:nc],psi_f[:,0:nf],dipole)

    return fc_gc[0,:],fc_fc[:,:],eg[0],ec[0:nc],ef[0:nf]

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

#----------------------------------------------------
def multd_0vc_fc(e_M,fc_M,nmodes,fc_indx_vc):
    fc=1e0
    e=0e0
    #print(fc_indx_vc)
    for i in range(nmodes):
        fc = fc * fc_M[i,fc_indx_vc[i]] 
        #print(fc_M[i,fc_indx_vc[i]])
        e = e + e_M[i,fc_indx_vc[i]] 
    return fc,e

#----------------------------------------------------
def multd_vcvf_fc(e_M,fc_M,nmodes,fc_indx_vc,fc_indx_vf):
    fc=1e0
    e=0e0
    #print(fc_indx_vc,fc_indx_vf)
    for i in range(nmodes):
#        for j in range(nmodes):
        fc = fc * fc_M[i,fc_indx_vc[i],fc_indx_vf[i]] 
        #print(fc_M[i,fc_indx_vc[i],fc_indx_vf[j]],fc_indx_vc[i],fc_indx_vf[j])
        e = e + e_M[i,fc_indx_vf[i]]
    return fc,e

#----------------------------------------------------
def read_all_pot(nmodes, fc_init_pot,fc_decay_pot,fc_fin_pot):
    R={};Vg={};Vd={};Vf={}
    for i in range(nmodes):
        R[i],Vg[i]=np.genfromtxt(fc_init_pot[i],skip_header=1,unpack=True,dtype=float)
        R[i],Vd[i]=np.genfromtxt(fc_decay_pot[i],skip_header=1,unpack=True,dtype=float)
        R[i],Vf[i]=np.genfromtxt(fc_fin_pot[i],skip_header=1,unpack=True,dtype=float)
    return R,Vg,Vd,Vf

#----------------------------------------------------
def plot_all_pot(nmodes,R,Vg,Vd,Vf):
    f,ax=plt.subplots(3,nmodes)
    for i in range(nmodes):
        ax[2,i].plot(R[i],Vg[i])
        ax[1,i].plot(R[i],Vd[i])
        ax[0,i].plot(R[i],Vf[i])
    plt.show()
    return

#------------------------------------------------------
def print_1d_fc(nmodes,fc_nvc,fc_nvf,fc_0vc,fc_vcvf,thsh):
    for i in range(nmodes):
        print('mode',i,',',fc_nvc[i],fc_nvf[i]);
        print('FC ampl., gs -> ce');
        for j in range(fc_nvc[i]):
            print('0 ->'+str(j),',',fc_0vc[i,j])

        print('FC ampl., ce -> f')
        for j in range(fc_nvc[i]):
            for k in range(fc_nvf[i]):
                print(str(j)+'->'+str(k),',',fc_vcvf[i,j,k])
        print()
    return

#------------------------------------------------------
# This prints the FC amplitudes in Spec format
#
def print_multd_fc(nmodes,fc_nvc,fc_nvf,fc_0vc,fc_vcvf,e_g,e_c,e_f,all_indx_vc,all_indx_vf,thsh):
    f_out_0vc=open('fc_0vc_temp','w')
    f_out_vcvf=open('fc_vcvf_temp','w')
    f_evc=open('evc_temp','w')
    f_evf=open('evf_temp','w')
    #print('dirty compatibility fix file',file=f_out_0vc)
    #print('dirty compatibility fix file',file=f_out_vcvf)


    for f_ev in [f_evc,f_evf]:
        print("Calculating eigenvalue(s) and eigenvector(s)\n",
              "for the initial state...\n",
              "Using complete matrix diagonalization procedure,\n",
              "\n","Eigenvalues:\n",
              '""""""""""""\n','---------------------------------',file=f_ev)
        e_0=0.0e0
        for e_i in e_g:
            e_0=e_0+e_i
        print('|     {}        |     {: 8.6E}    |'.format(0,e_0),file=f_ev)
        print('------------------------------------',file=f_ev)
        print('End of file, reading finish!',file=f_ev)

    for f_ev in [f_evc,f_evf]:
        print("Calculating eigenvalue(s) and eigenvector(s)\n",
              "from final state...\n",
              "Using complete matrix diagonalization procedure,\n",
              "\n","Eigenvalues:\n",
              '""""""""""""\n','---------------------------------',file=f_ev)
    for f_fc in [f_out_0vc,f_out_vcvf]:
        print("Spectrum:",file=f_fc)
        print(" =========",file=f_fc)
        print("     *E/a.u.*              *AMPT*              *FC*                *I->F*",file=f_fc)


    print('--------------------------------------')
    print(' 0  -> vc,       evc,          fc_amp')
    vc=0;vc_negl=0;
    for index in all_indx_vc:
        fc_t,e_t=multd_0vc_fc(e_c,fc_0vc,nmodes,index)
        if ( np.fabs(fc_t) >= thsh ):
            label=''.join(str(0) for a in index)+' -> '+''.join(str(a) for a in index)
            print('{}    {: 8.6E}      {: 8.6E}'.format(label,e_t,fc_t))
            print('    {: 8.6E}      {: 8.6E}      {: 8.6E}           {}'.format(e_t,fc_t,fc_t*fc_t,label),file=f_out_0vc)
            print('|     {}        |     {: 8.6E}    |'.format(vc,e_t),file=f_evc)
            print('------------------------------------',file=f_evc)
            vc=vc+1
        else:
            vc_negl=vc_negl+1
    # total number of states included
    vc_tot=vc

    print('--------------------------------------')
    print(' vc  -> vf,      evf,          fc_amp')
    vc=0;vf=0;
    for index_vc in all_indx_vc:
        fc_check,e_check=multd_0vc_fc(e_c,fc_0vc,nmodes,index_vc)
        if ( np.fabs(fc_check) >= thsh ):
            for index_vf in all_indx_vf:
                fc_t,e_t=multd_vcvf_fc(e_f,fc_vcvf,nmodes,index_vc,index_vf)
                label=''.join(str(b) for b in index_vc)+' -> '+''.join(str(a) for a in index_vf)
                print('{}    {: 8.6E}      {: 8.6E}'.format(label,e_t,fc_t))
                print('    {: 8.6E}      {: 8.6E}      {: 8.6E}         {}'.format(e_t,fc_t,fc_t*fc_t,label),file=f_out_vcvf)
                if(vc==0):
                    print('|     {}        |  {: 8.6E}      |'.format(vf,e_t),file=f_evf)
                    print('------------------------------------',file=f_evf)
                    vf=vf+1
            vc=vc+1

    print(vc_negl," states were found with fc amp below the threshold",thsh)
    print(vc_tot," states remain")
    print('\n The eSPec program finished successfully!',file=f_out_0vc)
    print('\n The eSPec program finished successfully!',file=f_out_vcvf)
    print("\n\nneglected: ",vc_negl,"\n",file=f_out_0vc)
    return vc_negl

#------------------------------------------------------
def print_1d_eigval(nmodes,fc_nvc,fc_nvf,e_g,e_c,e_f):
    for i in range(nmodes):
        print('mode',i,',',fc_nvc[i],fc_nvf[i]);
        print('eigenvalues ');
        print('e0 =',e_g[i],'\n')
        print('evc')
        for j in range(fc_nvc[i]):
            print(str(j),' | ',e_c[i,j])
        print()
        print('evf')
        for k in range(fc_nvf[i]):
            print(str(k),' | ',e_f[i,k])
        print()
    return

#----------------------------------------------------
def get_multd_fc(nmodes,inp_file,thsh=1e-8):
    mass=np.zeros(nmodes,dtype=float)
    fc_nvc=np.zeros(nmodes,dtype=int)
    fc_nvf=np.zeros(nmodes,dtype=int)
    fc_init_pot=[]
    fc_decay_pot=[]
    fc_fin_pot=[]

    f_inp=open(inp_file,'r')
    for line in f_inp:
        if 'fc_mass' in line:
            for i in range(nmodes):
                mass[i]=(1.836152667539e+3)*np.float(line.split()[i+1])
        elif 'fc_nvc' in line:
            for i in range(nmodes):
                fc_nvc[i]=np.float(line.split()[i+1])
        elif 'fc_nvf' in line:
            for i in range(nmodes):
                fc_nvf[i]=np.float(line.split()[i+1])
        elif 'fc_init_pot' in line:
            for i in range(nmodes):
                fc_init_pot.append(line.split()[i+1])
        elif 'fc_decay_pot' in line:
            for i in range(nmodes):
                fc_decay_pot.append(line.split()[i+1])
        elif 'fc_fin_pot' in line:
            for i in range(nmodes):
                fc_fin_pot.append(line.split()[i+1])
    # print(mass);print(fc_nvc);print(fc_nvf);print(fc_init_pot);print(fc_decay_pot);print(fc_fin_pot)

    # Here we read the potentials from the files
    R,V_g,V_c,V_f=read_all_pot(nmodes, fc_init_pot,fc_decay_pot,fc_fin_pot)
    #plot_all_pot(nmodes,R,V_g,V_c,V_f)

    # Here we generate the indexes for all
    #possible combinations of 1D vibrational levels
    all_indx_vc=gen_indexes(nmodes,fc_nvc)
    all_indx_vf=gen_indexes(nmodes,fc_nvf)
    
    # Here compute 1D FC for all modes--------------------------------------------------------------------------
    fc_0vc=np.zeros((nmodes,fc_nvc.max()),dtype=float)
    fc_vcvf=np.zeros((nmodes,fc_nvc.max(),fc_nvf.max()),dtype=float)
    e_g=np.zeros(nmodes,dtype=float)
    e_c=np.zeros((nmodes,fc_nvc.max()),dtype=float)
    e_f=np.zeros((nmodes,fc_nvf.max()),dtype=float)
    for i in range(nmodes):
        fc_0vc[i,0:fc_nvc[i]],fc_vcvf[i,0:fc_nvc[i],0:fc_nvf[i]],e_g[i],e_c[i,0:fc_nvc[i]],e_f[i,0:fc_nvf[i]]=\
        compute_1d_fc(mass[i],R[i],V_g[i],V_c[i],V_f[i],fc_nvc[i],fc_nvf[i])
    #------------------------------------------------------------------------------------------------------------
    
    # print computed 1d fc amplitudes
    print_1d_fc(nmodes,fc_nvc,fc_nvf,fc_0vc,fc_vcvf,thsh)
    print_1d_eigval(nmodes,fc_nvc,fc_nvf,e_g,e_c,e_f)

    # This functions computes all necessary multimode amplitudes FC amplitudes
    # and prints them into eSPec-like files for later use in the eSPec-Raman script
    vc_negl=print_multd_fc(nmodes,fc_nvc,fc_nvf,fc_0vc,fc_vcvf,e_g,e_c,e_f,all_indx_vc,all_indx_vf,thsh)

    #do_all_multd_fc(e_M,fc_M,fc_nv,all_indx)
    return
#-----------------------------------------------------------------------------------------------------------------
