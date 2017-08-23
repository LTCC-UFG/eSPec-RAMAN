import numpy as np
import itertools as it

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
