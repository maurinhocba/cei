
"""
Multiple Degree Of Freedom
Dynamics for linear systems
"""

import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt


# FUNCTIONS
def modes(Kc,Mc,ntype='max'):
    """
    Determine normal modes and retrive them ordered and normalized
    Input:  Kc,Mc,ntype='max'
            ntype is normalization type (criteria) - default 'max'
            means that the larger element in each mode is 1
            could also be: mass, stiff, norm
    Output: PHI, om
    """
    
    # find
    om, PHI = la.eig(Kc,Mc)
    # reorder
    om=np.sqrt(om.real, dtype=float)
    order=om.argsort()
    om=om[order]
    PHI=PHI[:, order]
    # normalize
    if ntype=='max':
        maxabs=np.argmax( np.abs(PHI), axis=1)
        for col in range(maxabs.shape[0]):
            PHI[:,col]=PHI[:,col] / PHI[maxabs[col],col]
            
    else:
        raise Exception('Normalization criteria {} not found.'.format(ntype))
    
    
    return PHI, om
    
    
def uncoup(PHI,Kc,Mc,Pc=None,u0=None,v0=None):
    """
    Determine constant matrix coeficients of the IVP (including
    the eq. of motion and the initial conditions) of a MDoF
    system in terms of the modal coordinates (q).
    Input:  PHI,Kc,Mc,Pc,u0,v0
    Output:
        eqmot: dict conaining Kb,Mb,Pb
        ic: dict containing q0,dq0
    """
    
    # eq. of motion
    Kb=np.diagonal( np.matmul( np.transpose(PHI), np.matmul(Kc,PHI) ) )
    Mb=np.diagonal( np.matmul( np.transpose(PHI), np.matmul(Mc,PHI) ) )
    if Pc is not None:
        Pb=         np.matmul( np.transpose(PHI),           Pc      )
    else:
        Pb=np.zeros(Mb.shape)
    eqmot={'Kb':Kb, 'Mb':Mb, 'Pb':Pb}
    # initial conditions
    q0=np.zeros(Pb.shape)
    dq0=np.zeros(Pb.shape)
    ic={'q0':q0, 'dq0':dq0}
    if u0 is not None:
        for i in range(u0.shape[0]):
            q0[i]=np.matmul( np.transpose(PHI[:,i]), np.matmul(Mc,u0) ) / Mb[i]
        
        ic['q0']=q0
        
    if v0 is not None:
        for i in range(u0.shape[0]):
            dq0[i]=np.matmul( np.transpose(PHI[:,i]), np.matmul(Mc,v0) ) / Mb[i]
        
        ic['dq0']=dq0
    
    
    return eqmot, ic


# RUN
if __name__ == '__main__':
    
    if True: # probar 'modes' y 'uncoup'
        Kc = np.array( [ [  285.71,  532.93, ],
                         [  532.93, 6231.86, ] ] )
        Mc=np.diag((1,1))
        PHI, om = modes(Kc,Mc)
        
        Pc=np.array([1000,0], dtype=float)
        u0=np.array([1,3], dtype=float)
        eqmot, ic = uncoup(PHI,Kc,Mc,Pc,u0)
        
        