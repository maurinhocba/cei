
"""
Simple Harmonic Oscillator
"""

import numpy as np
import matplotlib.pyplot as plt


# FUNCTIONS
def ufree_underd(m,c,k,u0,v0,t=0):
    """
    Free vibration response of a simple underdamped harmonic oscillator.
    Input: m,c,k,u0,v0,t=0
    Output:
        u:      dict with array with the evaluation of uf(t) at the
                instants in array t
        cts:    dict of constants determining the response
        msg:    text explaining the output (incluiding the values)
    """
    
    # SYSTEM PROPERTIES
    c_cr=2*np.sqrt(m*k) # critical damping
    xi=c/c_cr # [-] - damping ratio
    if xi >= 1:
        raise Exception('\'ufree_under\' is written for underdamped harmonic oscillators, \
                         thus xi=c/[2*sqrt(m*k)] should be lower than 1. \
                         The value of xi was: {}'.format(xi))
    else:
        om=np.sqrt(k/m) # [rad/s] - undamped angular frequency
        oD=om*np.sqrt(1-xi**2) # [rad/s] - damped angular frequency
        
    # RESPONSE
        A=u0
        B=(v0+u0*xi*om)/oD
        rho=np.sqrt(A**2 + B**2)
        theta=np.arccos(A/rho)
        
    # OUTPUT
        uf=np.exp(-xi*om*t)*rho*np.cos(oD*t-theta)
        u={'uf':uf}
        cts={'om':om,'omD':oD,'xi':xi,'A':A,'B':B,'rho':rho,'theta':theta}
        msg= \
        "The response is of the form:\n" + \
        "uf(t)=exp(-xi om t) * [A cos(omD t) + B sin(omD t)]\n" + \
        "The values of the constants are:\n" + \
        "    om  = {} [rad/s]\n".format(om) + \
        "    omD = {} [rad/s]\n".format(oD) + \
        "    xi  = {} [-]\n".format(xi) + \
        "    A   = {} \n".format(A)  + \
        "    B   = {} \n".format(B)  + \
        "This can be put in the form:\n" + \
        "uf(t)=exp(-xi om t) * rho cos(omD t - theta)\n" + \
        "The values of the constants for this version are:\n" + \
        "    rho   = {} \n".format(rho) + \
        "    theta = {} [rad]\n".format(theta)
        
    return u,cts,msg


def uharm_underd(m,c,k,u0,v0,P0,Om,t=0):
    """
    Driven vibration response of a simple underdamped harmonic oscillator
    excited by a harmonic force defined as P(t)=P0 sin(Om t) and, possibly,
    non-zero initial conditions.
    Input: m,c,k,u0,v0,P0,Om,t=0
    Output:
        u:      dict with arrays with the evaluations of ut(t) and us(t)
                at the instants in array t
        cts:    dict of constants determining the response
        msg:    text explaining the output (incluiding the values)
    """
    
    # SYSTEM PROPERTIES
    c_cr=2*np.sqrt(m*k) # critical damping
    xi=c/c_cr # [-] - damping ratio
    if xi >= 1:
        raise Exception('\'uharm_underd\' is written for underdamped harmonic oscillators, \
                         thus xi=c/[2*sqrt(m*k)] should be lower than 1. \
                         The value of xi was: {}'.format(xi))
    else:
        om=np.sqrt(k/m) # [rad/s] - undamped angular frequency
        oD=om*np.sqrt(1-xi**2) # [rad/s] - damped angular frequency
        
    # RESPONSE
    # steady state
        beta=Om/om
        gamma=1/np.sqrt( (1-beta**2)**2 + (2*beta*xi)**2 )
        rho=P0*gamma/k
        theta=np.arccos( (1-beta**2)*gamma )
    # transient
        At=u0+rho*np.sin(theta)
        Bt= xi*om/oD * At  +  1/oD * ( v0+rho*Om**np.cos(theta) )
        
    # OUTPUT
        ut=np.exp(-xi*om*t)*( At*np.cos(oD*t) + Bt*np.sin(oD*t) )
        us=rho*np.sin(Om*t-theta)
        u={'ut':ut,'us':us}
        cts={'om':om,'omD':oD,'xi':xi,'At':At,'Bt':Bt,'rho':rho,'theta':theta,'beta':beta,'gamma':gamma}
        msg= \
        "The response is of the form:\n" + \
        "u(t)=ut(t) + us(t)  ,\n" + \
        "where ut(t) is the transient solution\n" + \
        "ut(t)=exp(-xi om t) * [At cos(omD t) + Bt sin(omD t)]  ,\n" + \
        "and us(t) is the steady state solution\n" + \
        "us(t)=                 As cos(Om  t) + Bs sin(Om  t)   ,\n" + \
        "usually written as\n" + \
        "us(t)=P0/k gamma sin(Om t - theta)=rho sin(Om t - theta).\n" + \
        "The values of the constants are:\n" + \
        "    om    = {} [rad/s]\n".format(om) + \
        "    omD   = {} [rad/s]\n".format(oD) + \
        "    xi    = {} [-]\n".format(xi) + \
        "    beta  = {} [-]\n".format(beta) + \
        "    gamma = {} \n".format(gamma) + \
        "    rho   = {} \n".format(rho) + \
        "    theta = {} [rad]\n".format(theta) + \
        "    At    = {} \n".format(At)  + \
        "    Bt    = {} \n".format(Bt)  
        
    return u,cts,msg


def uimpl_underd(m,k,u0,v0,ltype,P0,t0,t=0):
    """
    Driven vibration response of a simple underdamped harmonic oscillator
    excited by an impulsive idealised force and, possibly, non-zero initial
    conditions.
    Input: m,k,u0,v0,ltype,P0,t0,t=0
    Output:
        u:      dict with arrays with the evaluations of uf(t) (in case
                that non-zero initial conditions are given) and tuples
                (ti,ui(ti)), ti being the instants for interval 'i' in the
                definition of the load, and ui(ti) the evaluations of u
                at those intervals
    """
    
    # SYSTEM PROPERTIES
    om=np.sqrt(k/m) # [rad/s] - undamped angular frequency
        
    # RESPONSE
    u={}
    # forced
    if ltype == 1:
        pass
        
        
    elif ltype == 4: # constant-positive slope; constant
        
        # inner functs
        def u1(m,om,P0,t0,t1):
            u1=(1/(om*m)) * P0/(om*t0)*(t1-np.sin(om*t1)/om)
            return u1
        
        def u2(m,om,P0,t0,t2):
            u2=(1/(om*m)) * P0/(om*t0)*(  t0  -  np.sin( om*(t2-t0) )/om  -  np.sin(om*t2)/om  )
            return u2
        
        # solutions
        if type(t)==float:
            if t<=t0:
                t1=t
                u1=u1(m,om,P0,t0,t1)
                u['int1']=t1,u1
                u['all']=t1,u1
            else:
                t2=t
                u2=u2(m,om,P0,t0,t2)
                u['int2']=t2,u2
                u['all']=t2,u2
                
                
        else: # should be np.array
            try:
                t1=t[ np.where(t<=t0) ]
                u1=u1(m,om,P0,t0,t1)
                u['int1']=t1,u1
            except: 
                print('\'uimpl_underd\' with \'ltype=4\' couldn\'t calculate \'u\' for \'interval 1\'.')
                                 
            try:
                t2=t[ np.where(t>t0) ] # note that I didn't use '>=' (so it is easier to sum uf with u1 and u2)
                u2=u2(m,om,P0,t0,t2)
                u['int2']=t2,u2
            except: 
                print('\'uimpl_underd\' with \'ltype=4\' couldn\'t calculate \'u\' for \'interval 2\'.')
                
            u['all'] = np.hstack(( u['int1'], u['int2']))
        
        
    elif ltype == 12: # constant-negative slope; zero
        
        # inner functs
        def u1(m,om,P0,t0,t1):
            u1=(1/(om*m)) * P0/om*( 1 - np.cos(om*t1) - t1/t0                        + np.sin(om*t1)/(om*t0) )
            return u1
        
        def u2(m,om,P0,t0,t2):
            u2=(1/(om*m)) * P0/om*(   - np.cos(om*t2) - np.sin( om*(t2-t0) )/(om*t0) + np.sin(om*t2)/(om*t0) )
            return u2
        
        # solutions
        if type(t)==float:
            if t<=t0:
                t1=t
                u1=u1(m,om,P0,t0,t1)
                u['int1']=t1,u1
                u['all']=t1,u1
            else:
                t2=t
                u2=u2(m,om,P0,t0,t2)
                u['int2']=t2,u2
                u['all']=t2,u2
                
        else: # should be np.array
            try:
                t1=t[ np.where(t<=t0) ]
                u1=u1(m,om,P0,t0,t1)
                u['int1']=t1,u1
            except: 
                print('\'uimpl_underd\' with \'ltype=12\' couldn\'t calculate \'u\' for \'interval 1\'.')
                                 
            try:
                t2=t[ np.where(t>t0) ] # note that I didn't use '>=' (so it is easier to sum uf with u1 and u2)
                u2=u2(m,om,P0,t0,t2)
                u['int2']=t2,u2
            except: 
                print('\'uimpl_underd\' with \'ltype=12\' couldn\'t calculate \'u\' for \'interval 2\'.')
                
            u['all'] = np.hstack(( u['int1'], u['int2']))
        
        
    elif ltype == 13: # hat: constant-positive slope; constant-negative slope; zero
        
        # inner functs
        def u1(m,om,P0,t0,t1):
            u1=(1/(om*m)) * P0/(om*t0)*(        t1                                                     - np.sin(om*t1)/om )
            return u1
        
        def u2(m,om,P0,t0,t2):
            u2=(1/(om*m)) * P0/(om*t0)*( 2*t0 - t2 + 2*np.sin(om*(t2-t0))/om                           - np.sin(om*t2)/om )
            return u2
        
        def u3(m,om,P0,t0,t3):
            u3=(1/(om*m)) * P0/(om*t0)*(             2*np.sin(om*(t3-t0))/om - np.sin(om*(t3-2*t0))/om - np.sin(om*t3)/om )
            return u3
        
        # solutions
        if type(t)==float:
            if t<=t0:
                t1=t
                u1=u1(m,om,P0,t0,t1)
                u['int1']=t1,u1
                u['all']=t1,u1
            elif t<=2*t0:
                t2=t
                u2=u2(m,om,P0,t0,t2)
                u['int2']=t2,u2
                u['all']=t2,u2
            else:
                t3=t
                u3=u3(m,om,P0,t0,t3)
                u['int3']=t3,u3
                u['all']=t3,u3
                
        else: # should be np.array
            try:
                t1=t[ np.where(t<=t0) ]
                u1=u1(m,om,P0,t0,t1)
                u['int1']=t1,u1
            except: 
                print('\'uimpl_underd\' with \'ltype=13\' couldn\'t calculate \'u\' for \'interval 1\'.')
                                 
            try:
                # t2=t[ np.where(t>t0 and t<=2*t0) ] DOES NOT WORK
                t2=t[ np.where(t>t0) ] # note that I didn't use '>=' (so it is easier to sum uf with u1 and u2)
                t2=t2[ np.where(t2<=2*t0) ]
                u2=u2(m,om,P0,t0,t2)
                u['int2']=t2,u2
            except: 
                print('\'uimpl_underd\' with \'ltype=13\' couldn\'t calculate \'u\' for \'interval 2\'.')
                                 
            try:
                t3=t[ np.where(t>2*t0) ] # note that I didn't use '>=' (so it is easier to sum uf with u1 and u2)
                u3=u3(m,om,P0,t0,t3)
                u['int3']=t3,u3
            except: 
                print('\'uimpl_underd\' with \'ltype=13\' couldn\'t calculate \'u\' for \'interval 3\'.')
                
            u['all'] = np.hstack(( u['int1'], u['int2'], u['int3']))
        
        
    else:
        raise Exception('Load type {} not found.'.format(ltype))
        
        
    # non-zero initial conditions
    if u0!=0 or v0!=0:
        uf,cts,msg = ufree_underd(m,0,k,u0,v0,t)
        u['uf']=t,uf['uf']
        u['all'][1]=uf['uf'][1]+u['all'][1]
        
        
    return u


# RUN
if __name__ == '__main__':
    
    if True: # probar 'ufree_underd' y 'uharm_underd' y 'uimpl_underd'
        # v. libre
        m=0.1
        k=500
        c=(m+k)*0.0001
        u0=10
        v0=0
        t=np.arange( 0,10,0.0001, dtype=float)
        uf,cts,msg = ufree_underd(m,c,k,u0,v0,t)
        print(cts)
        print(msg)
        # plt.plot(t,uf['uf'])
        
        #v. forzada - armónica
        P0=1e3
        Om=0.8*cts['om']
        u,cts,msg = uharm_underd(m,c,k,u0,v0,P0,Om,t)
        print(cts)
        print(msg)
        # plt.plot(t,u['ut'])
        # plt.plot(t,u['us'])
        # plt.plot(t,u['ut']+u['us'])
        
        #v. forzada - impulsiva
        
        P0=1e3
        t0=1.5
        ltype=4
        uimp = uimpl_underd(m,k,u0,v0,ltype,P0,t0,t)
        # plt.plot(uimp['uf'][0],uimp['uf'][1])
        plt.plot(uimp['int1'][0],uimp['int1'][1])
        plt.plot(uimp['int2'][0],uimp['int2'][1])
        
        
    if False: 
        pass







