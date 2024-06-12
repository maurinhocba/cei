

import numpy as np
import matplotlib.pyplot as plt
import sho
import mdof


# DATOS **************
# tiempo para evaluación de la respuesta
tt=np.linspace(0, 3, num=3000);

# caracterización dnámica de la estructura
# viga
h=10 # mm - altura de la sección
I=h**4/12
E=210e3
Kcv = np.array( [ [ 0,          0 ],
                  [ 0, E*I/26.7e6 ] ] )
# resorte
k=20
alpha=np.deg2rad(60)
g1= np.cos(alpha)
g2=-np.sin(alpha)
Kcr = np.array( [ [ k*g1**2, k*g1*g2 ],
                  [ k*g1*g2, k*g2**2 ] ] )
# todo
Kc = Kcv+Kcr
m=6e-3
Mc=np.diag((m,m))

# carga impulsiva
P_tipo=12
P_imp=np.array([0,-50], dtype=float)
t0=0.85
u0_imp=np.array([0,0], dtype=float)



# SOLUCIÓN **************
# modos
PHI, om = mdof.modes(Kc,Mc)

# CARGA IMPULSIVA
# sist. diagonalizado
eqmot_imp, ic_imp = mdof.uncoup(PHI,Kc,Mc,P_imp,u0_imp)
Tq1=2*np.pi/om[0] # período del modo 1
u_est_q1=eqmot_imp['Pb'][0]/eqmot_imp['Kb'][0] # desplazamiento estático de refe modo 1

# solución de los sociladores simples
q1_imp = sho.uimpl_underd( eqmot_imp['Mb'][0], eqmot_imp['Kb'][0], ic_imp['q0'][0], ic_imp['dq0'][0], P_tipo, eqmot_imp['Pb'][0], t0, t=tt)
q2_imp = sho.uimpl_underd( eqmot_imp['Mb'][1], eqmot_imp['Kb'][1], ic_imp['q0'][1], ic_imp['dq0'][1], P_tipo, eqmot_imp['Pb'][1], t0, t=tt)
# plt.plot( q1_imp['uf'][0], q1_imp['uf'][1])
# plt.plot( q1_imp['int1'][0], q1_imp['int1'][1])
# plt.plot( q1_imp['int2'][0], q1_imp['int2'][1])
# plt.plot( q1_imp['all'][0], q1_imp['all'][1])

# reconstrucción de la solución
#   coordenadas geométricas
q_imp=np.vstack(( q1_imp['all'][1], q2_imp['all'][1]))
u_imp=np.matmul( PHI, q_imp )
# plt.plot( tt, u_imp[0])
# plt.plot( tt, u_imp[1])


# FIGURA
fig, axs = plt.subplots( 2, 1, sharex=True)
# arriba - coordenadas modales
axs[0].plot( q1_imp['all'][0], q1_imp['all'][1])
axs[0].plot( q2_imp['all'][0], q2_imp['all'][1])
axs[0].grid()
axs[0].set_ylabel("q(t)")
axs[0].legend(['q1(t)', 'q2(t)'])

# abajo - coordenadas geométricas
axs[1].plot( tt, u_imp[0])
axs[1].plot( tt, u_imp[1])
axs[1].grid()
axs[1].set_ylabel("u [mm]")
axs[1].set_xlabel("t [s]")
axs[1].legend(['u1(t)', 'u2(t)'])

# guardar
# fig.savefig("sols.eps")

