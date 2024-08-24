
"""
Example #1 for dsm.py
RC Hibbeler - Análisis Estructural - 8ed - Ejemplo 16.1
In SI units
Origin of coord. sys. on the bottom-left corner
"""

import numpy as np
from dsm import Struc2D3Dof

# conversion factors
# input
feet2mm = 304.8
ksi2MPa = 6.89476
in22mm2 = 645.16
in42mm4 = 416231.426
kip2N   = 4448.22
# output
mm2in   = 0.0393701

# creat problem
str=Struc2D3Dof('simple frame')

# nodes
str.create_node(1, np.array([ 0*feet2mm, 20*feet2mm], dtype=float))
str.create_node(2, np.array([20*feet2mm, 20*feet2mm], dtype=float))
str.create_node(3, np.array([20*feet2mm,  0*feet2mm], dtype=float))

# elms
E=29e3*ksi2MPa
A=10*in22mm2
I=500*in42mm4
str.create_elm(1, [1,2],  E, A, I)
str.create_elm(2, [2,3],  E, A, I)

# loads
str.create_NL(2, 0, 5*kip2N)

# BCs
str.create_BC(1, 1, 0)
str.create_BC(3, 0, 0)
str.create_BC(3, 1, 0)
str.create_BC(3, 2, 0)

# solve
str.assemble()
str.impose_BCs()
str.solve()

# print
print( 'Unknwon displacements:')
print(f'> D1: node 2, u_x = {str.nodes[1].dofVal[0]*mm2in} in')
print(f'> D2: node 2, u_y = {str.nodes[1].dofVal[1]*mm2in} in')
print(f'> D3: node 2, phi = {str.nodes[1].dofVal[2]} rad')
print(f'> D4: node 1, u_x = {str.nodes[0].dofVal[0]*mm2in} in')
print(f'> D5: node 1, phi = {str.nodes[0].dofVal[2]} rad')
print( '')
print( 'Reactions:')
idx=np.where( str.reacts.row == str.nodes[0].dofLab[1] )[0][0]
print(f'> Q6: node 1, R_y = { str.reacts.vtr[idx]/kip2N } kip')
idx=np.where( str.reacts.row == str.nodes[2].dofLab[0] )[0][0]
print(f'> Q7: node 3, R_x = { str.reacts.vtr[idx]/kip2N } kip')
idx=np.where( str.reacts.row == str.nodes[2].dofLab[1] )[0][0]
print(f'> Q8: node 3, R_y = { str.reacts.vtr[idx]/kip2N } kip')
idx=np.where( str.reacts.row == str.nodes[2].dofLab[2] )[0][0]
print(f'> Q9: node 3, Mom = { str.reacts.vtr[idx]/kip2N*mm2in } kip*in')
print( '')
print( 'Elemental nodal loads in global coordinates:')
print( 'element #1')
print(f'> elem 1, node i, L_xg = {str.elmts[0].L_glo[0]/kip2N} kip')
print(f'> elem 1, node i, L_yg = {str.elmts[0].L_glo[1]/kip2N} kip')
print(f'> elem 1, node i, mome = {str.elmts[0].L_glo[2]/kip2N*mm2in} kip*in')
print(f'> elem 1, node j, L_xg = {str.elmts[0].L_glo[3]/kip2N} kip')
print(f'> elem 1, node j, L_yg = {str.elmts[0].L_glo[4]/kip2N} kip')
print(f'> elem 1, node j, mome = {str.elmts[0].L_glo[5]/kip2N*mm2in} kip*in')
print( 'element #2')
print(f'> elem 2, node i, L_xg = {str.elmts[1].L_glo[0]/kip2N} kip')
print(f'> elem 2, node i, L_yg = {str.elmts[1].L_glo[1]/kip2N} kip')
print(f'> elem 2, node i, mome = {str.elmts[1].L_glo[2]/kip2N*mm2in} kip*in')
print(f'> elem 2, node j, L_xg = {str.elmts[1].L_glo[3]/kip2N} kip')
print(f'> elem 2, node j, L_yg = {str.elmts[1].L_glo[4]/kip2N} kip')
print(f'> elem 2, node j, mome = {str.elmts[1].L_glo[5]/kip2N*mm2in} kip*in')
