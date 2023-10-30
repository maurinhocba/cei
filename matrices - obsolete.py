
import numpy as np

def mrpp(A,E,L,I,gamma1,gamma2):
    
    # Calcula la matriz de rigidez de una barra de pórtico plano, como se ve en Cálculo Estructural I.
    # El vector de desplazamientos asociado es:
    #             U_ij = [  u_i^x  u_i^y  phi_i^z  |  u_j^x  u_j^y  phi_j^z  ]^T
    
    k  =    E*A/L   ;
    k1 = 12*E*I/L**3;
    k2 =  6*E*I/L**2;
    k3 =  4*E*I/L   ;
    
    A = gamma1**2*k + gamma2**2*k1;
    B = gamma1*gamma2*(k-k1);
    C = gamma2*k2;
    D = gamma2**2*k + gamma1**2*k1;
    E = gamma1*k2;
    
    K = np.array( [ [ A,  B, -C,     -A, -B, -C   , ],
                    [ B,  D,  E,     -B, -D,  E   , ],
                    [-C,  E, k3,      C, -E,  k3/2, ],
                    #--------------------------
                    [-A, -B,  C  ,    A,  B,  C,    ],
                    [-B, -D, -E  ,    B,  D, -E,    ],
                    [-C,  E, k3/2,    C, -E, k3,    ] ] )
    
    return K

def mrep(E,G,L,I,J,gamma1,gamma2):
    
    # Calcula la matriz de rigidez de una barra de emparrillado plano, como se ve en Cálculo Estructural I.
    # El vector de desplazamientos asociado es:
    #             U_ij = [  u_i^z  phi_i^x  phi_i^y  |  u_j^z  phi_j^x  phi_j^y  ]^T
    
    k  = G*J/L;
    k1 = 12.0*E*I/L**3;
    k2 =  6.0*E*I/L**2;
    k3 =  4.0*E*I/L  ;
    
    A = gamma2*k2;
    B = -gamma1*k2;
    C = gamma1**2*k + gamma2**2*k3;
    D = gamma1*gamma2*( k - k3 );
    E = gamma2**2*k + gamma1**2*k3;
    F = -gamma1*gamma2*( k + k3/2 );
    G = -gamma1**2*k + gamma2**2*k3/2;
    H = -gamma2**2*k + gamma1**2*k3/2;
    
    K = np.array( [ [  k1,  A,  B,     -k1,  A,  B, ],
                    [   A,  C,  D,      -A,  G,  F, ],
                    [   B,  D,  E,      -B,  F,  H, ],
                    #---------------------------------
                    [ -k1, -A, -B,      k1, -A, -B, ],
                    [   A,  G,  F,      -A,  C,  D, ],
                    [   B,  F,  H,      -B,  D,  E, ] ] )
    
    return K
    
def mrotep1v2c(g1,g2):
    
    # Para un vector de 2 componentes v=[v1,v2]^T como u=[phi_x,phi_y]
    # Permite obtener la expresión en coord. globales a partir de la expresión en coordenadas locales
    # v_g=R*v_l siendo g1 y g2 los elementos del versor que define el eje local x
    
    R = np.array( [ [  g1,  -g2, ],
                    [  g2,   g1, ] ] )
    
    return R
    
def mrotep1v3c(g1,g2):
    
    # Para un vector de 3 componentes v=[v1,v2,v3]^T como u=[u_z,phi_x,phi_y]
    # Permite obtener la expresión en coord. globales a partir de la expresión en coordenadas locales
    # v_g=R*v_l siendo g1 y g2 los elementos del versor que define el eje local x
    
    R = np.array( [ [  1,   0,    0, ],
                    [  0,  g1,  -g2, ],
                    [  0,  g2,   g1, ] ] )
    
    return R
    
def mrotep2v3c(g1,g2):
    
    # Para un vector de 6 componentes v=[v1,v2,v3,v4,v5,v6]^T como u=[u_z^i,phi_x^i,phi_y^i , u_z^j,phi_x^j,phi_y^j]
    # Permite obtener la expresión en coord. globales a partir de la expresión en coordenadas locales
    # v_g=R*v_l siendo g1 y g2 los elementos del versor que define el eje local x
    
    R = np.array( [ [   1,  0,   0,      0,  0,  0, ],
                    [   0, g1, -g2,      0,  0,  0, ],
                    [   0, g2,  g1,      0,  0,  0, ],
                    #---------------------------------
                    [   0,  0,  0,       1,  0,   0, ],
                    [   0,  0,  0,       0, g1, -g2, ],
                    [   0,  0,  0,       0, g2,  g1, ] ] )
    
    return R