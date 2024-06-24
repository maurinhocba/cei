
"""
Direct Stifness Method
Mauro S. Maza - FCEFyN (UNC) - 2024-06-23
"""

"""
COSAS POR IMPLEMENTAR
- ElmFrame2D: asegurarse que "nodes" es una lista de dos nodos
- ElmFrame2D: asegurarse que "nodes" es una lista de punteros y no copias
- ElmFrame2D: incorporar la opción para articular extremos
- Node.add_dofs: asegurarse que el arreglo quede como vector columna
- Struc2D2Dof: asegurarse que "nodes" es una lista de punteros y no copias
- Struc2D2Dof: asegurarse que "elems" es una lista de punteros y no copias
- Struc2D2Dof: crear variables (o algún método) para llevar la cuenta de a qué GL corresponde cada fila/columna de cada matriz de rigidez
"""

import numpy as np
import matplotlib.pyplot as plt


# CLASSES
class Load:
    """
    general load
    """
    
    def __init__(self, extLab: int, type: int, **kwargs):
        # acá hay que implementar cosas con try except
        if type==1:
            if loaVal in kwars:
                self.loaVal = loaVal
            else:
                pass
        else:
            pass
    
    
class Node:
    """
    general node
    dimension defined when added to structure
    number of DoFs defined when added to structure (migth differ from dimension)
    """
    
    def __init__(self, extLab: int, coords: np.array([], dtype=float)):
        self.extLab = extLab
        self.coords = np.array([], dtype=float)
        
        self.intLab = int(0)                        # internal label, assigned when added to struct
        self.dofLab = np.array([], dtype=int)       # DoFs' labels, assigned when added to struct
        self.dofVal = np.array([], dtype=float)     # actual value of DoFs
        self.loads  = {}                            # dict {DOF: load}

    def add_dofs(self, dofLab: np.array([], dtype=int)):
        self.dofLab = dofLab
        self.dofVal = np.zeros( self.dofLab.shape )

    def add_load(self, dof: int, load: Load):
        self.loads[dof] = load

    def set_dofVal(self, dofVal: np.array([], dtype=float)):
        self.dofVal = dofVal


class ElmFrame2D:
    """
    2D Frame beam element.
    Not implemented yet:
        Might have one or both ends articulated;
        thus, might be a Truss bar.
    """
    
    def __init__(self, extLab: int, nodes: list, E: float, A: float, I: float, joints='none'):
        self.extLab = extLab
        self.nodes  = nodes     # list of nodes - should be pointers to node class instances
                                #   -->> además, habría que asegurarse que es una lista de dos nodos con dos coordenadas y tres GL <<--
        self.joints = joints    # one of 'none', 'first', 'second' and 'both'
        self.E      = E         # Young's modulus
        self.A      = A         # sectional area
        self.I      = I         # sectional inertia about z-axis
        
        self.L      = float(0)                      # length
        self.t      = np.array([], dtype=float)     # unit vector
        self.K_loc  = stiffM(self, local=True )     # stiffness matrix in local  coordinates
        self.K_glo  = stiffM(self, local=False)     # stiffness matrix in global coordinates
        
        self.dofLab = np.array([], dtype=int)       # DoFs' labels, assigned when added to struct
        
        self.length_and_uVectr()

    def length_and_uVectr(self):
        vect = self.nodes[1] - self.nodes[0]
        self.L = np.norm( vect )
        self.t = vect / self.L

    def stiffM(self, local=False):
        E = self.E
        A = self.A
        I = self.I
        L = self.L
        
        k  =    E*A/L   ;
        k1 = 12*E*I/L**3;
        k2 =  6*E*I/L**2;
        k3 =  4*E*I/L   ;
        
        if local:
            gamma1=1
            gamma2=0
        else:
            gamma1=self.t[0]
            gamma2=self.t[1]
        
        A = gamma1**2*k + gamma2**2*k1;
        B = gamma1*gamma2*(k-k1);
        C = gamma2*k2;
        D = gamma2**2*k + gamma1**2*k1;
        E = gamma1*k2;
        
        if joints=='none':
            # Stiffness matrix for displacements vector of the form:
            # U_ij = [  u_i^x  u_i^y  phi_i^z  |  u_j^x  u_j^y  phi_j^z  ]^T
            K = np.array( [ [ A,  B, -C,     -A, -B, -C   , ],
                            [ B,  D,  E,     -B, -D,  E   , ],
                            [-C,  E, k3,      C, -E,  k3/2, ],
                            #--------------------------
                            [-A, -B,  C  ,    A,  B,  C,    ],
                            [-B, -D, -E  ,    B,  D, -E,    ],
                            [-C,  E, k3/2,    C, -E, k3,    ] ] )
            return K
            
        else:
            print('ERROR: beams with joints not implemented yet - stiffness matrix determination of beam no. '+self.extLab)
            
    def set_DoF_labels(self):
        if joints=='none':
            # Stiffness matrix for displacements vector of the form:
            # U_ij = [  u_i^x  u_i^y  phi_i^z  |  u_j^x  u_j^y  phi_j^z  ]^T
            self.dofLab = np.block( [ self.nodes[0].dofLab(:) ]  ;  [ self.nodes[1].dofLab(:) ] )
        
        else:
            print('ERROR: beams with joints not implemented yet - stiffness matrix determination of beam no. '+self.extLab)
    

class Struc2D2Dof:
    """
    General bidimensional structure with 2 DoF per node
    """
    
    def __init__(self):
        self.nodes = [] # asegurarse de que esto es una lista de punteros
        self.elmts = [] # asegurarse de que esto es una lista de punteros
        self.ndofn = int(0)
        self.K_full = np.array([], dtype=float)
        self.K_wBCs = np.array([], dtype=float)
        self.K_cc   = np.array([], dtype=float)
        self.K_hc   = np.array([], dtype=float)
        self.K_hh   = np.array([], dtype=float)
        self.K_cond = np.array([], dtype=float)

    def add_node(self, node):
        if self.K_full.size > 0:
            print('ERROR: trying to add a node to a structure whose stiffness matriz has already been calculated.')
            exit()
        else:
            if node.coords.size==2:
                self.nodos.append(node)
                node.intLab=len(self.nodes)
                node.dofLab=np.arange( self.ndofn, self.ndofn+2 )
                self.ndofn = self.ndofn+2
            else:
                print('ERROR: trying to add a node with inconsistent number of coordinates (should be 2).')
                exit()

    def add_elm(self, elm):
        if self.K_full.size > 0:
            print('ERROR: trying to add an element to a structure whose stiffness matriz has already been calculated.')
            exit()
        else:
            if isinstance( ElmFrame2D, elm):
                # VERIFICAR SI LOS NODOS DEL ELEMENTO EXISTEN EN LA ESTRUCTURA
                self.elmts.append(elm)
                elm.set_DoF_labels()
            else:
                print('ERROR: element type not supported - element no. '+elm.extLab)

    def ensamblar_matriz_rigidez(self):
        n_dofs = sum([len(nodo.grados_libertad) for nodo in self.nodos])
        self.matriz_rigidez_global = np.zeros((n_dofs, n_dofs))
        # Ensamblaje de la matriz de rigidez global
        for elemento in self.elementos:
            k_local = elemento.rigidez
            # Mapeo de índices locales a globales
            # Agregar k_local a la posición correcta en matriz_rigidez_global
            pass

    def aplicar_condiciones_frontera(self, condiciones):
        # Modificación de la matriz y vector de cargas según las condiciones de frontera
        pass

    def resolver(self):
        # Implementación del solucionador del sistema de ecuaciones
        pass
