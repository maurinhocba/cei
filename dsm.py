
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

# IMPORTING ZONE
import numpy as np
# import matplotlib.pyplot as plt


# CLASSES
class Load:
    """
    general load
    """
    
    def __init__(self, extLab: int, type: int, **kwargs):
        # acá hay que implementar cosas con try except
        if type==1:
            if 'loaVal' in kwargs.keys():
                self.loaVal = kwargs['loaVal']
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
        self.coords = coords
        
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


class dpm:
    """
    dynamic property matrix
    has a matrix and
    a vector indicating the numbers of equations to which the rows and cols refer to
    """
    
    def __init__(self):
        self.mtx = np.array([], dtype=float)
        self.dof = np.array([], dtype=int)
        
    def set_mtx_dof(self, mtx, dof):
        if mtx.ndim != 2:
            print('ERROR: trying to set a dynamic property matrix with more or less than 2 dimensions.')
        else:
            if mtx.shape[0]!=mtx.shape[1]:
                print('ERROR: trying to set a dynamic property matrix with different number of rows and cols.')
            else:
                if dof.ndim!=1 or dof.shape[0]!=mtx.shape[0]:
                    print('ERROR: trying to set a dynamic property matrix with a wrong DoF vector.')
                else:
                    self.mtx = mtx
                    self.dof = dof


class ElmFrame2D:
    """
    2D Frame beam element.
    Not implemented yet:
        Might have one or both ends articulated;
        thus, might be a Truss bar.
    """
    
    def __init__(self, extLab: int, nodes: list, E: float, A: float, I: float, joints='none'):
        self.extLab = extLab
        self.nodes  = nodes     # list of external labels of nodes
        self.joints = joints    # one of 'none', 'first', 'second' and 'both'
        self.E      = E         # Young's modulus
        self.A      = A         # sectional area
        self.I      = I         # sectional inertia about z-axis
        
        self.ndofpn = []                        # list - number of DoFs per node
        self.dofpn  = ''                        # string - DoFs per node
        self.L      = float(0)                  # length
        self.t      = np.array([], dtype=float) # unit vector
        self.K_loc  = dpm()                     # stiffness matrix in local  coordinates
        self.K_glo  = dpm()                     # stiffness matrix in global coordinates
        self.dofLab = np.array([], dtype=int)   # DoFs' labels, assigned when added to struct
    
    def dofPerNode(self):
        if self.joints=='none':
            # Stiffness matrix for displacements vector of the form:
            # U_ij = [  u_i^x  u_i^y  phi_i^z  |  u_j^x  u_j^y  phi_j^z  ]^T
            self.ndofpn=[3, 3]
            self.dofpn='[  u_i^x  u_i^y  phi_i^z  |  u_j^x  u_j^y  phi_j^z  ]^T'
        elif self.joints=='first':
            # Stiffness matrix for displacements vector of the form:
            # U_ij = [  u_i^x  u_i^y  |  u_j^x  u_j^y  phi_j^z  ]^T
            self.ndofpn=[2, 3]
            self.dofpn='[  u_i^x  u_i^y  |  u_j^x  u_j^y  phi_j^z  ]^T'
        elif self.joints=='second':
            # Stiffness matrix for displacements vector of the form:
            # U_ij = [  u_i^x  u_i^y  phi_i^z  |  u_j^x  u_j^y  ]^T
            self.ndofpn=[3, 2]
            self.dofpn='[  u_i^x  u_i^y  phi_i^z  |  u_j^x  u_j^y  ]^T'
        elif self.joints=='second':
            # Stiffness matrix for displacements vector of the form:
            # U_ij = [  u_i^x  u_i^y  |  u_j^x  u_j^y  ]^T
            self.ndofpn=[2, 2]
            self.dofpn='[  u_i^x  u_i^y  |  u_j^x  u_j^y  ]^T'
        else:
            print('ERROR: undefined flag for nodes: joints='+self.joints)
    
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
        
        if self.joints=='none':
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
    
    def detStiff(self):
        self.length_and_uVectr()
        self.K_loc  = self.stiffM(self, local=True )     # stiffness matrix in local  coordinates
        self.K_glo  = self.stiffM(self, local=False)     # stiffness matrix in global coordinates
            
    """def set_DoF_labels(self):
        if joints=='none':
            # Stiffness matrix for displacements vector of the form:
            # U_ij = [  u_i^x  u_i^y  phi_i^z  |  u_j^x  u_j^y  phi_j^z  ]^T
            self.dofLab = np.block( [ self.nodes[0].dofLab(:) ]  ;  [ self.nodes[1].dofLab(:) ] )
        
        else:
            print('ERROR: beams with joints not implemented yet - stiffness matrix determination of beam no. '+self.extLab)
    """
    

class Struc2D3Dof:
    """
    General bidimensional structure with 3 DoF per node
    """
    
    def __init__(self, name='New Project'):
        selt.name  = name
        # self.nodes = []
        self.nodes = {}
        self.elmts = []
        self.ndofn = int(0)
        self.K_full = dpm()
        self.K_wBCs = dpm()
        self.K_cc   = dpm()
        self.K_hc   = dpm()
        self.K_hh   = dpm()
        self.K_cond = dpm()

    def add_node(self, node: Node):
        
        # VERIFICAR QUE LA ETIQUETA EXTERNA NO EXISTA
        
        if self.K_full.size > 0:
            print('ERROR: trying to add a node to a structure whose stiffness matrix has already been calculated.')
            exit()
        else:
            if node.coords.size==2:
                # self.nodos.append(node)
                self.nodes[node]=[]
                node.intLab=len(self.nodes)
            else:
                print('ERROR: trying to add a node with inconsistent number of coordinates (should be 2).')
                exit()

    def add_elm(self, elm):
        
        # VERIFICAR QUE LA ETIQUETA EXTERNA NO EXISTA
        
        if self.K_full.size > 0:
            print('ERROR: trying to add an element to a structure whose stiffness matrix has already been calculated.')
            exit()
        else:
            if isinstance( ElmFrame2D, elm):
                # verify that the nodes defining the beam exist in the structure
                verif=[]
                i=-1
                for targetNode in elm.nodes:
                    verif.append(False)
                    i=i+1
                    for possibleNode in self.nodes.keys():
                        if targetNode==possibleNode.extLab:
                            # self.nodes[possibleNode].append([elm.extLab, i, elm.ndofpn[i])
                            verif[-1]=True
                            break
                if all(verif):
                    self.elmts.append(elm)
                    elm.detStiff()
                else:
                    print('ERROR: nodes not in structure - beam no. '+elm.extLab)
                    exit()
            else:
                print('ERROR: element type not supported - element no. '+elm.extLab)

    def eqnums(self):
        # determine number of eq. per node and their labels
        for node, data in self.nodes.items():
            ndof=0
            for lst in data:
                if ndof<data[-1]:
                    ndof=data[-1]
            node.dofLab=np.arange( self.ndofn, self.ndofn+ndof-1 )
            self.ndofn = self.ndofn+ndof
        
    def ensamblar_matriz_rigidez(self):
        # determine number of eq. per node and their labels 
        self.eqnums()
        
        # prepare full matrix for assembling
        self.K_full.mtx = np.zeros((self.ndofn, self.ndofn))
        self.K_full.dof = np.arange(0, self.ndofn-1)
        
        # assembly
        for elm in self.elmts:
            # assing DoFs (eqs) numbers to element
            for targetNode in elm.nodes:
                for possibleNode in self.nodes.keys():
                    if targetNode==possibleNode.extLab:
                        elm.dofLab=np.concatenate((elm.dofLab,possibleNode.dofLab),axis=0)
                        break
            
            # assemble
            # iterate over element's DoFs
            it = np.nditer(elm.dofLab, flags=['f_index'])
            # and assemble the element matrix one row at a time
            for dof in it:
                self.K_full.mtx[dof,elm.dofLab]=elm.K_glo[it.index,:]

    def aplicar_condiciones_frontera(self, condiciones):
        # Modificación de la matriz y vector de cargas según las condiciones de frontera
        pass

    def resolver(self):
        # Implementación del solucionador del sistema de ecuaciones
        pass


# RUN
if __name__ == '__main__':
    
    if True:
        n1=Node(10, np.array([0,0  ], dtype=float))
        n2=Node(20, np.array([0,5.5], dtype=float))
        n3=Node(30, np.array([3,0  ], dtype=float))
        n4=Node(40, np.array([3,5.5], dtype=float))
        
        elm1=ElmFrame2D(1, [10,20], 210e3, 100, 10**4/12)
        elm2=ElmFrame2D(2, [10,40], 210e3, 100, 10**4/12)
        elm3=ElmFrame2D(3, [20,40], 210e3, 100, 10**4/12, joints='both')
    
    
    if False:
        pass
    
    
    