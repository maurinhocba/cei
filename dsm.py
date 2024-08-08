
"""
Direct Stifness Method
Mauro S. Maza - FCEFyN (UNC) - 2024-06-23
"""

"""
COSAS POR IMPLEMENTAR
- ElmFrame2D: incorporar la opción para articular extremos
- Node.add_dofs: asegurarse que el arreglo quede como vector columna    
        # IMPLEMENTAR MÉTODO PARA AGREGAR MUCHOS NODOS Y MUCHOS ELEMENTOS A LA VEZ
        # IMPLEMENTAR MÉTODO PARA ACTUALIZAR COORDENADAS DE UN PUNTO
        # PROBAR SI AL CAMBIAR LAS COORDENADA DE UN NODO, SE ACTUALIZA LA REFERENCIA EN EL ELEMENTO
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
        self.dofLab = np.array([], dtype=int)       # DoFs' labels, assigned befor assembly
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
    two vectors indicating the numbers of equations to which the rows and cols refer to
    """
    
    def __init__(self):
        self.mtx = np.array([], dtype=float)
        self.row = np.array([], dtype=int)
        self.col = np.array([], dtype=int)
        
    def set_mtx_dof(self, mtx, row, col):
        if mtx.ndim != 2:
            print('ERROR: trying to set a dynamic property matrix with more or less than 2 dimensions.')
        else:
            if row.ndim!=1 or row.shape[0]!=mtx.shape[0] or col.ndim!=1 or col.shape[0]!=mtx.shape[1]:
                print('ERROR: trying to set a dynamic property matrix with wrong DoF vectors.')
            else:
                self.mtx = mtx
                self.row = row
                self.col = col


class ElmFrame2D:
    """
    2D Frame beam element.
    Not implemented yet:
        Might have one or both ends articulated;
        thus, might be a Truss bar.
    """
    
    def __init__(self, extLab: int, nodes: list, E: float, A: float, I: float, joints='none'):
        self.extLab = extLab
        self.nodes  = nodes     # list of nodes
        self.joints = joints    # one of 'none', 'first', 'second' and 'both'
        self.E      = E         # Young's modulus
        self.A      = A         # sectional area
        self.I      = I         # sectional inertia about z-axis
        
        self.ndofpn = []                        # list - number of DoFs per node
        self.dofpn  = ''                        # string - DoFs per node
        self.L      = float(0)                  # length
        self.t      = np.array([], dtype=float) # unit vector
        self.K_loc  = np.array([], dtype=float) # stiffness matrix in local  coordinates
        self.K_glo  = np.array([], dtype=float) # stiffness matrix in global coordinates
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
            print(f'ERROR: "{self.joints}" is an undefined flag for joints')
    
    def length_and_uVectr(self):
        vect = self.nodes[1].coords - self.nodes[0].coords
        self.L = np.linalg.norm( vect )
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
            print( 'ERROR: beams with joints not implemented yet\n'
                  f'       stiffness matrix determination of beam no. {self.extLab}')
    
    def detStiff(self):
        self.length_and_uVectr()
        self.K_loc  = self.stiffM(local=True )     # stiffness matrix in local  coordinates
        self.K_glo  = self.stiffM(local=False)     # stiffness matrix in global coordinates
            
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
        self.name  = name
        self.desc  = ''
        self.nodes = {}         # {extLab:[node_object]}
        self.elmts = {}
        self.eqLab = {}         # {eqLab:[node_object local_DoF_0->2]}
        self.ndofn = int(0)
        self.K_full = dpm()
        self.K_wBCs = dpm()
        self.K_cc   = dpm()
        self.K_hc   = dpm()
        self.K_hh   = dpm()
        self.K_cond = dpm()

    def add_node(self, node: Node):
        if self.K_full.mtx.size > 0:
            print('WARNING: trying to add a node to a structure\n'
                  '         whose stiffness matrix has already been calculated.\n'
                  '         RECALC')
        else:
            if node.coords.size==2:
                if node.extLab in self.nodes.keys():
                    print( 'ERROR: trying to add a node with an existing external label\n'
                           '       node not added\n'
                          f'       node external label: {node.extLab}')
                else:
                    self.nodes[node.extLab]=[node]
                    node.intLab=len(self.nodes)
            else:
                print( 'ERROR: trying to add a node with inconsistent number of coordinates (should be 2)\n'
                       '       node not added\n'
                      f'       node external label: {node.extLab}')
    
    def create_node(self, extLab: int, coords: np.array([], dtype=float)):
        node=Node(extLab, coords)
        self.add_node(node)
    
    def add_elm(self, elm):
        if self.K_full.mtx.size > 0:
            print('WARNING: trying to add an element to a structure\n'
                  '         whose stiffness matrix has already been calculated.\n'
                  '         RECALC')
        else:
            if isinstance( elm, ElmFrame2D):
                # verify that the nodes defining the beam exist in the structure
                verif=[]
                for targetNode in elm.nodes:
                    if targetNode.extLab in self.nodes.keys():
                        verif.append(True)
                    else:
                        verif.append(False)
                        
                if all(verif):
                    if elm.extLab in self.elmts.keys():
                        print( 'ERROR: trying to add an elem with an existing external label\n'
                               '       elem not added\n'
                              f'       elem external label: {elm.extLab}')
                    else:
                        self.elmts[elm.extLab]=[elm]
                        elm.detStiff()
                else:
                    print( 'ERROR: trying to add an elem whose defining nodes are not in the structure\n'
                           '       elem not added\n'
                          f'       elem external label: {elm.extLab}')
            else:
                print( 'ERROR: element type not supported\n'
                      f'       elem external label: {elm.extLab}')

    def create_elm(self, extLab: int, node_labs: list, E: float, A: float, I: float, joints='none'):
        """
        arg "node_labs" here is a list of external labels of nodes
        which is different from arg "nodes" in ElmFrame2D.__init__() method
        """
        
        # verify that the nodes defining the elem exist in the structure
        # and create the list of node objects for ElmFrame2D.__init__()
        verif=[]
        nodes=[]
        for targetNode in node_labs:
            if targetNode in self.nodes.keys():
                verif.append(True)
                nodes.append( self.nodes[targetNode][0] )
            else:
                verif.append(False)
                
        if all(verif):
            elm=ElmFrame2D(extLab, nodes, E, A, I, joints)
            self.add_elm(elm)
        else:
            print( 'ERROR: trying to add an elem whose defining nodes are not in the structure\n'
                   '       elem not added\n'
                  f'       elem external label: {elm.extLab}')
    
    def eqnums(self):
        # determine an eq. label for each node's DoF
        for node in self.nodes.values():
            # assign eq. numbers to node's DoFs
            node.dofLab=np.arange( self.ndofn, self.ndofn+3 )
            # constructu dict. for labeling translation
            for idx, eqLab in enumerate(node.dofLab)
                self.eqLab[eqLab]=[node idx]
            # update number of DoFs in problem
            self.ndofn = self.ndofn+3
            
        # assign eq. labels to elements
        for elm in self.elmts:
            for node in elm.nodes:
                elm.dofLab=np.concatenate((elm.dofLab,targetNode.dofLab),axis=0)
            # for targetNode in elm.nodes:
                # for possibleNode in self.nodes.keys():
                    # if targetNode.extLab==possibleNode:
                        # elm.dofLab=np.concatenate((elm.dofLab,targetNode.dofLab),axis=0)
                        # break
                # else:
                    # continue
                # break
        
    def assemble_K(self):
        # determine number of eq. per node and their labels 
        self.eqnums()
        
        # prepare full matrix for assembling
        mtx = np.zeros((self.ndofn, self.ndofn))
        dof = np.arange(0, self.ndofn-1)
        self.K_full.set_mtx_dof(mtx, row=dof, col=dof)
        
        # assembly
        for elm in self.elmts: 
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
    
    # simple structure
    if True:
        # create structure object
        str=Struc2D3Dof('simple frame')
        
        # creat node and elem objects and add them to the structure
        if True: # in one step - recomended
            str.create_node(10, np.array([0,0  ], dtype=float))
            str.create_node(20, np.array([0,5.5], dtype=float))
            str.create_node(30, np.array([3,0  ], dtype=float))
            str.create_node(40, np.array([3,5.5], dtype=float))
            
            str.create_elm(1, [10,20], 210e3, 100, 10**4/12)
            str.create_elm(2, [10,40], 210e3, 100, 10**4/12)
            str.create_elm(3, [20,40], 210e3, 100, 10**4/12, joints='both')
            
        else: # in two steps - not recomended
            n1=Node(10, np.array([0,0  ], dtype=float))
            n2=Node(20, np.array([0,5.5], dtype=float))
            n3=Node(30, np.array([3,0  ], dtype=float))
            n4=Node(40, np.array([3,5.5], dtype=float))
            str.add_node(n1)
            str.add_node(n2)
            str.add_node(n3)
            str.add_node(n4)
        
            elm1=ElmFrame2D(1, [10,20], 210e3, 100, 10**4/12)
            elm2=ElmFrame2D(2, [10,40], 210e3, 100, 10**4/12)
            elm3=ElmFrame2D(3, [20,40], 210e3, 100, 10**4/12, joints='both')
            str.add_elm(elm1)
            str.add_elm(elm2)
            str.add_elm(elm3)
        
        # assemble global stiffness matrix
        
    
    if False:
        pass
    
    
    