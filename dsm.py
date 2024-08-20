
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
from typing import List
import copy
# import matplotlib.pyplot as plt


# CLASSES
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
        self.dofLab = np.array([], dtype=int)       # global DoFs' labels, assigned befor assembly
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
        self.row = np.array([], dtype=int)      # labels of global DoFs associated to each row of matrix
        self.col = np.array([], dtype=int)      # labels of global DoFs associated to each col of matrix
        
    def set_mtx_dof(self, mtx: np.array([], dtype=float),
                          row: np.array([], dtype=int  ),
                          col: np.array([], dtype=int  )):
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
    
    def __init__(self, extLab: int, nodes: List[Node], E: float, A: float, I: float, joints='none'):
        
        if not all(isinstance(node, Node) for node in nodes):
            raise ValueError( 'ERROR: ElmFrame2D 2nd arg. shall be a list of 2 Node objects\n'
                             f'       trying to create elem no. {extLab}')
            
        self.extLab = extLab
        self.nodes  = nodes     # list of node objects
        self.joints = joints    # one of 'none', 'first', 'second' and 'both'
        self.E      = E         # Young's modulus
        self.A      = A         # sectional area
        self.I      = I         # sectional inertia about z-axis
        
        self.intLab = int(0)                    # internal label, assigned when added to struct
        self.ndofpn = []                        # list - number of DoFs per node
        self.dofpn  = ''                        # string - DoFs per node
        self.L      = float(0)                  # length
        self.t      = np.array([], dtype=float) # unit vector
        self.K_loc  = np.array([], dtype=float) # stiffness matrix in local  coordinates
        self.K_glo  = np.array([], dtype=float) # stiffness matrix in global coordinates
        self.dofLab = np.array([], dtype=int)   # DoFs' labels, assigned on assembly
    
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
            print(f'ERROR: "{self.joints}" is an undefined flag for joints\n'
                  f'       determination of number of DoF per node of elem no. {self.extLab}')
    
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
            # if self.joints!='none', then self.dofLab should be updated (eliminating elements)
            print( 'ERROR: beams with joints not implemented yet\n'
                  f'       stiffness matrix determination of elem no. {self.extLab}')
    
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
    under static load system
    """
    
    def __init__(self, name='New Project'):
        self.name  = name
        self.desc  = ''
        self.nodes = []         # list of node objects
        self.elmts = []         # list of elem objects
        self.eqLab = {}         # {eqLab:[node_object, local_DoF_0->2]}
        self.ndofn = int(0)
        
        # stiffness matrices
        self.K_full = dpm()     # assembled global matrix (without BCs)
        self.K_wBCs = dpm()     # assembled global matrix with BCs
        self.K_BCLo = dpm()     # columns of K_full associated to (zero and non-zero) BCs (added to load vector)
        self.K_cc   = dpm()
        self.K_hc   = dpm()
        self.K_hh   = dpm()
        self.K_cond = dpm()
        
        # BCs aligned to global coordinate system
        self.BC_dict = {}                           # dictionary for BCs input {node_object:[local_Dof_label, value]}
        self.BC_labs = np.array([], dtype=int)      # labels of global DoFs with imposed value
        self.BC_vals = np.array([], dtype=float)    # values imposed
        
        # self.loads = []         # [ [node_extLab DoF_locLab value] ]

    def add_node(self, newNode: Node):
        # verify that the node does not exist yet
        if any(newNode is oldNode for oldNode in self.nodes):
            print( 'ERROR: trying to add a node that already exists in the structure\n'
                   '       node not added\n'
                  f'       node external label: {newNode.extLab}')
        else:
            if newNode.coords.size==2:
                # verify that the external label of the node does not exist yet
                if newNode.extLab in [oldNode.extLab for oldNode in self.nodes]:
                #if any(newNode.extLab == oldNode.extLab for oldNode in self.nodes):
                    print( 'ERROR: trying to add a node with an existing external label\n'
                           '       node not added\n'
                          f'       node external label: {newNode.extLab}')
                else:
                    newNode.intLab=len(self.nodes)
                    self.nodes.append(newNode)
                    if self.K_full.mtx.size > 0:
                        print('WARNING: trying to add a node to a structure\n'
                              '         whose global stiffness matrix has already been calculated.\n'
                              '         RECALC')
            else:
                print( 'ERROR: trying to add a node with inconsistent number of coordinates (should be 2)\n'
                       '       node not added\n'
                      f'       node external label: {newNode.extLab}')
    
    def create_node(self, extLab: int, coords: np.array([], dtype=float)):
        node=Node(extLab, coords)
        self.add_node(node)
    
    def add_elm(self, newElm):
        if isinstance( newElm, ElmFrame2D):
            # verify that the elem does not exist yet
            if any(newElm is oldElm for oldElm in self.elmts):
                print( 'ERROR: trying to add an elem that already exists in the structure\n'
                       '       elem not added\n'
                      f'       elem external label: {newElm.extLab}')
            else:
                # verify that the external label of the elem does not exist yet
                if newElm.extLab in [oldElm.extLab for oldElm in self.elmts]:
                    print( 'ERROR: trying to add an elem with an existing external label\n'
                           '       elem not added\n'
                          f'       elem external label: {newElm.extLab}')
                else:
                    # verify that the nodes defining the elem exist in the structure
                    if all(  any(newElmNode is struNode  for struNode in self.nodes)   for newElmNode in newElm.nodes):
                        newElm.intLab=len(self.elmts)
                        self.elmts.append(newElm)
                        if self.K_full.mtx.size > 0:
                            print('WARNING: trying to add an element to a structure\n'
                                  '         whose global stiffness matrix has already been calculated.\n'
                                  '         RECALC')
                    else:
                        print( 'ERROR: trying to add an elem whose defining nodes are not in the structure\n'
                               '       elem not added\n'
                              f'       elem external label: {newElm.extLab}')
        else:
            print( 'ERROR: element type not supported\n'
                  f'       elem external label: {newElm.extLab}')

    def create_elm(self, extLab: int, node_labs: list, E: float, A: float, I: float, joints='none'):
        """
        arg "node_labs" here is a list of external labels of nodes
        which is different from arg "nodes" in ElmFrame2D.__init__() method
        """
        
        # verify that the nodes defining the elem exist in the structure
        # and create the list of node objects for ElmFrame2D.__init__()
        verif=[]
        nodes=[]
        for targetLab in node_labs:
            for posibleNode in self.nodes:
                if targetLab==posibleNode.extLab:
                    verif.append(True)
                    nodes.append( posibleNode )
                    break
            else:
                verif.append(False)
                break
                
        if all(verif):
            elm=ElmFrame2D(extLab, nodes, E, A, I, joints)
            self.add_elm(elm)
        else:
            print( 'ERROR: trying to add an elem whose defining nodes are not in the structure\n'
                   '       elem not added\n'
                  f'       elem external label: {elm.extLab}')
    
    def create_BC(self, nodeLab: int, locDof: int, val: float):
        """
        nodeLab: external label of node
        locDof: local DoF label
        val: value imposed        
        """
        
        # verify that the node does exist
        for posibleNode in self.nodes:
            if nodeLab==posibleNode.extLab:
                if int in range(0,3):
                    self.BC_dict[posibleNode] = [locDof, val]   # dictionary for BCs input {node_object:[local_Dof_label, value]}
                    if self.K_wBCs.mtx.size > 0:
                        print('WARNING: trying to add a BC to a structure\n'
                              '         whose BCs has already been imposed.\n'
                              '         RECALC')
                    break
                else:
                    print( 'ERROR: trying to set a BC on a non-existent local Dof\n'
                           '       BC not set\n'
                          f'       target node external label and DoF: {nodeLab}, {locDof}')
        else:
            print( 'ERROR: trying to set a BC on a node that is not in the structure\n'
                   '       BC not set\n'
                  f'       target node external label: {nodeLab}')
    
    def eqnums(self):
        # determine an eq. label for each node's DoF
        for node in self.nodes:
            # assign eq. numbers to node's DoFs
            node.dofLab=np.arange( self.ndofn, self.ndofn+3 )
            # constructu dict. for labeling translation
            for idx, eqLab in enumerate(node.dofLab):
                self.eqLab[eqLab]=[node, idx]
            # update number of DoFs in problem
            self.ndofn = self.ndofn+3
            
        # initially assign eq. labels to elements
        # possibly updated at stiffness mtx determination
        # based on specific info of the element
        for elm in self.elmts:
            for node in elm.nodes:
                elm.dofLab=np.concatenate((elm.dofLab,node.dofLab),axis=0)
        
    def delete_unused_dof(self):
        '''
        assumes a symmetric stiffness matrix
        '''
        
        # determine rows and cols to delete
        # iterate over global sitffness matrix rows
        it = np.nditer(self.K_full.row, flags=['f_index'])
        ids2del=[]
        for row in it:
            if np.all(self.K_full.mtx[it.index, :] == 0):
                ids2del.append(it.index)
        # delete
        self.K_full.mtx = np.delete(self.K_full.mtx, ids2del, axis=0) # delete rows
        self.K_full.row = np.delete(self.K_full.row, ids2del, axis=0) # delete rows
        self.K_full.mtx = np.delete(self.K_full.mtx, ids2del, axis=1) # delete cols
        self.K_full.col = np.delete(self.K_full.col, ids2del, axis=0) # delete cols
        # update num of DoF
        self.ndofn=self.K_full.col.size
    
    def assemble_K(self):
        # determine number of eq. per node and their labels
        self.eqnums()
        
        # prepare full matrix for assembling
        mtx = np.zeros((self.ndofn, self.ndofn))
        dof = np.arange(0, self.ndofn)
        self.K_full.set_mtx_dof(mtx, row=dof, col=dof)
        
        # assembly
        for elm in self.elmts:
            # determine stiffness matrix
            elm.detStiff()
            # iterate over element's DoFs
            it = np.nditer(elm.dofLab, flags=['f_index'])
            # and assemble the element matrix one row at a time
            for row in it:
                self.K_full.mtx[row,elm.dofLab]=elm.K_glo[it.index,:]
                
        # delete rows and cols associated to DoFs with no stiffness
        self.delete_unused_dof()

    def impose_BCs(self):
        
        if not self.K_full.mtx.size > 0:
            print( 'ERROR: trying to set BCs to a structure whose stiffness matrix was not assembled\n'
                   '       BCs not set')
        else:
            # PREPARE ARRAYS
            self.BC_labs = numpy.empty(self.BC_dict.shape(), dtype=int)
            self.BC_vals = numpy.empty(self.BC_dict.shape(), dtype=float)
            for idx, (node, locDof_val) in enumerate(self.BC_dict.items()):
                self.BC_labs[idx] = node.dofLab[ locDof_val[0] ]
                self.BC_vals[idx] =              locDof_val[1]
            
            self.K_wBCs = copy.copy(self.K_full)    # shalow copy - is it necesary to make a deepcopy?
                
            # IMPOSE
            # find the indexes
            ids2del=[]
            idsProblem=[]
            for targetLab in self.BC_labs:
                if targetLab in self.K_wBCs.row:
                    idx=np.where(self.K_wBCs.row == targetLab)[0]
                    ids2del.append(idx)
                else:
                    print( 'WARNING: trying to set a BC to a DoF without associated stiffness\n'
                           '         the row and col for the DoF were eliminated during assembly\n'
                           '         BC not set\n'
                          f'         equation internal label: {targetLab}\n'
                          f'         node external label: {self.eqLab[targetLab][0].extLab}\n'
                          f'         local DoF label: {self.eqLab[targetLab][1]}')
                    idx=np.where(self.BC_labs == targetLab)[0]
                    idsProblem.append(idx)
            # delete BCs over DoFs without stiffness
            if len(idsProblem)!=0:
                self.BC_labs = np.delete(self.BC_labs, idsProblem, axis=0)
                self.BC_vals = np.delete(self.BC_vals, idsProblem, axis=0)
            
            # delete rows
            self.K_wBCs.mtx = np.delete(self.K_wBCs.mtx, ids2del, axis=0) # delete rows
            self.K_wBCs.row = np.delete(self.K_wBCs.row, ids2del, axis=0) # delete rows
            # save cols for loading
            self.K_BCLo.mtx = copy.copy( self.K_wBCs.mtx[:,self.BC_labs] )
            self.K_BCLo.row = copy.copy( self.K_wBCs.row )
            self.K_BCLo.col = copy.copy( self.BC_labs )
            # delete cols
            self.K_wBCs.mtx = np.delete(self.K_wBCs.mtx, ids2del, axis=1) # delete cols
            self.K_wBCs.col = np.delete(self.K_wBCs.col, ids2del, axis=0) # delete cols

    def resolver(self):
        # Implementación del solucionador del sistema de ecuaciones
        pass


# RUN
if __name__ == '__main__':
    
    # simple structure
    if True:
        # create structure object
        str=Struc2D3Dof('simple frame')
        
        # create node and elem objects and add them to the structure
        if True: # in one step - recomended
            str.create_node(10, np.array([0,0  ], dtype=float))
            str.create_node(20, np.array([0,5.5], dtype=float))
            str.create_node(30, np.array([3,0  ], dtype=float))
            str.create_node(40, np.array([3,5.5], dtype=float))
            
            str.create_elm(1, [10,20], 210e3, 100, 10**4/12)
            str.create_elm(2, [10,40], 210e3, 100, 10**4/12)
            str.create_elm(3, [20,40], 210e3, 100, 10**4/12)
            
        else: # in two steps - not recomended
            n1=Node(10, np.array([0,0  ], dtype=float))
            n2=Node(20, np.array([0,5.5], dtype=float))
            n3=Node(30, np.array([3,0  ], dtype=float))
            n4=Node(40, np.array([3,5.5], dtype=float))
            str.add_node(n1)
            str.add_node(n2)
            str.add_node(n3)
            str.add_node(n4)
        
            elm1=ElmFrame2D(1, [n1,n2], 210e3, 100, 10**4/12)
            elm2=ElmFrame2D(2, [n1,n4], 210e3, 100, 10**4/12)
            elm3=ElmFrame2D(3, [n2,n4], 210e3, 100, 10**4/12, joints='both')
            str.add_elm(elm1)
            str.add_elm(elm2)
            str.add_elm(elm3)
        
        # assemble global stiffness matrix
        str.assemble_K()
    
    if False:
        pass
    
    
    