
"""
Direct Stifness Method
Mauro S. Maza - FCEFyN (UNC) - 2024-06-23
"""
libName='DSM >>Direct Stifness Method<<'
libAuth='Mauro S. Maza'
libAffi='Facultad de Ciencias Exactas, Físicas y Naturales | Universidad Nacional de Córdoba'
libVers='vX.Y.Z-YYYYMMDD'


"""
COSAS POR IMPLEMENTAR
- ElmFrame2D: incorporar la opción para articular extremos
- Node.add_dofs: asegurarse que el arreglo quede como vector columna    
        # IMPLEMENTAR MÉTODO PARA AGREGAR MUCHOS NODOS Y MUCHOS ELEMENTOS A LA VEZ
        # IMPLEMENTAR MÉTODO PARA ACTUALIZAR COORDENADAS DE UN PUNTO
        # PROBAR SI AL CAMBIAR LAS COORDENADA DE UN NODO, SE ACTUALIZA LA REFERENCIA EN EL ELEMENTO
        
        # se puede mejorar la programación creando una función que ensamble cualquier cosa y usando eso múltiples veces
- on BC input, check if BC already exists for that node and dof
"""

# IMPORTING ZONE
import numpy as np
from typing import List
import copy
import pickle
from datetime import datetime
# import matplotlib.pyplot as plt


# PRIMITIVES
def matrix2FormattedTable( matrix: np.array([]), fmt: list ):
    """
    Accepts a NumPy array (ndim=1,2) and returns a formatted table in a string variable.
    The numbers in each column are formatted according to de fmt argument.
    fmt: list of list, each element being of the form [repetitions, format_string]
    The sum of "repetitions" must be equal to the number of columns in matrix.
    No spaces are added between elemnts by default.
    "matrix" might be of dtype=object in order to accept "None" elements.
    """
    # prepare formatting strings
    repSum=0
    fs=[] # expanded list of formatting strings
    for subList in fmt:
        repSum+=subList[0]
        for i in range(subList[0]):
            fs.append(subList[1])
    
    # translate
    if matrix.ndim==1  and  repSum==matrix.shape[0]:
        fT='' # formatted table
        for i in range(matrix.shape[0]):
            fi=fs[i]
            fT+=f'{ matrix[i] :{fi}}'
        return fT
    
    elif matrix.ndim==2  and  repSum==matrix.shape[1]:
        fT='' # formatted table
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                fT+=f'{ matrix[i,j]:{fs[j]}}'
            else:
                fT+='\n'
        return fT
        
    else:
        print('ERROR in matrix2FormattedTable input arguments')


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

    def add_dofs(self, dofLab: np.array([], dtype=int)):
        self.dofLab = dofLab
        self.dofVal = np.zeros( self.dofLab.shape )

    def set_dofVal(self, dofVal: np.array([], dtype=float)):
        self.dofVal = dofVal


class dpv:
    """
    dynamic property vector
    has a vector and
    another vector indicating the numbers of equations to which the rows refer to
    """
    
    def __init__(self):
        self.vtr = np.array([], dtype=float)
        self.row = np.array([], dtype=int)      # labels of global DoFs associated to each row of matrix
        
    def set_vtr_dof(self, vtr: np.array([], dtype=float),
                          row: np.array([], dtype=int  )):
        if vtr.ndim != 1:
            print('ERROR: trying to set a dynamic property vector with more or less than 1 dimensions.')
        else:
            if row.ndim!=1 or row.shape[0]!=vtr.shape[0]:
                print('ERROR: trying to set a dynamic property vector with a wrong DoF vector.')
            else:
                self.vtr = vtr
                self.row = row


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
        self.R      = np.array([], dtype=float) # rotation matrix: U_glo = R * U_loc
        self.K_loc  = np.array([], dtype=float) # stiffness matrix in local  coordinates
        self.K_glo  = np.array([], dtype=float) # stiffness matrix in global coordinates
        self.dofLab = np.array([], dtype=int)   # global DoFs' labels, assigned on assembly
        self.U_loc  = np.array([], dtype=float) # displacements of elem ends in local coordinates
        self.U_glo  = np.array([], dtype=float) # displacements of elem ends in global coordinates
        self.L_loc  = np.array([], dtype=float) # external loads on elem ends in local coordinates
        self.L_glo  = np.array([], dtype=float) # external loads on elem ends in global coordinates
    
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
            if not local:
                self.dofLab = np.block( [ self.nodes[0].dofLab[:] , self.nodes[1].dofLab[:] ] )
            
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
        self.K_loc = self.stiffM(local=True )     # stiffness matrix in local  coordinates
        self.K_glo = self.stiffM(local=False)     # stiffness matrix in global coordinates
    
    def detR(self):
        """
        Rotation matrix for a 6-element vector
            > U_glo = R * U_loc
            > L_glo = R * L_loc
        """
        
        g1=self.t[0]
        g2=self.t[1]
        r = np.array( [ [ g1, -g2, 0 ] ,
                        [ g2,  g1, 0 ] ,
                        [  0,   0, 1 ] ])
        o = np.zeros_like(r)
        self.R = np.block( [ [ r, o ] ,
                             [ o, r ] ])
        
    def elmSolve(self):
        # global coordinates
        self.U_glo = np.block( [ self.nodes[0].dofVal[:] , self.nodes[1].dofVal[:] ] )
        self.L_glo = np.matmul( self.K_glo, self.U_glo )
        
        # local coords.
        self.detR()
        self.U_loc = np.matmul( np.transpose(self.R), self.U_glo )
        self.L_loc = np.matmul( np.transpose(self.R), self.L_glo )
        
    def getL( self, nodeLocLab: int, locDof: int, local=True ):
        """
        retrieve elemental nodal load (in local or global coordinates)
        nodeLocLab runs from 0 to nnode-1
        locDof runs from 0 to nDof-1
        """
        idx=3*nodeLocLab+locDof
        if local:
            L=self.L_loc[idx]
        else:
            L=self.L_glo[idx]
        return L
        
    def reportBasic(self, titles=True, retType='string'):
        """
        Report element basic data
        titles: if True and retType='string', the first line of "tab" has the labels of the info in the second line
        retType: return type
                 might be "string" or "matrix"
                 if "matrix", "joint" info is not returned
        """
        mat1=np.block([ self.extLab, self.nodes[0].extLab, self.nodes[1].extLab ])
        mat2=np.block([ self.E, self.A, self.I, self.intLab, self.L, self.t ])
        if retType=='string':
            #      extLab    nodes_extLabs
            fmt1=[ [1,'8g'], [2,'8g'] ]
            tab1=matrix2FormattedTable( mat1, fmt1 )
            #      sec_props    intLab    L            t
            fmt2=[ [3,'13.5e'], [1,'8g'], [1,'13.5e'], [2,'9.4f'] ]
            tab2=matrix2FormattedTable( mat2, fmt2 )
            
            if titles:
                # the following must be in accordance with fmt1 and fmt2
                tab='  extLab  node_1  node_2  joints  E            A            I            intLab  L            gamma1   gamma2\n'
            else:
                tab=''
                
            tab += tab1 + f'{self.joints:>8}' + tab2
            return tab
        else:
            mat=np.block([mat1,mat2])
            return mat


class Struc2D3Dof:
    """
    General bidimensional structure with 3 DoF per node
    under static load system
    """
    
    def __init__(self, name='New Project'):
        # MAIN ATTRIBUTES
        self.name  = name
        self.desc  = ''
        self.nodes = []         # list of node objects
        self.elmts = []         # list of elem objects
        self.eqLab = {}         # translation from global DoF (eqLab) to node and local DoF {eqLab:[node_object, local_DoF_0->2]}
                                # the inverse is done with method self.getEqLab()
        self.ndofn = int(0)
        self.spDof = []         # suppressed DoFs without stiffness associated
        
        # stiffness matrices
        self.K_full = dpm()     # assembled global matrix (without BCs)
        self.K_wBCs = dpm()     # assembled global matrix with BCs (rows and cols deleted)
        self.K_BCLo = dpm()     # columns of (row reduced-)K_full associated to (zero and non-zero) BCs (added to load vector)
        self.K_cc   = dpm()
        self.K_hc   = dpm()
        self.K_hh   = dpm()
        self.K_cond = dpm()
        
        # load vectors
        self.reacts = dpv()     # rel. to K_full - reactions: K_full * U_full = L_full + reacts
        self.L_full = dpv()     # rel. to K_full - assembled global vector (without BCs)
        self.L_wBCs = dpv()     # rel. to K_wBCs - rows deleted and K_BCLo*BC_vals added
        
        # displacements vectors
        self.U_full = dpv()     # rel. to K_full
        self.U_wBCs = dpv()     # rel. to K_wBCs
        
        # OTHER ATTRIBUTES
        # nodal loads (aligned to global coordinate system)
        self.NL_dict = {}                           # dictionary for NLs input {internal_label:[node_object local_Dof_label, value]}
        
        # BCs aligned to global coordinate system
        self.BC_dict = {}                           # dictionary for BCs input {internal_label:[node_object local_Dof_label, value]}
        self.BC_labs = np.array([], dtype=int)      # labels of global DoFs with imposed value
        self.BC_vals = np.array([], dtype=float)    # values imposed

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
    
    def findNodes( self, node_labs: List[int] ):
        """
        find node objects from external label
        """
        
        # check that there are no repetitions in the list
        if len(node_labs) != len(set(node_labs)) :
            print( 'WARNING: node external labels repeated in list for findNodes()\n'
                  f'         nodes external labels: {node_labs}')
        
        nodeObjectsList=[]
        for targetLab in node_labs:
            # check if the node exists - try-except clause???
            if any(   targetLab==possibleNode.extLab   for possibleNode in self.nodes):
                # find the node
                for possibleNode in self.nodes:
                    if targetLab==possibleNode.extLab:
                        nodeObjectsList.append(possibleNode)
                        break
            else:
                print( 'ERROR: trying to find a node that is not in the structure\n'
                      f'       node external label: {targetLab}')
        
        return nodeObjectsList
    
    def add_elm(self, newElm):
        if isinstance( newElm, ElmFrame2D):
            # check that the elem does not exist yet
            if any(newElm is oldElm for oldElm in self.elmts):
                print( 'ERROR: trying to add an elem that already exists in the structure\n'
                       '       elem not added\n'
                      f'       elem external label: {newElm.extLab}')
            else:
                # check that the external label of the elem does not exist yet
                if newElm.extLab in [oldElm.extLab for oldElm in self.elmts]:
                    print( 'ERROR: trying to add an elem with an existing external label\n'
                           '       elem not added\n'
                          f'       elem external label: {newElm.extLab}')
                else:
                    # check that the nodes defining the elem exist in the structure
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

    def create_elm(self, extLab: int, node_labs: List[int], E: float, A: float, I: float, joints='none'):
        """
        arg "node_labs" here is a list of external labels of nodes
        which is different from arg "nodes" in ElmFrame2D.__init__() method
        """
        
        # create the list of node objects for ElmFrame2D.__init__()
        nodeObjectsList = self.findNodes(node_labs)
        # check that all the nodes defining the elem exist in the structure
        if len(node_labs)==len(nodeObjectsList):
            elm=ElmFrame2D(extLab, nodeObjectsList, E, A, I, joints)
            self.add_elm(elm)
        else:
            print( 'ERROR: trying to add an elem whose defining nodes are not in the structure\n'
                   '       elem not added\n'
                  f'       elem external label: {elm.extLab}')
                  
    def findElmts( self, elms_labs: List[int] ):
        """
        find element objects from external label
        """
        
        # check that there are no repetitions in the list
        if len(elms_labs) != len(set(elms_labs)) :
            print( 'WARNING: elems external labels repeated in list for findElmts()\n'
                  f'         elems external labels: {elms_labs}')
        
        elemObjectsList=[]
        for targetLab in elms_labs:
            # check if the node exists - try-except clause???
            if any(   targetLab==possibleElem.extLab   for possibleElem in self.elmts):
                # find the node
                for possibleElem in self.elmts:
                    if targetLab==possibleElem.extLab:
                        elemObjectsList.append(possibleElem)
                        break
            else:
                print( 'ERROR: trying to find an elem that is not in the structure\n'
                      f'       elem external label: {targetLab}')
        
        return elemObjectsList
    
    def create_NL(self, nodeLab: int, locDof: int, val: float):
        """
        nodeLab: external label of node
        locDof: local DoF label
                runs from 1 to nDof (engineering labeling)
        val: load value
        """
        
        # find the node
        nodeObjectsList = self.findNodes([nodeLab])
        # check that the node does exist
        if len(nodeObjectsList)==1:
            node=nodeObjectsList[0]
            locDof -= 1
            # check if the DoF is correct
            # try-except clause should check if locDof exists and inform
            if locDof in range(0,3):
                self.NL_dict[len(self.NL_dict)] = [node, locDof, val]   # dictionary for NLs input {internal_label:[node_object local_Dof_label, value]}
                if self.K_wBCs.mtx.size > 0:
                    print('WARNING: trying to add a nodal load to a structure\n'
                          '         whose BCs has already been imposed.\n'
                          '         RECALC')
            else:
                print( 'ERROR: trying to set a nodal load on a non-existent local Dof\n'
                       '       nodal load not set\n'
                      f'       target node external label and DoF: {nodeLab}, {locDof}')
        else:
            print( 'ERROR: trying to set a nodal load on a node that is not in the structure\n'
                   '       nodal load not set\n'
                  f'       target node external label: {nodeLab}')
    
    def create_BC(self, nodeLab: int, locDof: int, val: float):
        """
        nodeLab: external label of node
        locDof: local DoF label
                runs from 1 to nDof (engineering labeling)
        val: value imposed        
        """
        
        # find the node
        nodeObjectsList = self.findNodes([nodeLab])
        # check that the node does exist
        if len(nodeObjectsList)==1:
            node=nodeObjectsList[0]
            locDof -= 1
            # check if the DoF is correct
            # try-except clause should check if locDof exists and inform
            if locDof in range(0,3):
                self.BC_dict[len(self.BC_dict)] = [node, locDof, val]   # dictionary for BCs input {internal_label:[node_object local_Dof_label, value]}
                if self.K_wBCs.mtx.size > 0:
                    print('WARNING: trying to add a BC to a structure\n'
                          '         whose BCs has already been imposed.\n'
                          '         RECALC')
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
        Find global DoFs without associated stiffness
        and delete corresponding rows and cols in "_full"
        arrays of the structure and in BCs.
        As the "_wBCs" arrays start as copies of the "_full" arrays,
        they inherit the absence of this "suppressed DoFs".
        No change is made to the NL_dict attribute.
        No change is made to the elemental matrices or vectors.
        Elemental displacement vectors are completed with zeros.
        
        This method assumes a symmetric assembled stiffness matrix.
        '''
        
        # DETERMINE ROWS AND COLS TO DELETE
        # iterate over global sitffness matrix rows
        it = np.nditer(self.K_full.row, flags=['f_index'])
        ids2del=[]
        for row in it:
            if np.all(self.K_full.mtx[it.index, :] == 0):
                ids2del.append(it.index)
        self.spDof=self.K_full.row[ids2del].tolist()
        
        # DELETE
        # stiffness matrix
        self.K_full.mtx = np.delete(self.K_full.mtx, ids2del, axis=0) # delete rows
        self.K_full.row = np.delete(self.K_full.row, ids2del, axis=0) # delete rows
        self.K_full.mtx = np.delete(self.K_full.mtx, ids2del, axis=1) # delete cols
        self.K_full.col = np.delete(self.K_full.col, ids2del, axis=0) # delete cols
        # load vector
        self.L_full.vtr = np.delete(self.L_full.vtr, ids2del, axis=0) # delete rows
        self.L_full.row = np.delete(self.L_full.row, ids2del, axis=0) # delete rows
        # disp. vector
        self.U_full.vtr = np.delete(self.U_full.vtr, ids2del, axis=0) # delete rows
        self.U_full.row = np.delete(self.U_full.row, ids2del, axis=0) # delete rows
        # boundary conditions
        #   find BCs to delete
        keys2del=[]
        for bcLab, node_locDof_val in self.BC_dict.items():
            node=node_locDof_val[0]
            locDof=node_locDof_val[1]
            if node.dofLab[ locDof ] in self.spDof:
                keys2del.append(bcLab)
        #   delete
        for key in keys2del:
            del self.BC_dict[key]
        
        # UPDATE NUM OF DOF
        self.ndofn=self.K_full.col.size
    
    def assemble(self):
        # determine number of eq. per node and their labels
        self.eqnums()
        
        # INITIALIZE
        dof = np.arange(0, self.ndofn)
        # stiffness matrix
        mtx = np.zeros((self.ndofn, self.ndofn))
        self.K_full.set_mtx_dof(mtx, row=dof, col=dof)
        # load vector
        vtr = np.zeros((self.ndofn))
        self.L_full.set_vtr_dof(vtr, row=dof)
        self.reacts.set_vtr_dof(vtr, row=dof)
        # disp. vector
        self.U_full.set_vtr_dof(vtr, row=dof)
        
        # ASSEMBLY
        # stiffness matrix
        for elm in self.elmts:
            # determine stiffness matrix
            elm.detStiff()
            # iterate over element's DoFs
            it = np.nditer(elm.dofLab, flags=['f_index'])
            # and assemble the element matrix one row at a time
            for row in it:
                self.K_full.mtx[row,elm.dofLab] += elm.K_glo[it.index,:]
        # load vector
        for node_locDof_val in self.NL_dict.values():
            node   = node_locDof_val[0]
            locDof = node_locDof_val[1]
            val    = node_locDof_val[2]
            idx = node.dofLab[ locDof ]
            self.L_full.vtr[idx] += val
        
        # DELETE rows and cols associated to DoFs with no stiffness
        self.delete_unused_dof()

    def impose_BCs(self):
        
        if not self.K_full.mtx.size > 0:
            print( 'ERROR: trying to set BCs to a structure whose stiffness matrix was not assembled\n'
                   '       BCs not set')
        else:
            # PREPARE ARRAYS
            self.BC_labs = np.empty(len(self.BC_dict), dtype=int)
            self.BC_vals = np.empty(len(self.BC_dict), dtype=float)
            for idx, (bcLab, node_locDof_val) in enumerate(self.BC_dict.items()):
                node  =node_locDof_val[0]
                locDof=node_locDof_val[1]
                val   =node_locDof_val[2]
                self.BC_labs[idx] = node.dofLab[ locDof ]
                self.BC_vals[idx] = val
            
            self.K_wBCs = copy.deepcopy(self.K_full)
            self.L_wBCs = copy.deepcopy(self.L_full)
            self.U_wBCs = copy.deepcopy(self.U_full)
                
            # IMPOSE
            # find the indexes
            ids2del=[]
            for targetLab in self.BC_labs:
                idx=np.where(self.K_wBCs.row == targetLab)[0][0]
                ids2del.append(idx)
            
            # stiffness matrix
            # delete rows
            self.K_wBCs.mtx = np.delete(self.K_wBCs.mtx, ids2del, axis=0) # delete rows
            self.K_wBCs.row = np.delete(self.K_wBCs.row, ids2del, axis=0) # delete rows
            # save cols for loading
            self.K_BCLo.mtx = copy.copy( self.K_wBCs.mtx[:,ids2del] )
            self.K_BCLo.row = copy.copy( self.K_wBCs.row )
            self.K_BCLo.col = copy.copy( self.BC_labs )
            # delete cols
            self.K_wBCs.mtx = np.delete(self.K_wBCs.mtx, ids2del, axis=1) # delete cols
            self.K_wBCs.col = np.delete(self.K_wBCs.col, ids2del, axis=0) # delete cols
            
            # load vector
            self.L_wBCs.vtr = np.delete(self.L_wBCs.vtr, ids2del, axis=0) # delete rows
            self.L_wBCs.row = np.delete(self.L_wBCs.row, ids2del, axis=0) # delete rows
            self.L_wBCs.vtr += np.matmul( self.K_BCLo.mtx, self.BC_vals)
            
            # disp. vector
            self.U_wBCs.vtr = np.delete(self.U_wBCs.vtr, ids2del, axis=0) # delete rows
            self.U_wBCs.row = np.delete(self.U_wBCs.row, ids2del, axis=0) # delete rows
            self.U_full.vtr[ids2del]=self.BC_vals[:] # this works because ids2del is determined here, some rows above

    def solve(self):
        # UNKNOWN DISP.
        # assembled displacement vector
        self.U_wBCs.vtr = np.linalg.solve(self.K_wBCs.mtx, self.L_wBCs.vtr)
        ids = [np.where(self.U_full.row == dofLab)[0][0] for dofLab in self.U_wBCs.row]
        self.U_full.vtr[ids]=self.U_wBCs.vtr[:]
        
        # REACTIONS
        self.reacts.vtr = np.matmul( self.K_full.mtx, self.U_full.vtr) - self.L_full.vtr
        
        # NODES
        for node in self.nodes:
            node.dofVal=np.zeros( node.dofLab.shape, dtype=float)
            for locIdx, dofLab in enumerate(node.dofLab):
                if dofLab in self.U_full.row:
                    gloIdx = np.where(self.U_full.row == dofLab)[0][0]
                    node.dofVal[locIdx]=self.U_full.vtr[gloIdx]
        # ELEMENTS
        for elm in self.elmts:
            elm.elmSolve()
            
    def getEqLab( self, nodeExtLab: int, locDof: int ):
        """
        determine global Dof label (eqLab) for node external label and local DoF
        locDof runs from 0 to nDof-1
        """
        # check if the node exists
        if any(   nodeExtLab==possibleNode.extLab   for possibleNode in self.nodes):
            # find the node
            for possibleNode in self.nodes:
                if nodeExtLab==possibleNode.extLab:
                    # try-except clause should determine if locDof exists and inform
                    globalDof = possibleNode.dofLab[locDof]
                    return globalDof
                    break
        else:
            print( 'ERROR: trying to find a global DoF label for a node that is not in the structure\n'
                  f'       node external label: {nodeExtLab}')
                  
    def getU( self, nodeExtLab: int, locDof: int ):
        """
        retrieve global displacement
        locDof runs from 0 to nDof-1
        """
        globalDof = self.getEqLab( nodeExtLab, locDof )
        U = self.U_full.vtr[globalDof]
        return U
                  
    def getR( self, nodeExtLab: int, locDof: int ):
        """
        retrieve reaction
        locDof runs from 0 to nDof-1
        """
        globalDof = self.getEqLab( nodeExtLab, locDof )
        R = self.reacts.vtr[globalDof]
        return R
                  
    def ask4U( self, nodeExtLab: int, locDof: int ):
        """
        retrieve global displacement
        locDof runs from 1 to nDof (engineering labeling)
        """
        U = self.getU( nodeExtLab, locDof-1 )
        return U
                  
    def ask4R( self, nodeExtLab: int, locDof: int ):
        """
        retrieve reaction
        locDof runs from 1 to nDof (engineering labeling)
        """
        R = self.getR( nodeExtLab, locDof-1 )
        return R
                  
    def ask4elem( self, elemExtLab: int ):
        """
        retrieve element object with specific external label
        """
        elmList = self.findElmts( [elemExtLab] )
        return elmList[0]
                  
    def ask4ENL( self, elemExtLab: int, nodeLocLab: int, locDof: int, local=True ):
        """
        retrieve elemental nodal load (in local or global coordinates)
        nodeLocLab runs from 1 to nnode (engineering labeling)
        locDof runs from 1 to nDof (engineering labeling)
        """
        elmList = self.findElmts( [elemExtLab] )
        elm=elmList[0]
        ENL=elm.getL( nodeLocLab-1, locDof-1, local )
        return ENL
        
    def save2bin(self, fileName='', path=''):
        """
        fileName: string without extension
        path: string
              if '', file will be saved in local directory
              might be
                - a full path; example: 'E:\GitHub\cei'
                - or a relative path; example: 'subfolder\subsufolder'
              note that in both cases there is no initial or ending '\' character
              the directory will not be created (in case it does not exist)
        """
        if fileName=='':
            fileName=self.name
        if path!='':
            fileName=path+'\\'+fileName
        
        with open( fileName+'.pkl', 'wb') as binFile:
            pickle.dump(self, binFile)
            
    
    def report(self, fileName='', path=''):
        """
        fileName: string without extension
        path: string
              if '', file will be saved in local directory
              might be
                - a full path; example: 'E:\GitHub\cei'
                - or a relative path; example: 'subfolder\subsufolder'
              note that in both cases there is no initial or ending '\' character
              the directory will not be created (in case it does not exist)
        """
        def dict2Arr(d,nDof):
            """
            NOT WORKING PROPERLY
            NOT BEING USED
            d: dict with data as NL_dict =  {internal_label:[node_object local_Dof_label, value]}.
            nDof: number of DoFs that should be considered for the returned matrix/table.
            "None" is set to the elements corresponding to DoFs without load or BC set.
            The returned array cannot be prosseced by matrix2FormattedTable()
            as for the 2024-08-30 version (becuase of the "None" elements).
            """
            # find nodes external labels and number of columns
            labs=set() # empty set
            maxDof=0
            for l in d.values():
                labs.add( l[0].extLab )
                maxDof = max( maxDof, l[1] )
            
            # construct the array
            arr=np.zeros( ( len(labs), maxDof+2 ), dtype=object)
            arr.fill(None)
            labs=list(labs)
            labs.sort()
            i=-1
            for lab in labs: # one row for each node with data
                i+=1
                arr[i,0]=lab
                for l in d.values():
                    if lab==l[0].extLab:
                        dof=l[1]
                        val=l[2]
                        if arr[i,dof+1] is None:
                            arr[i,dof+1]=val
                        else:
                            arr[i,dof+1]+=val
            return arr
            
        if fileName=='':
            fileName=self.name
        if path!='':
            fileName=path+'\\'+fileName
            
        with open( fileName+'.rep', 'w') as repFile:
            
            # HEADLINE
            headline = '\n\nCreated with '+libName+' '+libVers+ \
                         '\nby '+libAuth+ \
                         '\nat '+libAffi+'\n\n\n'
            repFile.write(headline)
            
            # PROBLEM DESCRIPTION
            repFile.write('Problem name: '+self.name+'\n')
            DaT = datetime.now()
            DaT = DaT.strftime("%d/%m/%Y %H:%M:%S")
            repFile.write('Date: '+DaT+'\n')
            if self.desc!='':
                repFile.write('Problem description:\n    '+self.desc+'\n\n\n')
            else:
                repFile.write('Problem description: no description.\n\n\n')
                
            # PROBLEM DATA
            repFile.write('PROBLEM DATA\n')
            
            # nodes
            repFile.write( '> Nodes labels, coordinates, internal labels and global DoF (eq. numbers)\n')
            repFile.write(f'>> Number of nodes: {len(self.nodes)}\n')
            repFile.write( '>> Table\n')
            # prepare data
            n=self.nodes[0]
            matrix=np.block([ n.extLab, n.coords, n.intLab, n.dofLab ])
            for n in self.nodes[1:]:
                nMtx = np.block([ n.extLab, n.coords, n.intLab, n.dofLab ])
                matrix=np.block([[matrix],[nMtx]])
            # the next works only for 2 coord and 3 DoF nodes
            repFile.write('  extLab   coord_1        coord_2       intLab    DoF1    DoF2    DoF3\n')
            #      extLab     coords       intLab    dofLab
            fmt=[  [1,'8g'], [2,'15.5e'], [1,'8g'], [3,'8g']  ]
            table=matrix2FormattedTable( matrix, fmt )
            repFile.write(table)
            repFile.write('\n')
            
            # elements
            repFile.write( '> Element global data\n')
            repFile.write(f'>> Number of elements: {len(self.elmts)}\n')
            repFile.write( '>> Table\n')
            # write titles
            e=self.elmts[0]
            repFile.write( e.reportBasic()+'\n' )
            # write the rest
            for e in self.elmts[1:]:
                repFile.write( e.reportBasic(titles=False)+'\n' )
            repFile.write('\n')
            
            # loads
            repFile.write( '> Loads\n')
            repFile.write( '>> Loads input\n')
            repFile.write(f'>>> Number of loads: {len(self.NL_dict)}\n')
            repFile.write( '>>> Table\n')
            repFile.write('    node     DoF    load\n')
            #     node      DoF       value
            fmt=[ [1,'8g'], [1,'8g'], [1,'15.5e'] ]
            for l in self.NL_dict.values():
                mat=np.array( [ l[0].extLab, l[1], l[2] ] )
                line = matrix2FormattedTable( mat, fmt )
                repFile.write( line + '\n' )
            repFile.write('\n')
            repFile.write('>> Loads processed\n')
            repFile.write('(TO BE IMPLEMENTED)\n')
            # the next works only for 2 coord and 3 DoF nodes
            #     node      values
            # fmt=[ [1,'8g'], [3,'15.5e'] ]
            # mat=dict2Arr(self.NL_dict)
            # table=matrix2FormattedTable( mat, fmt )
            # repFile.write('    node   L_in_DoF_1     L_in_DoF_2     L_in_DoF_3\n')
            # repFile.write(table)
            repFile.write('\n')
            
            # BCs
            repFile.write( '> Boundary conditions\n')
            repFile.write( '>> BCs input\n')
            repFile.write(f'>>> Number of loads: {len(self.NL_dict)}\n')
            repFile.write( '>>> Table\n')
            repFile.write('    node     DoF    load\n')
            #     node      DoF       value
            fmt=[ [1,'8g'], [1,'8g'], [1,'15.5e'] ]
            for l in self.NL_dict.values():
                mat=np.array( [ l[0].extLab, l[1], l[2] ] )
                line = matrix2FormattedTable( mat, fmt )
                repFile.write( line + '\n' )
            repFile.write('\n')
            repFile.write('>> BCs processed\n')
            repFile.write('(TO BE IMPLEMENTED)\n')
            repFile.write('\n')
                
            # INTERNAL DATA
            repFile.write('INTERNAL DATA - (TO BE IMPLEMENTED)\n')
            repFile.write('> Elemental matrices (global coordinates - with associated global DoFs)\n')
            repFile.write('> Data relative to suppressed DoFs\n')
            repFile.write('> Assembled matrix (K_full) and external load vector (L_full) (with associated global DoFs)\n')
            repFile.write('> Assembled matrix (K_wBCs) and external load vector (L_wBCs) after BCs (with associated global DoFs)\n')
            repFile.write('> Loads related to BCs (K_BCLo*BC_vals=vector) (with associated global DoFs)\n')
            repFile.write('\n')
                
            # RESULTS
            repFile.write('RESULTS - (TO BE IMPLEMENTED)\n')
            repFile.write('> Displacements (for every node and DoF - table-formatted)\n')
            repFile.write('> Reactions (for every node and DoF - table-formatted)\n')
            repFile.write('> Elemental nodal loads (table-formatted)\n')
            repFile.write('\n')
            
            
# RUN
if __name__ == '__main__':
    
    """
    Example #1 for dsm.py
    RC Hibbeler - Análisis Estructural - 8ed - Ejemplo 16.1
    In SI units
    Origin of coord. sys. on the bottom-left corner
    """
    if True:
        # import numpy as np
        # from dsm import Struc2D3Dof

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
        str.create_node(extLab=1, coords=np.array([ 0*feet2mm, 20*feet2mm], dtype=float))
        str.create_node(extLab=2, coords=np.array([20*feet2mm, 20*feet2mm], dtype=float))
        str.create_node(extLab=3, coords=np.array([20*feet2mm,  0*feet2mm], dtype=float))

        # elms
        E=29e3*ksi2MPa
        A=10*in22mm2
        I=500*in42mm4
        str.create_elm(extLab=1, node_labs=[1,2], E=E, A=A, I=I, joints='none')
        str.create_elm(extLab=2, node_labs=[2,3], E=E, A=A, I=I, joints='none')

        # loads
        str.create_NL(nodeLab=2, locDof=1, val=5*kip2N)

        # BCs
        str.create_BC(nodeLab=1, locDof=2, val=0) # node 1 - pinned (vertical)
        str.create_BC(nodeLab=3, locDof=1, val=0) # node 3 - fixed
        str.create_BC(nodeLab=3, locDof=2, val=0) # node 3 - fixed
        str.create_BC(nodeLab=3, locDof=3, val=0) # node 3 - fixed

        # solve
        str.assemble()
        str.impose_BCs()
        str.solve()

        # print
        print( 'Unknwon displacements:')
        print(f'> D1: node 2, u_x = {str.ask4U(2,1)*mm2in:9.6f} in')
        print(f'> D2: node 2, u_y = {str.ask4U(2,2)*mm2in:9.6f} in')
        print(f'> D3: node 2, phi = {str.ask4U(2,3)      :9.6f} rad')
        print(f'> D4: node 1, u_x = {str.ask4U(1,1)*mm2in:9.6f} in')
        print(f'> D5: node 1, phi = {str.ask4U(1,3)      :9.6f} rad')
        print( '')
        print( 'Reactions:')
        print(f'> Q6: node 1, R_y = {str.ask4R(1,2)/kip2N      :8.2f} kip')
        print(f'> Q7: node 3, R_x = {str.ask4R(3,1)/kip2N      :8.2f} kip')
        print(f'> Q8: node 3, R_y = {str.ask4R(3,2)/kip2N      :8.2f} kip')
        print(f'> Q9: node 3, Mom = {str.ask4R(3,3)/kip2N*mm2in:8.2f} kip*in')
        print( '')
        print( 'Elemental nodal loads in global coordinates:')
        print( 'element #1')
        print(f'> elem 1, node i, L_xg = {str.ask4ENL( elemExtLab=1, nodeLocLab=1, locDof=1, local=False )/kip2N      :8.2f} kip')
        print(f'> elem 1, node i, L_yg = {str.ask4ENL( elemExtLab=1, nodeLocLab=1, locDof=2, local=False )/kip2N      :8.2f} kip')
        print(f'> elem 1, node i, mome = {str.ask4ENL( elemExtLab=1, nodeLocLab=1, locDof=3, local=False )/kip2N*mm2in:8.2f} kip*in')
        print(f'> elem 1, node j, L_xg = {str.ask4ENL( elemExtLab=1, nodeLocLab=2, locDof=1, local=False )/kip2N      :8.2f} kip')
        print(f'> elem 1, node j, L_yg = {str.ask4ENL( elemExtLab=1, nodeLocLab=2, locDof=2, local=False )/kip2N      :8.2f} kip')
        print(f'> elem 1, node j, mome = {str.ask4ENL( elemExtLab=1, nodeLocLab=2, locDof=3, local=False )/kip2N*mm2in:8.2f} kip*in')
        print( 'element #2')
        print(f'> elem 2, node i, L_xg = {str.ask4ENL( elemExtLab=2, nodeLocLab=1, locDof=1, local=False )/kip2N      :8.2f} kip')
        print(f'> elem 2, node i, L_yg = {str.ask4ENL( elemExtLab=2, nodeLocLab=1, locDof=2, local=False )/kip2N      :8.2f} kip')
        print(f'> elem 2, node i, mome = {str.ask4ENL( elemExtLab=2, nodeLocLab=1, locDof=3, local=False )/kip2N*mm2in:8.2f} kip*in')
        print(f'> elem 2, node j, L_xg = {str.ask4ENL( elemExtLab=2, nodeLocLab=2, locDof=1, local=False )/kip2N      :8.2f} kip')
        print(f'> elem 2, node j, L_yg = {str.ask4ENL( elemExtLab=2, nodeLocLab=2, locDof=2, local=False )/kip2N      :8.2f} kip')
        print(f'> elem 2, node j, mome = {str.ask4ENL( elemExtLab=2, nodeLocLab=2, locDof=3, local=False )/kip2N*mm2in:8.2f} kip*in')
        
        # save the object as binary file
        str.save2bin() # using the object's "name" attribute for the file name
        str.save2bin( fileName='alternative name' )
        
        # read the objetc from binary file
        with open('alternative name.pkl', 'rb') as binFile:
            readObject = pickle.load(binFile)
            if not str is readObject:
                print('\n\nNow you have a deepcopy of the object loaded in memory\n\n')
        
        # print a report
        str.report()
    
    # simple structure
    if False:
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
    
    # try matrix2FormattedTable()
    if False:
        # array
        integers=np.array([[1, 2, 3],[4, 5, 6]], dtype=int)
        floatings=np.array([[1.1, 2.2, 3.33],[4, 5.888, 6.6]], dtype=float)
        matrix=np.block([integers,floatings]) # integers are converted to floats
        # format
        fmt=[  [1,'g'], [2,'4g'], [3,'10.2e']  ]
        table=matrix2FormattedTable( matrix, fmt )
        print(table)
    
    
    