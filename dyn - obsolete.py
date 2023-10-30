

def solreg (P,K,w,Om,zita,t):
    
    """
    Función para la determinación de la Respuesta en Régimen de un Oscilador 
    Simple excitado por una carga armónica.
    
    """

    import numpy as np

    beta  = Om/w # Relación entre Frecuencia de la carga y Frecuencia de la Estructura
    gamma = 1/np.sqrt((1-beta**2)**2+(2*zita*beta)**2) # Factor de Amplificación Dinámico
    rho   = P*gamma/K # Amplitud de la Respuesta
    
    # Determinación del ángulo de Desfasaje de la Respuesta
    if zita==0 and beta>1:
        theta=np.pi
    else:
        theta = np.arccos((1-beta**2)*gamma)
    
    # Respuesta Dinámica    
    U     = rho*np.sin(Om*t-theta)
    return U

def statConds(K,DynDoFs):
    
    """
    condensación estática de la matriz de rigidez
    
    input ---------------------------------------------------------------------
    K:       numpy array - matriz de regidez con condiciones de vínculo impuestas
    DynDoFs: list - grados de libertad dinámicos
             nombrados según las filas y columnas correspondientes de K
           
    uotput --------------------------------------------------------------------
    Kc: numpy array - matriz de rigidez condensada
    """
    
    import numpy as np
    
    
    # determine geometrical degrees of freedom
    for i in range(K.shape[0]):
        try:
            b=DynDoFs.index(i)
        except ValueError:
            if 'GeoDoFs' in globals() or 'GeoDoFs' in locals():
                GeoDoFs.append(i)
            else:
                GeoDoFs = [i]
        else:
            "Do nothing"
    
    DynDoFs = np.asarray(DynDoFs, dtype=int)
    GeoDoFs = np.asarray(GeoDoFs, dtype=int)
    
    
    row = np.transpose(np.kron(np.ones((DynDoFs.size,1),dtype=int),DynDoFs))
    col =              np.kron(np.ones((DynDoFs.size,1),dtype=int),DynDoFs)
    K11 = K[row, col]
    
    row = np.transpose(np.kron(np.ones((GeoDoFs.size,1),dtype=int),DynDoFs))
    col =              np.kron(np.ones((DynDoFs.size,1),dtype=int),GeoDoFs)
    K12 = K[row, col]
    
    row = np.transpose(np.kron(np.ones((GeoDoFs.size,1),dtype=int),GeoDoFs))
    col =              np.kron(np.ones((GeoDoFs.size,1),dtype=int),GeoDoFs)
    K22 = K[row, col]
    
    Kc=K11-np.matmul(np.matmul(K12,np.linalg.inv(K22)),np.transpose(K12))
    
    return Kc
