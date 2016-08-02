import numpy as np

def read_bin(name, ndims = 2):
    '''
        reading from a binary file in 3 dimensions
    '''
    
    with open(name, mode="rb") as file:
        file_raw = file.read()
        
    if ndims == 2:
        n,p = np.fromstring(file_raw[:4],dtype=np.uint16)
        file_trans = np.fromstring(file_raw[4:], dtype=np.uint8)
        f = np.reshape(file_trans, (n,p), order = "F")
        
    else: 
        n,p,q = np.fromstring(file_raw[:6],dtype=np.uint16)
        file_trans = np.fromstring(file_raw[6:], dtype=np.uint8)
        f = np.reshape(file_trans, (n,p,q), order = "F")
    
    return f
