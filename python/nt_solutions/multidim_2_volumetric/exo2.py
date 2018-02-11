M1 = np.copy(MWT)

for j in range(int(np.log2(n)),0,-1):
    p = int(n/2**j)
    sel = np.arange(0,p)
    sel1 = np.arange(0,2*p)
    selw = np.arange(p,2*p)
    even = np.arange(0,2*p,2)
    odd = np.arange(1,2*p,2)
    # average/ difference along X
    A = M1[np.ix_(sel,sel1,sel1)]
    D = M1[np.ix_(selw,sel1,sel1)]
    M1[np.ix_(even,sel1,sel1)] = (A + D)/np.sqrt(2)
    M1[np.ix_(odd,sel1,sel1)] = (A - D)/np.sqrt(2)
    # average/ difference along Y
    A = M1[np.ix_(sel1,sel,sel1)]
    D = M1[np.ix_(sel1,selw,sel1)]
    M1[np.ix_(sel1,even,sel1)] = (A + D)/np.sqrt(2)
    M1[np.ix_(sel1,odd,sel1)] = (A - D)/np.sqrt(2)
    # average/ difference along Z
    A = M1[np.ix_(sel1,sel1,sel)]
    D = M1[np.ix_(sel1,sel1,selw)]
    M1[np.ix_(sel1,sel1,even)] = (A + D)/np.sqrt(2)
    M1[np.ix_(sel1,sel1,odd)] = (A - D)/np.sqrt(2)