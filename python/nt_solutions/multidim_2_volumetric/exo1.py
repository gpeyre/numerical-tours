MW = np.copy(M)
for j in range(1,int(np.log2(n))+1):
    p = int(n/2**(j-1))
    sel = np.arange(0,p)
    even = np.arange(0,p,2)
    odd = np.arange(1,p,2)
    # average/ difference along X
    MW[np.ix_(sel,sel,sel)] = np.concatenate(((MW[np.ix_(even,sel,sel)] + MW[np.ix_(odd,sel,sel)])/np.sqrt(2), (MW[np.ix_(even,sel,sel)] - MW[np.ix_(odd,sel,sel)])/np.sqrt(2)),0)
    # average/ difference along Y
    MW[np.ix_(sel,sel,sel)] = np.concatenate(((MW[np.ix_(sel,even,sel)] + MW[np.ix_(sel,odd,sel)])/np.sqrt(2), (MW[np.ix_(sel,even,sel)] - MW[np.ix_(sel,odd,sel)])/np.sqrt(2)),1)
    # average/ difference along Z
    MW[np.ix_(sel,sel,sel)] = np.concatenate(((MW[np.ix_(sel,sel,even)] + MW[np.ix_(sel,sel,odd)])/np.sqrt(2), (MW[np.ix_(sel,sel,even)] - MW[np.ix_(sel,sel,odd)])/np.sqrt(2)),2)