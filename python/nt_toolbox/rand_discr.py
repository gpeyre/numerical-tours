import numpy as np
from numpy import random

def rand_discr(p, m = 1):
    """
        rand_discr - discrete random generator 
        
          y = rand_discr(p, n);
        
          y is a random vector of length n drawn from 
          a variable X such that
              p(i) = Prob( X=i )
        
          Copyright (c) 2004 Gabriel Peyré
    """

    # makes sure it sums to 1
    p = p/np.sum(p)
    
    n = len(p)
    coin = random.rand(m)
    cumprob = np.append(0,+ np.cumsum(p))
    sample = np.zeros(m)
    
    for j in range(n):
        ind = [(coin > cumprob[j]) & (coin <= cumprob[j+1])]
        sample[ind] = j
        
    return sample