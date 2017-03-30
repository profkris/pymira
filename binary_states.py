# -*- coding: utf-8 -*-
"""
Created on Tue Mar 07 09:04:50 2017

@author: simon
"""

import numpy as np
vals = np.asarray([14.7,-14.7,10.,-10.,-20.54,20.54])
nEl = tmp.shape[0]

mnState = None
mn = 1e6

from itertools import product
for i,state in enumerate(product([-1,1], repeat=nEl)): 
    prod = np.sum([v*s for (v,s) in zip(vals,state)])
    print i,prod
    if prod<mn:
        mnState = state
        mn = prod
    
print(mnState)