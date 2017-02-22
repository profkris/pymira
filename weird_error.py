# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 08:03:04 2017

@author: simon
"""

class A(object):
    
    def __init__(self, attr=[]):
        self.attr = attr
        
A1 = A()
A2 = A()
tmp = []
A1.attr = tmp
A2.attr = tmp
A1.attr.append(['new entry'])
print(A1.attr)
print(A2.attr)