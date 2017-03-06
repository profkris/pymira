# -*- coding: utf-8 -*-
"""
Created on Thu Mar 02 18:17:21 2017

@author: simon
"""

class Front(object):
    
    def __init__(self,nodes,delay=None,Q=None):
        if type(nodes) is not list:
            nodes = [nodes]
        self.current_nodes = nodes
        
        self.front_size = len(self.current_nodes)
        
        if type(delay) is not list:
            delay = [delay]
        self.delay = delay
        assert len(delay)==self.front_size
        
        if type(Q) is not list:
            Q = [Q]
        self.Q = Q
        assert len(Q)==self.front_size
        
        self.current_from = []
        
        self.previous_nodes = []
        self.previous_delays = []
        self.previous_Qs = []
        
        self.next_nodes = []
        self.next_delays = []
        self.next_Qs = []
        
        self.node_history = []
        self.nstep = 0
        
    def get_current_front(self):
        return self.current_nodes, self.delay, self.Q
    
    def step_front(self,nodes,delay=None,Q=None):
        self.next_nodes.extend(nodes)
        self.next_delays.extend(delay)
        self.next_Qs.extend(Q)
        
    def complete_step(self):
        self.previous_nodes = self.current_nodes
        self.previous_delay = self.delay
        self.previous_Qs = self.Q
        
        self.current_nodes = self.next_nodes
        self.delay = self.next_delays
        self.Q = self.next_Qs
        
        self.next_nodes = []
        self.next_delays = []
        self.next_Qs = []
        
        self.front_size = len(self.current_nodes)
        self.nstep += 1