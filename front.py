# -*- coding: utf-8 -*-
"""
Created on Thu Mar 02 18:17:21 2017

@author: simon
"""

class Front(object):
    
    def __init__(self,nodes,delay=None,Q=None,distance=None):
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
        
        if type(distance) is not list:
            distance = [distance]
        self.distance = distance
        assert len(distance)==self.front_size
        
        self.current_from = []
        
        self.previous_nodes = []
        self.previous_delays = []
        self.previous_distances = []
        self.previous_Qs = []
        
        self.next_nodes = []
        self.next_delays = []
        self.next_distances = []
        self.next_Qs = []
        
        self.node_history = []
        self.nstep = 0
        
    def get_current_front(self):
        return self.current_nodes, self.delay, self.Q, self.distance
    
    def step_front(self,nodes,delay=None,Q=None,distance=None):
        self.next_nodes.extend(nodes)
        self.next_delays.extend(delay)
        self.next_Qs.extend(Q)
        self.next_distances.extend(distance)
        
    def complete_step(self):
        self.previous_nodes = self.current_nodes
        self.previous_delay = self.delay
        self.previous_Qs = self.Q
        self.previous_distances = self.distance
        
        self.current_nodes = self.next_nodes
        self.delay = self.next_delays
        self.Q = self.next_Qs
        self.distance = self.next_distances
        
        self.next_nodes = []
        self.next_delays = []
        self.next_Qs = []
        self.next_distances = []
        
        self.front_size = len(self.current_nodes)
        self.nstep += 1