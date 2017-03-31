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
        
        self.front_size = len(nodes)        
        self.capacity = 1000*self.front_size
        
        if type(delay) is not list:
            delay = [delay]
        self.delay = delay
        
        if type(Q) is not list:
            Q = [Q]
        self.Q = Q
        
        if type(distance) is not list:
            distance = [distance]
        self.distance = distance
        
        self.current_from = []
        
        self.previous_nodes = []
        self.previous_delays = []
        self.previous_distances = []
        self.previous_Qs = []
        
        # Preallocate memory
        self.next_front_size = 0
        self.next_nodes = [None]*self.capacity
        self.next_delays = [None]*self.capacity
        self.next_distances = [None]*self.capacity
        self.next_Qs = [None]*self.capacity
        
        self.node_history = []
        self.nstep = 0
        
    def _increase_capacity(self,min_size=None,reset=False):
        while True:
            self.capacity *= 10.
            if min_size is None or self.capacity>min_size:
                break

        if not reset:            
            ext_size = self.capacity-self.front_size
        else:
            ext_size = self.capacity
            
        self.next_nodes.extend([None]*ext_size)
        self.next_delays.extend([None]*ext_size)
        self.next_Qs.extend([None]*ext_size)
        self.next_distances.extend([None]*ext_size)
        
    def get_current_front(self):
        return self.current_nodes, self.delay, self.Q, self.distance
    
    def step_front(self,nodes,delay=None,Q=None,distance=None):
        
        n_to_add = len(nodes)
        if n_to_add>self.capacity:
            self._increase_capacity(min_size=self.front_size+n_to_add)
            
        self.next_nodes[self.next_front_size:self.next_front_size+n_to_add] = nodes
        self.next_delays[self.next_front_size:self.next_front_size+n_to_add] = delay
        self.next_Qs[self.next_front_size:self.next_front_size+n_to_add] = Q
        self.next_distances[self.next_front_size:self.next_front_size+n_to_add] = distance
        self.next_front_size += n_to_add
        
    def complete_step(self):
        self.previous_nodes = self.current_nodes
        self.previous_delay = self.delay
        self.previous_Qs = self.Q
        self.previous_distances = self.distance
        
        self.current_nodes = self.next_nodes[0:self.next_front_size]
        self.delay = self.next_delays[0:self.next_front_size]
        self.Q = self.next_Qs[0:self.next_front_size]
        self.distance = self.next_distances[0:self.next_front_size]
        self.front_size = len(self.current_nodes)
        
        self.next_nodes = [None] * self.capacity
        self.next_delays = [None] * self.capacity
        self.next_Qs = [None] * self.capacity
        self.next_distances = [None] * self.capacity
        self.next_front_size = 0

        self.nstep += 1