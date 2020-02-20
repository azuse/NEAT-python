import numpy as np
from typing import Optional, Sequence, List

class Node:
    nodeId:int
    nodeType:str        # input output hidden
    incoming:List[int] = []  # in come link id
    outcoming:List[int] = [] # out come link id

    activation_count:int = 0 # keeps track of which activation the node is currently in
    activation:float = 0 # The total activation entering the NNode 
    last_activation:float = 0 # Holds the previous step's activation for recurrency
    last_activation2:float = 0 # Holds the activation BEFORE the prevous step's
    
    activesum:float = 0 # The incoming activity before being processed 
    activation:float = 0 # The total activation entering the NNode 
    active_flag:bool = False

    override:bool = False
    override_value:float = 0

    nodeId_old:int
    
    def new_node(self, nodeId, nodeType):
        newnode = Node()
        newnode.nodeId = nodeId
        newnode.nodeType = nodeType
        newnode.incoming = []
        newnode.outcoming = []
        return newnode

    def sensor_load(self, value:float):
        if self.nodeType == "input":
            self.last_activation2 = self.last_activation
            self.last_activation = self.activation
            self.activation = value
            self.activation_count += 1
            return True
        else:
            return False
        
    def get_active_out(self):
        if self.activation_count > 0:
            return self.activation
        else:
            return 0

    def get_active_out_td(self):
        if self.activation_count > 1:
            return self.last_activation
        else:
            return 0
        
    def active_override(self):
        self.activation = self.override_value
        self.override = False

    def override_output(self, new_output:float):
        self.override_value = new_output
        self.override = True
        
    def copy(self):
        newnode = Node()
        newnode.nodeId = self.nodeId
        newnode.nodeType = self.nodeType
        newnode.incoming = self.incoming.copy()
        newnode.outcoming = self.outcoming.copy()
        
        newnode.activation_count = self.activation_count
        newnode.activation = self.activation
        newnode.last_activation = self.last_activation
        newnode.last_activation2 = self.last_activation2

        newnode.activesum = self.activesum
        newnode.active_flag = self.active_flag

        newnode.override = self.override
        newnode.override_value = self.override_value


        return newnode