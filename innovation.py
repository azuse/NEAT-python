import numpy as np
from typing import Optional, Sequence

# Innovation class is a class used to store and track all the mutation history in the genome
class Innovation:
    innovation_type:str = "" # "node" "link"
    
    node_in_id:int = 0
    node_out_id:int = 0

    innovation_num1:int = 0
    innovation_num2:int = 0# if the innovation is a node, need two innovation numbers for the node and the link to the node

    new_weight:float = 0
    new_traitnum:int = 0

    newnode_id:int = 0

    old_innovation_num:int = 0 #If a new node was created, this is the innovnum of the gene's link it is being stuck inside

    recur_flag:bool = False

    # construct new node
    def new_node(self, nin:int, nout:int, num1:int, num2:int, newid:int, oldinnov:int):
        newinnov = Innovation()
        newinnov.innovation_type = "node"
        newinnov.node_in_id = nin
        newinnov.node_out_id = nout
        newinnov.innovation_num1 = num1
        newinnov.innovation_num2 = num2
        newinnov.newnode_id = newid
        newinnov.old_innovation_num = oldinnov
        
        return newinnov
        
    # construct new link
    def new_link(self, nin:int, nout:int, num1:int, weight:float, recur:bool=False):
        newinnov = Innovation()        
        newinnov.innovation_type="link"
        
        newinnov.node_in_id=nin
        newinnov.node_out_id=nout
        
        newinnov.innovation_num1 = num1
        
        newinnov.new_weight = weight

        newinnov.recur_flag = recur

        return newinnov