import numpy as np
from typing import Optional, Sequence


class Gene:
    traitnum:int
    innovation_num:int
    mutation_num:int = 0
    enabled:bool
    frozen:bool # if a gene is frozen its weight can't be changed


    # gene is a link between nodes
    inNodeId:int
    outNodeId:int
    weight:float
    is_recurrent:bool
    time_delay:bool = False
    added_weight:float = 0 # the amount of weight adjsutment

    def add_gene_with_trait(self, traitnum, weight, inodeId, onodeId, recur_flag, innov, mnum):
        newgene = Gene()
        newgene.weight = weight
        newgene.inNodeId = inodeId
        newgene.outNodeId = onodeId
        newgene.is_recurrent = recur_flag
        newgene.innovation_num = innov
        newgene.mutation_num = mnum
        
        newgene.enabled = True
        newgene.frozen = False

        newgene.traitnum = traitnum
        return newgene


    def add_gene_with_no_trait(self, weight, inodeId, onodeId, recur_flag, innov, mnum):
        newgene = Gene()
        newgene.weight = weight
        newgene.inNodeId = inodeId
        newgene.outNodeId = onodeId
        newgene.is_recurrent = recur_flag
        newgene.innovation_num = innov
        newgene.mutation_num = mnum
        
        newgene.enabled = True
        newgene.frozen = False

        return newgene

    def add_gene_from_file(self):
        pass
    def print_gene_to_file(self):
        pass

    def copy(self):
        newgene = Gene()
        newgene.innovation_num = self.innovation_num
        newgene.enabled = self.enabled
        newgene.frozen = self.frozen

        newgene.inNodeId = self.inNodeId
        newgene.outNodeId = self.outNodeId
        newgene.weight = self.weight
        newgene.is_recurrent = self.is_recurrent
        newgene.time_delay = self.time_delay
        newgene.added_weight = self.added_weight

        return newgene

