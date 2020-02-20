import numpy as np
from typing import Optional, Sequence

# Trait: a trait contain a list of parameter. each node points to a trait and traits can evolve on their own
class Trait:
    trait_id:int
    params:list
    def new_trait(self, id:int, *args):
        newtrait = Trait()
        newtrait.trait_id = id
        newtrait.params = []
        for arg in args:
            newtrait.params.append(arg)
        return newtrait