import numpy as np
from typing import Optional, Sequence
from genome import Genome
from network import Network
# from species import Species

class Organism:
    fitness:float = 0.0
    orig_fitness:float = 0.0 # a fitness won't change in adjustments
    winner:bool = False
    
    net:Network     # Organism's phenotype
    genome:Genome   # Organism's genotype
    # species:Species

    expected_offspring:int = 0 # number of children organism may have
    generation:int  # which generation this organism is from

    # TODO
