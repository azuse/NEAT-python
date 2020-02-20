import numpy as np
from typing import Optional, Sequence, List
from innovation import Innovation
from gene import Gene
from genome import Genome
from organism import Organism
from neat import Neat
# from population import Population

class Species:
    id:int
    age:int
    age_of_last_improvement:int = 0
    
    avg_fitness:float
    max_fitness:float
    max_fitness_ever:float = 0
    min_fitness:float

    genomes:List[Genome]
    # organisms:List[Organism] # TODO

    obliterate:bool = False

    expected_offspring:int = 0

    def new_species(self, id):
        newspecies = Species()
        newspecies.id = id
        newspecies.age = 1

        newspecies.avg_fitness = 0.0
        newspecies.max_fitness = 0.0
        newspecies.min_fitness = 0.0

        newspecies.genomes = []
        # newspecies.organisms = []

        return newspecies

    def add_genome(self, genome):
        self.genomes.append(genome)
        return self
    
    def remove_genome(self, genome):
        try:
            self.genomes.remove(genome)
            # print("SUCCESS: genome remove")
        except:
            # print("ERROR: genome remove error")
            pass
        return self

    def adjust_fitness(self):
        curgenome:Genome

        num_parents:int
        count:int

        age_debt:int

        age_debt = (self.age - self.age_of_last_improvement + 1) - Neat.dropoff_age

        if age_debt == 0:
            age_debt = 1

        for curgenome in self.genomes:
            curgenome.orig_fitness = curgenome.fitness
            if age_debt >= 1 or self.obliterate:
                curgenome.fitness *= 0.01

            # Give a fitness boost up to some young age (niching)
            if self.age <= 10:
                curgenome.fitness *= Neat.age_significance

            # Do not allow negative fitness
            if curgenome.fitness < 0:
                curgenome.fitness = 0.00001

            # Share fitness with the species
            curgenome.fitness = curgenome.fitness/len(self.genomes)

        # Sort the population and mark for death those after survival_thresh*pop_size
        def compare(x:Genome, y:Genome):
            if x.fitness > y.fitness:
                return -1
            elif x.fitness < y.fitness:
                return 1
            else:
                return 0
        from functools import cmp_to_key
        self.genomes = sorted(self.genomes, key=cmp_to_key(compare))

        # Update age_of_last_improvement here
        if self.genomes[0].orig_fitness > self.max_fitness_ever:
            self.max_fitness_ever = self.genomes[0].orig_fitness

        # Decide how many get to reproduce based on survival_thresh*pop_size
        # Adding 1.0 ensures that at least one will survive
        from math import floor
        num_parents = floor(Neat.survival_thresh * len(self.genomes) + 1)

        # Mark the champion
        self.genomes[0].champion = True
        # Mark for death those who are ranked too low to be parents
        for i in range(num_parents, len(self.genomes)):
            self.genomes[i].eliminate = True

    def count_offspring(self, skim):
        curgenome:Genome
        e_o_intpart:int
        e_o_fracpart:float
        skim_intpart:float
        self.expected_offspring = 0

        for curgenome in self.genomes:
            from math import floor
            e_o_intpart = floor(curgenome.expected_offspring)
            from math import fmod
            e_o_fracpart = fmod(curgenome.expected_offspring, 1.0)
            
            self.expected_offspring += e_o_intpart

            # Skim off the fractional offspring
            skim += e_o_fracpart

            # NOTE:  Some precision is lost by computer Must be remedied later
            if skim > 1.0:
                skim_intpart = floor(skim)
                self.expected_offspring += skim_intpart
                skim -= skim_intpart

        return skim

    def reproduce(self, generation:int, pop, sorted_species:List):
        count:int
        curgenome:Genome
        poolsize:int # number of genome in oldgeneration
        genomenum:int
        genomecount:int

        mon:Genome
        dad:Genome
        baby:Genome

        new_genome:Genome

        curspecies:Species
        newspecies:Species
        compgenome:Genome

        randspecies:Species
        randmult:float
        randspeciesnum:int
        spcount:int
        cursp:Species

        pause:int

        outside:bool
        found:bool
        champ_done:bool = False

        thechamp:Genome

        giveup:int
        
        mut_struct_baby:bool
        mate_baby:bool

        mut_power:float = Neat.weight_mut_power

        total_fitness:float = 0.0
        marble:float # The marble will have a number between 0 and total_fitness
        spin:float # 0Fitness total while the wheel is spinning

        # check for mistake
        if self.expected_offspring > 0 and len(self.genomes) == 0:
            print("ERROR: Attempt to reproduce out of empty species")
            return False

        poolsize = len(self.genomes) - 1
        thechamp = self.genomes[0]

        # create the designated number of offspring for the Species
        # one at a time
        from math import floor
        self.expected_offspring = floor(self.expected_offspring)
        for count in range(0, self.expected_offspring):
            mut_struct_baby = False
            mate_baby = False

            outside = False

            if self.expected_offspring > Neat.pop_size:
                print("Alert: expected offspring > Neat.popszie")
                # input()
            
            # if we have a super champ, finish off some special clones
            if thechamp.super_champ_offspring > 0:
                mom = thechamp
                new_genome = mom.copy() # mom.duplicate(count) count is id, my genome has no id implement so

                if thechamp.super_champ_offspring == 1:
                    pass
                
                if thechamp.super_champ_offspring > 1:
                    if np.random.uniform(0,1) <0.8 or Neat.mutate_add_link_prob == 0.0:
                        new_genome.mutate_link_weights(mut_power)
                    else:
                        new_genome.mutate_add_link(Neat.newlink_tries, pop.innovations, pop.innovations[-1].innovation_num1 + 1)
                        mut_struct_baby = True

                # baby = new Organism(0.0, new_genome, generation) # TODO
                baby = new_genome
                

                # pop_champ # TODO

            elif not champ_done and self.expected_offspring > 5:
                mom = thechamp
                new_genome = mom.copy()
                baby = new_genome
                champ_done = True
            # First, decide whether to mate or mutate If there is only one organism in the pool, then always mutate
            elif np.random.uniform(0,1) < Neat.mutate_only_prob or poolsize == 0:
                # choose a random parent
                genomenum = np.random.randint(0, poolsize+1)
                curgenome = self.genomes[genomenum]

                mom = curgenome

                new_genome = mom.copy()

                # do mutation depending on probabilityies of various mutations 
                if np.random.uniform(0,1) < Neat.mutate_add_node_prob:
                    new_genome.mutate_add_node(innovations=pop.innovations, curinnov_num=len(pop.innovations))
                    mut_struct_baby = True
                elif np.random.uniform(0,1) < Neat.mutate_add_link_prob:
                    new_genome.mutate_add_link(innovations=pop.innovations, curinnov_num=len(pop.innovations), tries=Neat.newlink_tries)
                    mut_struct_baby = True
                else:
                    # not structural mutation
                    if np.random.uniform(0,1) < Neat.mutate_link_weights_prob:
                        new_genome.mutate_link_weights(Neat.weight_mut_power)
                baby = new_genome
            else:
                # Otherwise we should mate 
                # mating TODO
                # choose random mom
                genomenum = np.random.randint(0, len(self.genomes))
                curgenome = self.genomes[genomenum]
                mom = curgenome

                # choose random dad(within species)
                genomenum = np.random.randint(0, len(self.genomes))
                curgenome = self.genomes[genomenum]
                dad = curgenome

                # //Perform mating based on probabilities of differrent mating types # TODO
                # just mate
                new_genome = mom.mate_multipoint(g=dad, genomeId=count, fitness1=mom.orig_fitness, fitness2=dad.orig_fitness, interspec_flag=outside)
                mate_baby = True

                # determin wether to mutate baby's genome
                if np.random.uniform(0,1) > Neat.mate_only_prob or dad == mom or dad.compatibility(mom) == 0:
                    # do mutation depending on probabilityies of various mutations 
                    if np.random.uniform(0,1) < Neat.mutate_add_node_prob:
                        new_genome.mutate_add_node(innovations=pop.innovations, curinnov_num=len(pop.innovations))
                        mut_struct_baby = True
                    elif np.random.uniform(0,1) < Neat.mutate_add_link_prob:
                        new_genome.mutate_add_link(innovations=pop.innovations, curinnov_num=len(pop.innovations), tries=Neat.newlink_tries)
                        mut_struct_baby = True
                    else:
                        # not structural mutation
                        if np.random.uniform(0,1) < Neat.mutate_link_weights_prob:
                            new_genome.mutate_link_weights(Neat.weight_mut_power)
                    baby = new_genome
                else:
                    baby = new_genome

            # add baby to its proper species
            baby.mut_struct_baby = mut_struct_baby
            baby.mate_baby = mate_baby

            if len(pop.species) == 0:
                # create first species
                newspecies = Species().new_species(pop.last_species + 1)
                pop.last_species += 1
                pop.species.append(newspecies)
                newspecies.add_genome(baby)
                baby.species = newspecies

            else:
                for curspecies in pop.species:
                    found = False
                    if len(curspecies.genomes)> 0:
                        compgenome = curspecies.genomes[0]
                        if baby.compatibility(compgenome) < Neat.compat_thresh:
                            curspecies.add_genome(baby)
                            baby.species = curspecies
                            found = True
                            break
                
                if not found:
                    newspecies = Species().new_species(pop.last_species + 1)
                    pop.last_species += 1
                    pop.species.append(newspecies)
                    newspecies.add_genome(baby)
                    baby.species = curspecies

            
        return True



        


                





