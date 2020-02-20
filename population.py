import numpy as np
from typing import Optional, Sequence, List
from innovation import Innovation
from organism import Organism
from genome import Genome
from species import Species
from neat import Neat

class Population:
    # organisms:List[Organism] = [] # TODO
    species:List[Species]
    genomes:List[Genome]
    innovations:List[Innovation]

    winnergen:int = 0
    highest_fitness:float = 0.0
    last_fitness:float = 0.0
    highest_last_changed:int = 0

    cur_node_id:int = 0
    cur_innov_num:int = 0

    last_species:int = 0
    def __init__(self):
        self.species:List[Species] = []
        self.genomes:List[Genome] = []
        self.innovations:List[Innovation] = []

    def create_population(self, genome:Genome, size):
        # create species include the genome
        newpopulation = Population()
        newpopulation.spawn(genome, size)
        return newpopulation

    # mutate and spawn new genomes
    def spawn(self, genome:Genome, size):
        new_genome:genome
        for count in range(0,size):
            new_genome = genome.copy()
            # mutate
            new_genome.mutate_link_weights(power=Neat.weight_mut_power, mode="COLDGAUSSPOINT")
            # randomize traits
            # TODO

            # create new organism
            # TODO
            self.genomes.append(new_genome)

        # keep record of current node id and innovation num
        self.cur_innov_num = new_genome.get_last_gene_innovnum
        self.cur_node_id = new_genome.get_last_node_id

        # seperate new population into species
        self.speciate()

        return True

    # speciate the genomes based on distence
    def speciate(self):
        curgenome:Genome  # For stepping through Population
        curspecies:Species  # Steps through species
        comgenome:Genome  # Organism for comparison 
        newspecies:Species  # For adding a new species

        counter:int = 0  # Species counter


        # Step through all existing organisms
        for curgenome in self.genomes:
            # search a species for each genome
            if len(self.species) == 0:
                # Create the first species
                counter += 1
                newspecies = Species().new_species(counter)
                self.species.append(newspecies)
                newspecies.add_genome(curgenome) # add current genome
                curgenome.species = newspecies
            else:
                for curspecies in self.species:
                    if len(curspecies.genomes) > 0:
                        for comgenome in curspecies.genomes:
                            if curgenome.compatibility(comgenome) < Neat.compat_thresh:
                                # Found compatible species, so add this organism to it
                                curspecies.add_genome(curgenome)
                                curgenome.species = curspecies
                                comgenome = 0 # mark search over
                                break
                            else:
                                # Keep searching for a matching species
                                break
                        
                        if comgenome == 0:
                            break
                # If we didn't find a match, create a new species
                if comgenome != 0:
                    counter += 1
                    newspecies = Species().new_species(counter)
                    self.species.append(newspecies)
                    newspecies.add_genome(curgenome)
                    curgenome.species = newspecies

        self.last_species = counter # Keep track of highest species
        return True
        
    def epoch(self, gen:int):
        curspecies:Species
        deadspecies:Species  # For removing empty Species

        curgenome:Genome
        deadgenome:Genome

        curinnov:Innovation
        deadinnov:Innovation

        total:float = 0 # compute average fitness on all genomes
        overall_average:float = 0 # average modified fitness on all genomes

        gencount:int = 0


        skim:float
        total_expected:int 
        total_genomes:int = len(self.genomes)
        max_expected:int
        best_species:Species
        final_expected:int = 0

        pause:int

        # right to make babies can be stolen from bad species and give to good ones
        NUM_STOLEN:int = Neat.babies_stolen
        one_fifth_stolen:int
        one_tenth_stolen:int

        stored_species:List[Species] # species stored by max fit genome in species
        stolen_babies:int # baby stolen from bad species and give to champion

        half_pop:int

        best_species_num:int # use to debug
        best_ok:bool

        num_species_target:int = 4 # try to keep species at this number
        num_species:int = len(self.species)
        compat_mod:float = 0.3 # modify compat thresh to control speciation
        
        def compare(a:Species, b:Species):
            if a.genomes[0].orig_fitness > b.genomes[0].orig_fitness:
                return -1
            elif a.genomes[0].orig_fitness < b.genomes[0].orig_fitness:
                return 1
            else:
                return 0

        from functools import cmp_to_key
        sorted_species:List[Species] = sorted(self.species, key=cmp_to_key(compare))
        
        # Flag the lowest performing species over age 20 every 30 generations 
        for curspecies in reversed(sorted_species):
            if curspecies.age >= 20:
                break

        if gen % 30 == 0:
            curspecies.obliterate = True

        print("Number of species: {}".format(num_species))
        print("compat_threth: {}".format(Neat.compat_thresh))

        for curspecies in self.species:
            curspecies.adjust_fitness()

        # Go through the organisms(genomes) and add up their fitnesses to compute the overall average
        for curgenome in self.genomes:
            total += curgenome.fitness
        
        if total_genomes == 0:
            overall_average = 0
        else:
            overall_average = total/total_genomes
        print("Generation = {} ,Overall average = {}".format(gen, overall_average))

        # Now compute expected number of offspring for each individual organism(genome)
        for curgenome in self.genomes:
            curgenome.expected_offspring = curgenome.fitness/overall_average

        # Now add those offspring up within each Species to get the number of offspring per Species
        skim = 0
        total_expected = 0
        for curspecies in self.species:
            skim = curspecies.count_offspring(skim)
            total_expected += curspecies.expected_offspring

        # Need to make up for lost foating point precision in offspring assignment If we lost precision, give an extra baby to the best Species
        if total_expected < total_genomes:
            # Find the Species expecting the most
            max_expected = 0
            final_expected = 0
            for curspecies in self.species:
                if curspecies.expected_offspring >= max_expected:
                    max_expected = curspecies.expected_offspring
                    best_species = curspecies
                final_expected += curspecies.expected_offspring
            # Give the extra offspring to the best species
            best_species.expected_offspring += 1
            final_expected += 1

            # If we still arent at total, there is a problem
            # Note that this can happen if a stagnant Species
            # dominates the population and then gets killed off by its age
            # Then the whole population plummets in fitness
            # If the average fitness is allowed to hit 0, then we no longer have 
            # an average we can use to assign offspring.
            if final_expected < total_expected:
                print("Population diead!")
                for curspecies in self.species:
                    curspecies.expected_offspring = 0

                best_species.expected_offspring = total_genomes  # give all offspring to the best species

        # sort species by max fitness 
        def compare_species(x:Species, y:Species):
            if x.genomes[0].orig_fitness > y.genomes[0].orig_fitness:
                return -1
            elif x.genomes[0].orig_fitness < y.genomes[0].orig_fitness:
                return 1
            else:
                return 0
                
        from functools import cmp_to_key
        sorted_species:List[Species] = sorted(sorted_species, key=cmp_to_key(compare_species))

        best_species_num = sorted_species[0].id

        # print debug info
        for curspecies in reversed(sorted_species):
            print("Species id: {} (Size: {}) , orig fitness: {} , last improved: {}".format(curspecies.id, len(curspecies.genomes), curspecies.genomes[0].orig_fitness, curspecies.age - curspecies.age_of_last_improvement))

        # plot best genome
        Neat().plot_genome(sorted_species[0].genomes[0], gen=gen)
        Neat().print_genome_to_file(sorted_species[0].genomes[0])
        # Check for Population-level stagnation
        curspecies = sorted_species[0]
        if curspecies.genomes[0].orig_fitness > self.highest_fitness:
            # the population improved
            self.highest_fitness = curspecies.genomes[0].orig_fitness
            self.highest_last_changed = 0
            print("NEW POPULATION FITNESS RECORD: {}".format(self.highest_fitness))
        else:
            # the population stagnated
            self.highest_last_changed += 1
            print("{} generations since last population fitness record".format(self.highest_last_changed))
        self.last_fitness = curspecies.genomes[0].orig_fitness

        # Check for stagnation- if there is stagnation, perform delta-coding
        if self.highest_last_changed > Neat.dropoff_age + 5:
            print("Performing delta coding")

            half_pop = Neat.pop_size/2

            curspecies = sorted_species[0]
            curspecies.genomes[0].super_champ_offspring = half_pop
            curspecies.expected_offspring = half_pop
            curspecies.age_of_last_improvement = curspecies.age

            if len(sorted_species) > 2:
                curspecies = sorted_species[1]
                curspecies.genomes[0].super_champ_offspring = Neat.pop_size - half_pop
                curspecies.expected_offspring = Neat.pop_size - half_pop
                curspecies.age_of_last_improvement = curspecies.age

                # Get rid of all species under the first 2
                for i in range(2,len(sorted_species)):
                    sorted_species[i].expected_offspring = 0

            elif len(sorted_species) == 2:
                curspecies = sorted_species[1]
                curspecies.genomes[0].super_champ_offspring += Neat.pop_size - half_pop
                curspecies.expected_offspring = Neat.pop_size - half_pop
        
        elif Neat.babies_stolen > 0:
            # stole babies
            stolen_babies = 0
            for curspecies in reversed(sorted_species):
                if stolen_babies < NUM_STOLEN:
                    # print("Consider stealing Species {}, age: {}, expected_offspring: {}".format(curspecies.id, curspecies.age, curspecies.expected_offspring))
                    if curspecies.age > 5 and curspecies.expected_offspring > 2:
                        # print("Stealing! ")

                        if curspecies.expected_offspring - 1 >= NUM_STOLEN-stolen_babies:
                            # This species has enough to finish off the stolen pool
                            curspecies.expected_offspring -= (NUM_STOLEN - stolen_babies)
                            stolen_babies = NUM_STOLEN
                        else:
                            # Not enough here to complete the pool of stolen
                            stolen_babies += curspecies.expected_offspring - 1
                            curspecies.expected_offspring = 1
            
            one_fifth_stolen = Neat.babies_stolen / 5
            one_tenth_stolen = Neat.babies_stolen / 10

            i = 0 # current species pointer
            # Don't give to dying species even if they are champs
            while i != len(sorted_species) and sorted_species[i].age_of_last_improvement > Neat.dropoff_age:
                i += 1

            # Concentrate A LOT on the number one speciesd
            if i != len(sorted_species) and stolen_babies >= one_fifth_stolen:
                sorted_species[i].genomes[0].super_champ_offspring = one_fifth_stolen
                sorted_species[i].expected_offspring += one_fifth_stolen
                stolen_babies -= one_fifth_stolen
                i += 1
            
            # Don't give to dying species even if they are champs
            while i != len(sorted_species) and sorted_species[i].age_of_last_improvement > Neat.dropoff_age:
                i += 1

            
            if i != len(sorted_species) and stolen_babies >= one_fifth_stolen:
                sorted_species[i].genomes[0].super_champ_offspring = one_fifth_stolen
                sorted_species[i].expected_offspring += one_fifth_stolen
                stolen_babies -= one_fifth_stolen
                i += 1

            # Don't give to dying species even if they are champs
            while i != len(sorted_species) and sorted_species[i].age_of_last_improvement > Neat.dropoff_age:
                i += 1

            if i != len(sorted_species) and stolen_babies >= one_tenth_stolen:
                sorted_species[i].genomes[0].super_champ_offspring = one_tenth_stolen
                sorted_species[i].expected_offspring += one_tenth_stolen
                stolen_babies -= one_tenth_stolen
                i += 1

            # Don't give to dying species even if they are champs
            while i != len(sorted_species) and sorted_species[i].age_of_last_improvement > Neat.dropoff_age:
                i += 1

            while i != len(sorted_species) and stolen_babies > 0:
                if np.random.uniform(0,1) > 0.1:
                    if stolen_babies > 3:
                        sorted_species[i].genomes[0].super_champ_offspring = 3
                        sorted_species[i].expected_offspring += 3
                        stolen_babies -= 3
                    else:
                        sorted_species[i].genomes[0].super_champ_offspring = stolen_babies
                        sorted_species[i].expected_offspring += stolen_babies
                        stolen_babies = 0

                    i += 1
                    # Don't give to dying species even if they are champs
                    while i != len(sorted_species) and sorted_species[i].age_of_last_improvement > Neat.dropoff_age:
                        i += 1
                    
            # If any stolen babies aren't taken, give them to species #1's champ
            if stolen_babies > 0:
                curspecies = sorted_species[0]
                curspecies.genomes[0].super_champ_offspring += stolen_babies
                curspecies.expected_offspring += stolen_babies
                stolen_babies = 0

        # Kill off all Organisms marked for death.  The remainder will be allowed to reproduce.
        for curgenome in self.genomes:
            if curgenome.eliminate:
                for curspecies in self.species:
                    try: 
                        curspecies.genomes.remove(curgenome)
                    except:
                        pass
                deadgenome = curgenome
                self.genomes.remove(curgenome)
            
        print("REPRODUCING")
        last_id = self.species[0].id
        for curspecies in self.species: 
            curspecies.reproduce(generation=gen, pop=self, sorted_species=sorted_species)

            # Set the current species to the id of the last species checked
            # TODO
        
        print("Reproducing completed")

        # destroy and remove the old generation from the organism and species
        for curgenome in self.genomes:
            for curspecies in self.species:
                try:
                    curspecies.remove_genome(genome=curgenome)
                except:
                    pass
        self.genomes = [] # remove all old genome

        # remove empty species and age ones
        i = 0
        while i < len(self.species):
            curspecies = self.species[i]
            if len(curspecies.genomes) == 0:
                self.species.remove(curspecies)
                continue
            else:
                # Age any Species that is not newly created in this generation #TODO
                curspecies.age += 1

                # Go through the organisms of the curspecies and add them to the master list
                for curgenome in curspecies.genomes:
                    # curgenome.id = genomecount++ # my genome has no id
                    self.genomes.append(curgenome)
                
                i += 1
        

        print ("Epoch complete")
        return True






