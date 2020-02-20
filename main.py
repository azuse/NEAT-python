import numpy as np
from typing import Optional, Sequence, List

from node import Node
from gene import Gene
from innovation import Innovation
from trait import Trait
from species import Species
from genome import Genome
from population import Population
from organism import Organism
from neat import Neat

def plot_genome(genome:Genome):
    import networkx as nx
    G = nx.DiGraph()
    for gene in genome.genes:
        nx.add_path(G, [gene.inNodeId, gene.outNodeId])
    import matplotlib.pyplot as plt
    nx.draw_networkx(G)
    plt.show()



# run one epoch for one population on pole experiment
# and create population's next generation
def pole_epoch(pop:Population, gen:int):
    # TODO

    # evaluate each organism(genome)
    win:bool = False
    winnernum:int
    for i, curgenome in enumerate(pop.genomes):
        if(pole_evaluate(curgenome)):
            win = True
        if i % 100 == 0:
            print("epoch process: {}".format(i/len(pop.genomes)))
    

    # Average and max their fitnesses for dumping to file and snapshot
    # TODO

    # Take a snapshot of the population, so that it can be visualized later on
    # TODO

    # Only print to file every print_every generations
    # TODO

    if win:
        for curgenome in pop.genomes:
            if curgenome.winner:
                # winnernum = curgenome.genome_id # TODO
                print("winner is decided!")
            

    # Create the next generation
    pop.epoch(gen)
    

    # compute average fitness

    # take a snapshot of the population for visualize TODO

    if win: 
        return 1
    else:
        return 0

# evaluate a genome's fitness
def pole_evaluate(genome:Genome):
    numnodes:int = len(genome.nodes)
    # how many visits is allowed in every activation
    thresh:int = numnodes * 2

    # TODO organism and trait
    
    # cross MAX_STEPS means it's the winner
    MAX_STEPS = 100000

    # balance the pole
    genome.fitness = go_crat(genome, MAX_STEPS, thresh)

    # print("genome fitness {}".format(genome.fitness))

    if genome.fitness > MAX_STEPS:
        genome.winner = True
        return True
    else:
        genome.winner = False
        return False
    
# pole simulator function
def go_crat(genome:Genome, max_steps:int, thresh:int):
    x:float # cart position
    x_dot:float # cart velocity
    theta:float # pole angle radians
    theta_dot:float # pole angular velocity

    steps:int = 0
    y:int

    random_start:int = 1 # start the pole in random position

    input:List[float] = [.0,.0,.0,.0,.0]
    out1:float
    out2:float

    twelve_degrees = 0.2094384

    if random_start:
        x = np.random.uniform(-2.4, 2.4)
        x_dot = np.random.uniform(-1,1)
        theta = np.random.uniform(-0.2, 0.2)
        theta_dot = np.random.uniform(-1.5, 1.5)
    else:
        x = x_dot = theta = theta_dot = 0.0

    while steps < max_steps:
        input[0] = 1 # bias
        input[1] = (x + 2.4) / 4.8
        input[2] = (x + .75) / 1.5
        input[3] = (theta + twelve_degrees) / .41
        input[4] = (theta_dot + 1) / 2

        genome.load_sensors(input)

        # activate the net (genome)
        # If it loops, exit returning only fitness of 1 step
        if not genome.activate():
            return 1

        # decide which way to push via which output unit is greater
        try:
            out1 = genome.outputs[0].activation
            out2 = genome.outputs[1].activation
        except:
            return 1

        if out1 > out2:
            y = 0
        else:
            y = 1

        # Apply action to the simulated cart-pole
        (x, x_dot, theta, theta_dot) = cart_pole(y, x, x_dot, theta, theta_dot)

        # Check for failure.  If so, return steps
        if x<-2.4 or x>2.4 or theta < -twelve_degrees or theta > twelve_degrees:
            return steps
        else:
            steps += 1
    
    return steps

# cart_pole:  Takes an action (0 or 1) and the current values of the
# four state variables and updates their values by estimating the state
# TAU seconds later.
def cart_pole(action, x, x_dot, theta, theta_dot):
    xacc:float
    thetaacc:float
    force:float
    costheta:float
    sintheta:float
    tmp:float

    GRAVITY = 9.8
    MASSCART = 1.0
    MASSPOLE = 0.1
    TOTAL_MASS = (MASSPOLE + MASSCART)
    LENGTH = 0.5 # actually half the pole's length 
    POLEMASS_LENGTH = (MASSPOLE * LENGTH)
    FORCE_MAG = 10.0
    TAU = 0.02 # seconds between state updates 
    FOURTHIRDS = 1.3333333333333

    force = FORCE_MAG if action > 0 else -FORCE_MAG
    costheta = np.cos(theta)
    sintheta = np.sin(theta)

    tmp = (force + POLEMASS_LENGTH * theta_dot * theta_dot * sintheta) / TOTAL_MASS

    thetaacc = (GRAVITY * sintheta - costheta * tmp) / (LENGTH * (FOURTHIRDS - MASSPOLE * costheta * costheta / TOTAL_MASS))

    xacc = tmp - POLEMASS_LENGTH * thetaacc * costheta / TOTAL_MASS

    # Update the four state variables, using Euler's method.
    x += TAU * x_dot
    x_dot += TAU * xacc
    theta += TAU * theta_dot
    theta_dot += TAU * thetaacc

    return (x, x_dot, theta, theta_dot)
    
                
    
# Network Parameters #
first_population_size = Neat.pop_size
# main
if __name__ == "__main__":
    # create a genome
    newgenome = Genome()
    # create full connected network in the genome
    newgenome.create_genome_full_connected(input_num=5, output_num=2, hidden_num=0)
    # plot_genome(newgenome)
    # cerate a population using the genome


    pop = Population().create_population(newgenome, first_population_size)
    # evolve through generation
    log_best_fitness = []
    gen = 0
    while gen <= 10:
        print("Generation: {}".format(gen))
        ret = pole_epoch(pop=pop, gen=gen)
        log_best_fitness.append(pop.last_fitness)
        if ret == 1:
            print("a generation has won the match!!")
            break
        gen += 1

    # plot fitness history
    # import matplotlib.pyplot as plt 
    # fig, ax = plt.subplots()
    
    # plt.plot(log_best_fitness, label="10 generation evolve example")
    # plt.yscale("symlog")
    # plt.xlabel("generation")
    # plt.ylabel("fitness")
    # from matplotlib.ticker import ScalarFormatter
    # for axis in [ax.xaxis, ax.yaxis]:
    #     axis.set_major_formatter(ScalarFormatter())
    # for i, number in enumerate(log_best_fitness):
    #     ax.annotate(number, (i, number))
    # plt.legend()
    # plt.show()

        

