import numpy as np
from typing import Optional, Sequence, List


class Neat:
    # parameters
    trait_param_mut_prob = 0.5 
    trait_mutation_power = 1.0
    linktrait_mut_sig = 1.0
    nodetrait_mut_sig = 0.5
    weight_mut_power = 1
    recur_prob = 0.05
    disjoint_coeff = 1.0
    excess_coeff = 1.0
    mutdiff_coeff = 7.0
    compat_thresh = 9.0
    age_significance = 1.0
    survival_thresh = 0.4
    mutate_only_prob = 0.9
    mutate_random_trait_prob = 0.1
    mutate_link_trait_prob = 0.1
    mutate_node_trait_prob = 0.1
    mutate_link_weights_prob = 0.5
    mutate_toggle_enable_prob = 0.0
    mutate_gene_reenable_prob = 0.00
    mutate_add_node_prob = 0.3
    mutate_add_link_prob = 0.3
    interspecies_mate_rate = 0.01
    mate_multipoint_prob = 0.6
    mate_multipoint_avg_prob = 0.4
    mate_singlepoint_prob = 0.0
    mate_only_prob = 0.2
    recur_only_prob = 0.2
    pop_size = 1000
    dropoff_age = 15
    newlink_tries = 20
    print_every = 60
    babies_stolen = 0
    num_runs = 1

    def plot_genome(self, genome, gen=0):

        import matplotlib.pyplot as plt
        import networkx as nx
        from networkx.drawing.nx_agraph import graphviz_layout

        G = nx.DiGraph()
        ed = []

        for gene in genome.genes:
            if gene.enabled:
                ed.append([gene.inNodeId, gene.outNodeId, round(gene.weight, 4)])

        G.add_weighted_edges_from(ed)
        pos = graphviz_layout(G, prog='dot', args="-Grankdir=LR")
        name = "Generation: {}, Fitness: {}, Orig_fitness: {}".format(gen, genome.fitness, genome.orig_fitness)
        nx.draw(G,pos,with_labels=True,font_weight='bold',label=name)
        labels = nx.get_edge_attributes(G,'weight')
        nx.draw_networkx_edge_labels(G, with_labels=True, pos=pos, font_weight='bold', edge_labels=labels, font_size=6,label=name)
        plt.legend(loc='upper right')
        plt.show()

    def print_genome_to_file(self, genome):
        import time
        f = open("genome.{}.log".format(time.time()), "w")
        f.write("fitness orig_fitness\n")
        f.write("{} {}\n".format(genome.fitness, genome.orig_fitness))
        f.write("node\n")

        for node in genome.nodes:
            f.write("{} {} {} {}\n".format(node.nodeId, node.nodeType, node.incoming, node.outcoming))
        
        f.write("gene\n")

        for gene in genome.genes:
            f.write("{} {} {} {} {} {}\n".format(gene.innovation_num, gene.inNodeId, gene.outNodeId, gene.weight, gene.is_recurrent, gene.enabled))

        f.close()
    
    def read_genome_from_file(self, filename):
        f = open(filename, "r")
        f.readline()
        line = f.readline()
        var = line.split(" ")
        fitness = var[0]
        orig_fitness = var[1]
        f.readline()
        nodes = []
        while line != "gene":
            line = f.readline() 
            if line == "gene":
                break
            var = line.split(" ")
            from node import Node
            newnode = Node().new_node(var[0], var[1])
            nodes.append(newnode)
        genes = []
        while line != "":
            line = f.readline()
            if line == "":
                break
            var = line.split(" ")
            from gene import Gene
            newgene = Gene().add_gene_with_no_trait()
            genes.append(newgene)

        
            




