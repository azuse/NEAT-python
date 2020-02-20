import numpy as np
from typing import Optional, Sequence, List
from gene import Gene
from innovation import Innovation
from node import Node
from trait import Trait
from neat import Neat

class Genome:
    traits:List[Trait]
    nodes:List[Node]
    genes:List[Gene]
    recurrent_probability = 0.5
    loop_probability = 0.5

    inputs:List[Node]
    hiddens:List[Node]
    outputs:List[Node]

    fitness:float = 0.0 # A measure of fitness for the Organism
    orig_fitness:float = 0.0 # A fitness measure that won't change during adjustments
    winner:bool = False # is this organism winner?
    champion:bool = False # champion in the species epoch
    eliminate:bool = False  # eliminate this genome in the species epoch

    expected_offspring:int = 0
    super_champ_offspring:int = 0

    # species:Species # which species this genome is belong to

    mut_struct_baby:bool = False # is the genome a mutated baby in reproducing
    mate_baby:bool = False  # is the genome a mated baby in reproducing
    def __init__(self):
        self.traits:List[Trait] = []
        self.nodes:List[Node] = []
        self.genes:List[Gene] = []
        self.inputs:List[Node] = []
        self.hiddens:List[Node] = []
        self.outputs:List[Node] = []
    
    def create_genome_full_connected(self, input_num, output_num, hidden_num):
        cur_node_id = 0

        # create input nodes
        for i in range(0, input_num):
            newnode = Node().new_node(nodeId=cur_node_id, nodeType="input")
            self.inputs.append(newnode)
            self.nodes.append(newnode)
            cur_node_id += 1

        # create hidden nodes
        for i in range(0, hidden_num):
            newnode = Node().new_node(nodeId=cur_node_id, nodeType="hidden")
            self.hiddens.append(newnode)
            self.nodes.append(newnode)
            cur_node_id += 1

        # create output nodes
        for i in range(0, output_num):
            newnode = Node().new_node(nodeId=cur_node_id, nodeType="output")
            self.outputs.append(newnode)
            self.nodes.append(newnode)
            cur_node_id += 1

        if hidden_num == 0:
            # no hidden nodes for genome, connect input to output
            count = 1
            for input in self.inputs:
                for output in self.outputs:
                    newgene = Gene().add_gene_with_no_trait(weight=0, inodeId=input.nodeId, onodeId=output.nodeId, recur_flag=False, innov=count, mnum=0)
                    self.genes.append(newgene)
                    input.outcoming.append(self.genes.index(newgene))
                    output.incoming.append(self.genes.index(newgene))
                    count+=1
            return 0
        else:
            # create full connect links for input->hidden
            count = 1
            for input in self.inputs:
                for hidden in self.hiddens:
                    newgene = Gene().add_gene_with_no_trait(weight=0, inodeId=input.nodeId, onodeId=hidden.nodeId, recur_flag=False, innov=count, mnum=0)
                    self.genes.append(newgene)
                    input.outcoming.append(self.genes.index(newgene))
                    hidden.incoming.append(self.genes.index(newgene))
                    count+=1
            
            # create full connect links for hidden->output
            for hidden in self.hiddens:
                for output in self.outputs:
                    newgene = Gene().add_gene_with_no_trait(weight=0, inodeId=hidden.nodeId, onodeId=output.nodeId, recur_flag=False, innov=count, mnum=0)
                    self.genes.append(newgene)
                    hidden.outcoming.append(self.genes.index(newgene))
                    output.incoming.append(self.genes.index(newgene))
                    count+=1

            # recurrently connect all hidden nodes
            for hidden1 in self.hiddens:
                for hidden2 in self.hiddens:
                    newgene = Gene().add_gene_with_no_trait(weight=0, inodeId=hidden1.nodeId, onodeId=hidden2.nodeId, recur_flag=True, innov=count, mnum=0)
                    self.genes.append(newgene)
                    hidden1.outcoming.append(self.genes.index(newgene))
                    hidden2.incoming.append(self.genes.index(newgene))
                    count+=1

            return 0
        
        


    # mutate links' weight, add gaussian noise to the weight
    def mutate_link_weights(self, power, mode=""):
        gene:Gene
        for gene in self.genes:
            if gene.frozen:
                # don't change frozen genes
                continue
            randnum = np.random.normal(-1,1) * power
            gausspoint = 0.5
            coldgausspoint = 0.3
            randomchoice = np.random.uniform(0,1)
            if mode == "GAUSSPOINT":
                gene.weight += randnum
            elif mode == "COLDGAUSSPOINT":
                gene.weight = randnum
            elif randomchoice > gausspoint:
                gene.weight += randnum
            elif randomchoice > coldgausspoint:
                gene.weight = randnum 

            # cap the weights at 8.0
            if gene.weight > 8.0:
                gene.weight = 8.0
            elif gene.weight < -8.0:
                gene.weight = -8.0
            

            # record innovation (?)
            gene.mutation_num = gene.weight
            

    # linkMutate: add link(gene)  # TODO delete link
    def mutate_add_link(self, tries:int, innovations:List[Innovation], curinnov_num:int):
        found = False

        # decide whether is recurrent
        newGenome = self
        if np.random.normal(0,1) < self.recurrent_probability:
            recurrent = True
        else:
            recurrent = False

        # skip input nodes of the network
        first_noninput=0
        while 1:
            if(self.nodes[first_noninput].nodeType == "input"):
                first_noninput += 1
            else:
                break
        
        trycount = 0
        if recurrent:
            # add a recurrent link
            while trycount < tries:
                # whether make a loop recurrent
                if np.random.normal(0,1) < self.loop_probability:
                    loop = True
                else:
                    loop = False

                if loop:
                    node1 = self.nodes[np.random.randint(first_noninput,len(self.nodes))]
                    node2 = node1
                else:
                    node1 = self.nodes[np.random.randint(0, len(self.nodes))]
                    node2 = self.nodes[np.random.randint(first_noninput, len(self.nodes))]

                # check if the gene exists
                for gene in self.genes:
                    if gene.inNodeId == node1.nodeId and gene.outNodeId == node2.nodeId:
                        trycount += 1
                        continue
                
                # check if the gene is recurrent (prevent recurrent loop)
                deepcount = 0
                # threshold = len(self.nodes) * len(self.nodes)
                threshold = 5
                recurrent_flag = self.is_recurrent(node1, node2, threshold, deepcount)

                # consider output out of outputs recurrent
                if node1.nodeType == "output" or node2.nodeType == "output":
                    recurrent_flag = True

                if not recurrent_flag:
                    trycount += 1
                    continue
                else:
                    trycount = tries
                    found = True
                    break
        else:
            # add a nonrecurrent gene
            while trycount < tries:
                node1 = self.nodes[np.random.randint(0, len(self.nodes))]
                node2 = self.nodes[np.random.randint(first_noninput, len(self.nodes))]

                # check if the link exists
                for gene in self.genes:
                    if gene.inNodeId == node1.nodeId and gene.outNodeId == node2.nodeId:
                        trycount += 1
                        continue

                # check if the link is recurrent (prevent recurrent loop)
                deepcount = 0
                # threshold = len(self.nodes) * len(self.nodes)
                threshold = 5
                recurrent_flag = self.is_recurrent(node1, node2, threshold, deepcount)

                # consider output out of outputs recurrent
                if node1.nodeType == "output" or node2.nodeType == "output":
                    recurrent_flag = True

                if recurrent_flag:
                    trycount += 1
                    continue
                else:
                    trycount = tries
                    found = True
                    break

        if found:
            # make sure recurrent 
            if recurrent:
                recurrent_flag = True
            
            # check is the innovation occured in the population
            new_innovation = True
            for innovation in innovations:
                if innovation.innovation_type == "link" and innovation.node_in_id == node1.nodeId and innovation.node_out_id == node2.nodeId and recurrent_flag == innovation.recur_flag:
                    # create a new gene with old innovation
                    newgene = Gene().add_gene_with_no_trait(weight=innovation.new_weight, inodeId=node1.nodeId, onodeId=node2.nodeId, recur_flag=recurrent_flag, innov=innovation.innovation_num1, mnum=0)
                    new_innovation = False
                    break
            
            if new_innovation:
                # choose a random trait
                    # TODO 
                # if len(self.traits) == 0:
                #     traitnum = 0        
                # else:
                #     traitnum = np.random.randint(0, len(self.traits))

                # choose a new weight
                newweight:float = np.random.normal(-1,1)
                # create new gene
                newgene = Gene().add_gene_with_no_trait(weight=newweight, inodeId=node1.nodeId, onodeId=node2.nodeId, recur_flag=recurrent_flag, innov=curinnov_num, mnum=0)
                
                # add to innovation
                newinnovation = Innovation().new_link(nin=node1.nodeId, nout=node2.nodeId, num1=curinnov_num, weight=newweight, recur=recurrent_flag)
                innovations.append(newinnovation)
            
            # add gene to genome
            self.add_gene(newgene)

            return True

        else:
            return False


    # nodeMutate: add node  # TODO delete node
    def mutate_add_node(self, innovations:List[Innovation], curinnov_num:int):
        # find a random gene
        thegene:Gene
        found:bool = False
        for tries in range(0,100):
            genenum = np.random.randint(0, len(self.genes))
            if self.genes[genenum].enabled and self.nodes[self.genes[genenum].inNodeId].nodeType != "bias":
                thegene = self.genes[genenum]
                found = True
                break

        if not found:
            return False

        thegene.enabled = False

        # check innovations for same innovation
        done = False
        for innovation in innovations:
            if innovation.innovation_type == "node" and innovation.node_in_id == thegene.inNodeId and innovation.node_out_id == thegene.outNodeId and innovation.old_innovation_num == thegene.innovation_num:
                # the same innovation
                inNodeId = thegene.inNodeId

                # create new node
                new_node = Node().new_node(nodeId=innovation.newnode_id, nodeType="hidden") ### is innovation.newnode_id unique?? is it the index of newnode in list?? no!!
                # change nodeId to match index
                new_node.nodeId = len(self.nodes)
                self.nodes.append(new_node)

                new_gene1 = Gene().add_gene_with_no_trait(weight=1.0, inodeId=thegene.inNodeId, onodeId=new_node.nodeId, recur_flag=thegene.is_recurrent, innov=innovation.innovation_num1, mnum=0)
                new_gene2 = Gene().add_gene_with_no_trait(weight=1.0, inodeId=new_node.nodeId, onodeId=thegene.outNodeId, recur_flag=False, innov=innovation.innovation_num2, mnum=0)
                self.genes.append(new_gene1)
                self.genes.append(new_gene2)
                done = True
                break
        
        if not done:
            # new innovation here
            new_node = Node().new_node(nodeId=len(self.nodes), nodeType="hidden")
            self.nodes.append(new_node)
            
            innov1 = curinnov_num
            innov2 = curinnov_num + 1
            new_gene1 = Gene().add_gene_with_no_trait(weight=1.0, inodeId=thegene.inNodeId, onodeId=new_node.nodeId, recur_flag=thegene.is_recurrent, innov=innov1, mnum=0)
            new_gene2 = Gene().add_gene_with_no_trait(weight=1.0, inodeId=new_node.nodeId, onodeId=thegene.outNodeId, recur_flag=False, innov=innov2, mnum=0)
            self.genes.append(new_gene1)
            self.genes.append(new_gene2)

            new_innov = Innovation().new_node(nin=thegene.inNodeId, nout=thegene.outNodeId, num1=innov1, num2=innov2, newid=new_node.nodeId, oldinnov=thegene.innovation_num)
            innovations.append(new_innov)
            done = True
        
        return True


        



    # add a gene to the genome
    def add_gene(self, newgene):
        self.genes.append(newgene)

    # iteration function to caculate if the link is recurrent
    def is_recurrent(self, innode, outnode, threshold, deepcount):
        deepcount += 1
        if deepcount > threshold:
            return False
        
        if innode == outnode:
            return True
        else:
            for geneId in innode.incoming:
                gene = self.genes[geneId]
                if not gene.is_recurrent:
                    if self.is_recurrent(self.nodes[gene.inNodeId], outnode, threshold, deepcount):
                        return True
            return False
        
    # copy this genome, return a new genome instance
    def copy(self):
        newgenome = Genome()
        newgenome.traits = self.traits.copy()
        for node in self.nodes:
            newnode = node.copy()
            newgenome.nodes.append(newnode)
        for gene in self.genes:
            newgene = gene.copy()
            newgenome.genes.append(newgene)
        for node in newgenome.nodes:
            if node.nodeType == "input":
                newgenome.inputs.append(node)
            elif node.nodeType == "hidden":
                newgenome.hiddens.append(node)
            elif node.nodeType == "output":
                newgenome.outputs.append(node)

        newgenome.fitness = self.fitness
        newgenome.orig_fitness = self.orig_fitness
        newgenome.winner = False
        newgenome.champion = False
        newgenome.expected_offspring = 0
        
        return newgenome

    # get last node id of genes
    def get_last_node_id(self):
        last_node:Node = self.ndoes[-1]
        return last_node.nodeId

    # get last gene's innovation number
    def get_last_gene_innovnum(self):
        last_gene:Gene = self.genes[-1]
        return last_gene.innovation_num

    def load_sensors(self, input:List[float]):
        i = 0
        for node in self.nodes:
            if node.nodeType == "input" and i <= 4:
                node.sensor_load(input[i])
                i += 1
            
    def activate(self):
        curnode:Node
        curlink:Gene

        add_amount:float # For adding to the activesum
        onetime:bool # Make sure we at least activate once
        abortcount:int = 0 # Used in case the output is somehow truncated from the network

        # Keep activating until all the outputs have become active 
        # (This only happens on the first activation, because after that they
        # are always active)

        onetime = False

        while(self.outputsoff() or not onetime):
            abortcount += 1

            if abortcount == 20:
                return False

            # For each node, compute the sum of its incoming activation 
            for curnode in self.nodes:
                if curnode.nodeType != "input":
                    curnode.activesum = 0
                    curnode.active_flag = False

                    # For each incoming connection, add the activity from the connection to the activesum 
                    for incomingId in curnode.incoming:
                        curlink = self.genes[incomingId]
                        if not curlink.time_delay:
                            add_amount = curlink.weight * self.nodes[curlink.inNodeId].get_active_out()
                            if self.nodes[curlink.inNodeId].active_flag or self.nodes[curlink.inNodeId].nodeType == "input":
                                curnode.active_flag =True
                            curnode.activesum += add_amount
                        else:
                            add_amount = curlink.weight * self.nodes[curlink.inNodeId].get_active_out_td()
                            curnode.activesum += add_amount
                        
            # Now activate all the non-sensor nodes off their incoming activation 
            for curnode in self.nodes:
                if curnode.nodeType != "input":
                    # Only activate if some active input came in
                    if curnode.active_flag:
                        # Keep a memory of activations for potential time delayed connections
                        curnode.last_activation2 = curnode.last_activation
                        curnode.last_activation = curnode.activation

                        # If the node is being overrided from outside,
                        # stick in the override value
                        if curnode.override:
                            curnode.active_override()
                        else:
                            # Now run the net activation through an activation function
                            slop = 1
                            curnode.activation = 1 / ( 1 + np.exp( -slop * curnode.activesum ))

                        curnode.activation_count += 1

            onetime = True

            # TODO adaptable
        return True

    # If not all output are active then return true
    def outputsoff(self):
        for curnode in self.nodes:
            if curnode.nodeType == "output" and curnode.activation_count == 0:
                return True
        
        return False

    # calculate the compatilbility of the variable genome and this genome
    def compatibility(self, genome):
            # Genes
            p1gene:Gene  # self
            p2gene:Gene  # outside

            # innovation numbers
            p1innov:int
            p2innov:int

            # intermediate value
            mut_diff:float

            # counters
            num_disjoint:float = 0
            num_excess:float = 0
            mut_diff_total:float = 0
            num_matching:float = 0

            max_genome_size:int = max(len(genome.genes), len(self.genes))
            # size of the larger genome

            # move through the genes in each potential parent
            # until both genome end
            p1 = 0
            p2 = 0
            while not p1 == len(self.genes) and not p2 == len(genome.genes):
                if p1 == len(self.genes):
                    p2 += 1
                    num_excess += 1
                elif p2 == len(self.genes):
                    p1 += 1
                    num_excess += 1
                else:
                    # extract current numbers
                    p1gene = self.genes[p1]
                    p2gene = genome.genes[p2]
                    p1innov = p1gene.innovation_num
                    p2innov = p2gene.innovation_num

                    if p1innov == p2innov:
                        num_matching += 1
                        mut_diff = abs(p1gene.mutation_num -p2gene.mutation_num)
                        mut_diff_total += mut_diff

                        p1 += 1
                        p2 += 1
                    if p1innov < p2innov:
                        p1 += 1
                        num_disjoint += 1
                    elif p1innov > p2innov:
                        p2+=1
                        num_disjoint += 1
            
            ret = Neat.disjoint_coeff * num_disjoint + Neat.excess_coeff * num_excess + Neat.mutdiff_coeff * mut_diff_total / num_matching
            return ret
    
    # mate a genome with this genome
    def mate_multipoint(self, g, genomeId:int, fitness1:float, fitness2:float, interspec_flag:bool):
        newtraits:List[Trait] # TODO
        newnodes:List[Node] = []
        newgenes:List[Gene] = []
        newgenome:Genome

        curgene2:Gene

        # iterators for moving through two parents genes
        p1gene:int = 0
        p2gene:int = 0
        p1innov:float
        p2innov:float
        chosengene:Gene
        inode:Node
        onode:Node
        disable:bool = False

        # decide which genome is better
        # TODO
        p1better:bool = True

        # move through each genome
        while not (p1gene == len(self.genes) and p2gene == len(g.genes)):
            skip:bool = False

            if p1gene == len(self.genes):
                chosengene = g.genes[p2gene]
                inode = g.nodes[chosengene.inNodeId]
                onode = g.nodes[chosengene.outNodeId]
                p2gene += 1
                if p1better:
                    skip = True
            elif p2gene == len(g.genes):
                chosengene = self.genes[p1gene]
                inode = self.nodes[chosengene.inNodeId]
                onode = self.nodes[chosengene.outNodeId]
                p1gene += 1
                if not p1better:
                    skip = True
            else:
                p1innov = self.genes[p1gene].innovation_num
                p2innov = g.genes[p2gene].innovation_num

                if p1innov == p2innov:
                    if np.random.uniform(0,1) < 0.5:
                        chosengene = self.genes[p1gene]
                        inode = self.nodes[chosengene.inNodeId]
                        onode = self.nodes[chosengene.outNodeId]
                    else:
                        chosengene = g.genes[p2gene]
                        inode = g.nodes[chosengene.inNodeId]
                        onode = g.nodes[chosengene.outNodeId]

                    # if one is disabled, offspring maybe disabled
                    if not self.genes[p1gene].enabled or not g.genes[p2gene].enabled:
                        if np.random.uniform(0,1) < 0.75:
                            disable =True
                    p1gene += 1
                    p2gene += 1
                elif p1innov < p2innov:
                    chosengene = self.genes[p1gene]
                    inode = self.nodes[chosengene.inNodeId]
                    onode = self.nodes[chosengene.outNodeId]
                    p1gene += 1
                    if p1better:
                        skip = True
                elif p1innov > p2innov:
                    chosengene = g.genes[p2gene]
                    inode = g.nodes[chosengene.inNodeId]
                    onode = g.nodes[chosengene.outNodeId]
                    p2gene += 1
                    if not p1better:
                        skip = True

            # for interspecies mating allow all genes grow through
            #TODO

            # check if chosengene conficts with an already chosen gene
            for curgene2 in newgenes:
                if curgene2.inNodeId == chosengene.inNodeId and curgene2.outNodeId == chosengene.outNodeId:
                    skip = True
                    break

            if not skip:
                # add chosengene to the baby

                # add trait #TODO

                # add nodes
                # check if inode and onode exists in newgenome
                inode_repeat = False
                for node in newnodes:
                    if node.nodeId_old == inode.nodeId and node.nodeType == inode.nodeType:
                        inode_repeat = True
                        gene_in_id = node.nodeId
                if not inode_repeat:
                    newnode1 = Node().new_node(inode.nodeId, inode.nodeType)
                    newnode1.nodeId_old = inode.nodeId
                    newnode1.nodeId = len(newnodes)
                    newnodes.append(newnode1)
                    gene_in_id = newnode1.nodeId
                # else:
                    # print("ignore one repeat in node, node Type: {}".format(inode.nodeType))
                
                
                onode_repeat = False
                for node in newnodes:
                    if node.nodeId_old == onode.nodeId and node.nodeType == inode.nodeType:
                        onode_repeat = True
                        gene_out_id = node.nodeId
                    
                if not onode_repeat:
                    newnode2 = Node().new_node(onode.nodeId, onode.nodeType)
                    newnode2.nodeId_old = onode.nodeId
                    newnode2.nodeId = len(newnodes)
                    newnodes.append(newnode2)
                    gene_out_id = newnode2.nodeId
                # else:
                    # print("ignore one repeat out node, node Type: {}".format(onode.nodeType))
                    
                
                # add gene
                newgene = Gene().add_gene_with_no_trait(weight=chosengene.weight, inodeId=gene_in_id, onodeId=gene_out_id, recur_flag=chosengene.is_recurrent, innov=chosengene.innovation_num, mnum=0)
                if disable:
                    newgene.enabled = False
                newgenes.append(newgene)
                newnodes[gene_in_id].outcoming.append(newgenes.index(newgene))
                newnodes[gene_out_id].incoming.append(newgenes.index(newgene))


        new_genome = Genome().create(newnodes=newnodes, newgenes=newgenes)

        # return baby genome
        return new_genome

    def create(self, newnodes:List[Node], newgenes:List[Gene]):
        newgenome = Genome()
        for node in newnodes:
            newnode = node.copy()
            newgenome.nodes.append(newnode)
        for gene in newgenes:
            newgene = gene.copy()
            newgenome.genes.append(newgene)
        for node in newgenome.nodes:
            if node.nodeType == "input":
                newgenome.inputs.append(node)
            elif node.nodeType == "hidden":
                newgenome.hiddens.append(node)
            elif node.nodeType == "output":
                newgenome.outputs.append(node)
        return newgenome


                





        

        


                        
                    


        







        
        
