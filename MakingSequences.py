import collections
import dendropy
import math
import numpy as np
import msprime
import MakingDendroPyTrees


def makeSequencesFromScratch(birth_rate, death_rate, n_populations, extra_populations, seqLenInput, initial_size_Var, mu_Var, recombination_rate):
    tree = MakingDendroPyTrees.makeDendroTree(birth_rate, death_rate, n_populations, extra_populations)
    for edge in tree.postorder_edge_iter():
        edge.length /= mu_Var
    return makeSequencesFromTree(tree, seqLenInput, initial_size_Var, mu_Var, recombination_rate), tree

def makeSequencesFromTree(tree, seqLenInput, initial_size_Var, mu_Var, recombination_rate):
    rawSeqs = getRawSequencesFromTree(tree, seqLenInput, initial_size_Var, mu_Var, recombination_rate)
    return fixedUp(rawSeqs, seqLenInput)

def getRawSequencesFromTree(tree, seqLenInput, initial_size_Var, mu_Var, recombination_rate):

    demography = msprime.Demography()
    real_dem = demography.from_species_tree(str(tree), initial_size_Var)
    print(real_dem)

    TtoIdMap = {}
    for pop in real_dem.populations:
        if pop.name[0] != 'p':
            TtoIdMap[pop.id] = pop.name

    samplesDict = {}
    for i in range(len(tree.leaf_nodes())):
        # keyStr = 'T' + str(i+1)
        samplesDict[i] = 1
        # samplesDict[keyStr] = 1

    ts = msprime.sim_ancestry(samples=samplesDict,
                              sequence_length=seqLenInput,
                              recombination_rate=recombination_rate,
                              demography=real_dem,
                              ploidy=1)
    mts = msprime.sim_mutations(ts, rate=mu_Var, model='JC69')

    alignments = {}
    i = 0
    for t in mts.alignments():
        alignments[TtoIdMap[i]] = t
        i += 1

    return alignments

def fixedUp(alignments, seqLenInput):
    """
    Takes alignments with 'N' for nonmutated spots and turns them into strings.
    First makes a random ancestral alignment, and then fills in the 'N' spots
    IN EACH SEQUENCE of alignments with the corresponding ancestral alignment spot
    """

    digitToGene = {0:'A', 1:'T', 2:'C', 3:'G'}
    digitSeq = np.random.randint(0,4,int(seqLenInput))
    randomSeq = [digitToGene[digit] for digit in digitSeq]

    newAlignments = {}

    for alignment in alignments:
        newSequence = ''
        for letterInx in range(len(alignments[alignment])):
            if alignments[alignment][letterInx] == 'N':
                newSequence += randomSeq[letterInx]
            elif alignments[alignment][letterInx] != 'N':
                newSequence += alignments[alignment][letterInx]
        newAlignments[alignment] = newSequence
    return newAlignments
