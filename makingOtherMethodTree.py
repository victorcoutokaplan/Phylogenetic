from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo, AlignIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd
import dendropy




def inferTreeFromSeqs(sequenceDict, mu_Var):
    writeAlignObjectFromSeq(sequenceDict)
    align = AlignIO.read('sequencesFile.phy','phylip')

    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    distMatrix = calculator.get_distance(align)
    print(distMatrix)

    # # Create a DistanceTreeConstructor object
    # constructor = DistanceTreeConstructor()
    # # Construct the phlyogenetic tree using UPGMA algorithm
    # UGMATree = constructor.upgma(distMatrix)
    # # Construct the phlyogenetic tree using NJ algorithm
    # NJTree = constructor.nj(distMatrix)
    #


    numSeqs = len(sequenceDict.keys())
    fullMatrix = np.zeros((numSeqs,numSeqs))
    for i in range(numSeqs):
        for j in range(i):
            fullMatrix[i][j] = distMatrix.matrix[i][j]
            fullMatrix[j][i] = fullMatrix[i][j]

    df = pd.DataFrame(fullMatrix)
    df.to_csv('otherMethodDistanceMatrix.csv')

    pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
        src=open('otherMethodDistanceMatrix.csv'),
        delimiter=",")
    nj_tree = pdm.nj_tree()

    for leaf in nj_tree.leaf_node_iter():
        leaf.taxon.label = list(sequenceDict.keys())[int(leaf.taxon.label)]

    for edge in nj_tree.postorder_edge_iter():
        if edge.length != None:
            edge.length /= mu_Var

    return nj_tree


def writeAlignObjectFromSeq(seqs):
    toWrite = []
    for seqIndx in seqs:
        toWrite.append(SeqRecord(Seq(seqs[seqIndx]),  id=str(seqIndx), name=str(seqIndx)))
    align = Align.MultipleSeqAlignment(toWrite)
    AlignIO.write(align, 'sequencesFile.phy', 'phylip')
    return None
