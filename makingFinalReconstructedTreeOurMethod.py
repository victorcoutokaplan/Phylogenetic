import dendropy
import numpy as np
import pandas as pd
import makingPairwiseDistancesOurMethod


def makeManyMatricesOfDistances(number_cuts, sequenceDict):

    lenSeq = len(sequenceDict['T1'])
    mats = []

    for i in range(number_cuts):
        cutSeqs = {}
        for key in sequenceDict.keys():
            cutSeqs[key] = sequenceDict[key][int(i * lenSeq/number_cuts) : int((i+1)*lenSeq/number_cuts)]
        ourMatrix = makingPairwiseDistancesOurMethod.makeTimeAndSizeMatrix(cutSeqs, lenSeq/number_cuts)
        mats.append(ourMatrix)

    return mats

def makeFinalDistMatrixFromSmallerOnes(mats):

    countMat = np.zeros([5,5])

    for mat in mats:
        distMat = mat[0]
        for i in range(len(mats[0][0][0])):
            for j in range(len(mats[0][0][0])):
                if distMat[i][j] > -0.05:
                    countMat[i][j] += 1

    print(countMat)

    finalDistMat = np.zeros([5,5])

    for mat in mats:
        distMat = mat[0]
        for i in range(len(mats[0][0][0])):
            for j in range(len(mats[0][0][0])):
                if distMat[i][j] > 0:
                    finalDistMat[i][j] += distMat[i][j] / countMat[i][j]

    df = pd.DataFrame(finalDistMat)

    for i in range(len(mats[0][0][0])):
        for j in range(i+1):
            if i != j:
                if finalDistMat[i][j] == 0.0:
                    print("fail")

    return df


def makeTree(sequenceDict, mu_Var, DistMatrix):
    df = pd.DataFrame(DistMatrix)
    df.to_csv('ourMethodDistanceMatrix.csv')

    pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
        src=open('ourMethodDistanceMatrix.csv'),
        delimiter=",")

    nj_tree = pdm.nj_tree()

    for leaf in nj_tree.leaf_node_iter():
        leaf.taxon.label = list(sequenceDict.keys())[int(leaf.taxon.label)]

    for edge in nj_tree.postorder_edge_iter():
        if edge.length != None:
            edge.length /= mu_Var



    return nj_tree
