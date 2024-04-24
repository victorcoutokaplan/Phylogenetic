import dendropy
import numpy as np
import pandas as pd

def makeTree(DistMatrix):
    df = pd.DataFrame(DistMatrix)
    df.to_csv('ourMethodDistanceMatrix.csv')

    pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
        src=open('ourMethodDistanceMatrix.csv'),
        delimiter=",")

    nj_tree = pdm.nj_tree()
    return nj_tree
