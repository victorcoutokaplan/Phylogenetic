import numpy as np
import scipy.optimize
import scipy.special
import math

def makeTimeAndSizeMatrix(dictionaryOfSequences, sequenceLengths):
    '''
    Takes in a dictionary of sequences. Pairwise computes estimated t and theta.
    Returns a matrix of the (t, theta) infered pairs, for each spot i,j where
    the infered values are from comparing dictionaryOfSequences[i] and
    dictionaryOfSequences[j]
    '''

    dictionaryOfSequencesRelabeled = {}

    for intIndex in range(len(dictionaryOfSequences.keys())):
        dictionaryOfSequencesRelabeled[intIndex] = dictionaryOfSequences[list(dictionaryOfSequences.keys())[intIndex]]

    numberLeafNodes = len(dictionaryOfSequencesRelabeled.keys())
    finalMatrixTime = np.zeros((numberLeafNodes,numberLeafNodes))
    finalMatrixPop = np.zeros((numberLeafNodes,numberLeafNodes))

    # get distances
    # distList = [7,8,9,10,11,12]
    distList = [9,10,11,12]
    distances = computeDistancesForKCounts(distList, dictionaryOfSequencesRelabeled)

    for seqKey1 in range(numberLeafNodes):
        for seqKey2 in range(seqKey1 + 1, numberLeafNodes):
            print('working on inferring', seqKey1, seqKey2)
            #get curves
            curvesDict = findCurves(distList, distances, (seqKey1,seqKey2), sequenceLengths)

            inferred = minimizeWithNewVersion_10_11_12_withFailGoingTo11_12Then9_10_11_12(curvesDict, max(np.random.normal(loc=0.002, scale=0.0005), 0.0002) ) # random starting choice
            finalMatrixTime[(seqKey1, seqKey2)] = inferred[1]
            finalMatrixTime[(seqKey2, seqKey1)] = inferred[1]
            finalMatrixPop[(seqKey1, seqKey2)] = inferred[0]
            finalMatrixPop[(seqKey2, seqKey1)] = inferred[0]
    return finalMatrixTime, finalMatrixPop

#=============================================
# GET DISTANCES

def computeDistancesForKCounts(kcounts, sequenceDictionary):
    numberLeafNodes = len(sequenceDictionary.keys())
    distances = {}

    calcPossibleDict = {}
    for k in kcounts:
        calcPossibleDict[k] = calcPossible(k)

    kmerVectDict = {}
    for k in kcounts:
        for seqKey in range(numberLeafNodes):
            print('in getting kmer vectors, working on k:', k, 'seq', seqKey)
            kmerVectDict[(k, seqKey)] = getKMerVect(sequenceDictionary[seqKey], k, calcPossibleDict)

    for k in kcounts:
        distancesForK = np.zeros((numberLeafNodes, numberLeafNodes))
        for seqKey1 in range(numberLeafNodes):
            for seqKey2 in range(seqKey1 + 1, numberLeafNodes):
                s1v = kmerVectDict[(k, seqKey1)]
                s2v = kmerVectDict[(k, seqKey2)]
                distancesForK[(seqKey1, seqKey2)] = np.linalg.norm(s1v - s2v) ** 2
                distancesForK[(seqKey2, seqKey1)] = distancesForK[(seqKey1, seqKey2)]

        distances[k] = distancesForK
    return distances

def getKMerVect (seq, k, calcPossibleDict):
    ngrams = find_ngrams(seq, k)
    possibilities = calcPossibleDict[k]
    ngrams += possibilities
    blownUpVect = np.unique(ngrams , return_counts = True)[1]
    blownUpVect -= np.full(4**k, 1)
    return blownUpVect

def find_ngrams(input_list, n):
    ngrams = zip(*[input_list[i:] for i in range(n)])
    return [''.join(ngram) for ngram in ngrams]

def calcPossible(k):
    if (k==0):
        return ['']
    else:
        last = calcPossible(k-1)
        a = ['A' + entry for entry in last]
        c = ['C' + entry for entry in last]
        t = ['T' + entry for entry in last]
        g = ['G' + entry for entry in last]
        return a + c + t + g


#=============================================
# GET CURVES

def findCurves(distList, distances, LeafIds, sequenceLength):
    result = {}

    for kNum in distList:
        ridgePoints = findRidge(kNum, distances[kNum][LeafIds] , sequenceLength)

        coeffs = np.polyfit(ridgePoints[0],ridgePoints[1], 15) # variable 15, should not be

        result[kNum] = np.polynomial.Polynomial(np.flip(coeffs))

    scaledResult = {}

    for key in result:
        scaledResult[key] = 100*result[key] # scaled still

    return scaledResult


def findRidge(k, commonDist, sequenceLength):
    final = []
    for thetaValue in np.linspace(0.0000001,.01,5000): #note bound set here should be variable
        final.append(    findCorrectTFromTheta(k, commonDist, thetaValue, sequenceLength)   )

    return np.array([np.linspace(0.0000001,.01,5000), np.array(final)]) #note bound set here should be variable

def findCorrectTFromTheta(k, commonDist, thetaValue, sequenceLength):

    curr_f = lambda t: expectedDist(k, t, thetaValue, sequenceLength) - commonDist

    return scipy.optimize.fsolve(curr_f, x0=thetaValue)[0]

#this is our expected dist function:
def expectedDist(k, t, theta, m):
    firstPart = 2*(m-k+1)
    param = 3/(8*theta)
    summation = 1 - (1/(4**k)) * scipy.special.hyp2f1(-k, param, 1 + param, -3*math.exp((-4*t)/3))
    return firstPart * summation


#=============================================
# INFERENCE

def minimizeWithNewVersion_10_11_12_withFailGoingTo11_12Then9_10_11_12(curves, startingTheta):
    '''
    current state. Note it MUST take in curves for k = 9,10,11,12 so only four.
    '''

    #trying 10_11_12
    curves_to_try = [curves[10], curves[11], curves[12]]
    candidateResult = minimizeWithNewVersion(curves_to_try, startingTheta)
    if candidateResult[0] > 0.0001 and  candidateResult[1] > 0.0001:
        return candidateResult

    #trying 11,12
    curves_to_try = [curves[11], curves[12]]
    candidateResult = minimizeWithNewVersion(curves_to_try, startingTheta)
    if candidateResult[0] > 0.0001 and  candidateResult[1] > 0.0001:
        return candidateResult

    #trying 9_10_11_12
    curves_to_try = [curves[9], curves[10], curves[11], curves[12]]
    candidateResult = minimizeWithNewVersion(curves_to_try, startingTheta)
    if candidateResult[0] > 0.0001 and  candidateResult[1] > 0.0001:
        return candidateResult

    # none converged
    return -1,-1

def forThetaFindMinCurvesTVictorVersion(theta, curves):
    def bayesianDist(candidateT, curves):
        exponent = (sum([(curve(theta) - candidateT)**2 for curve in curves]))
        return -np.exp(-exponent)

    newGuess = scipy.optimize.minimize(bayesianDist, max(np.random.normal(loc=0.004, scale=0.0005), 0.0004) , args=curves)['fun'] #RANDOM STARTING GUESS
    return np.sqrt(np.sqrt(newGuess+1))

def minimizeWithNewVersion(curves, startingTheta):
    finalTheta = scipy.optimize.minimize(forThetaFindMinCurvesTVictorVersion, startingTheta, args=(curves))
    aveCurve = sum(curves)/len(curves)
    finalTheta = finalTheta['x'][0]
    return finalTheta, aveCurve(finalTheta)/100

#=============================================
