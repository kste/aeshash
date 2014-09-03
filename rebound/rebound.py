'''
Created on Jan 9, 2014

@author: Stefan Koelbl
'''


from random import sample
from gost.gost import GOST

DDT = []
SBOX_MATCH = []
LinearStepList = []
InverseLinearStepList = []
possibleOutputDifferences = []
possibleInputDifferences = []

def propagateDifferencesThroughSout(state):
    """
    Returns a state containing a possible output difference after applying
    the S-Box for the given input difference.
    """
    result = GOST()
    for x in range(8):
        for y in range(8):
            result.setValue(x, y, sample(getValidDiffsForOutputDiff(state.getValue(x, y)), 1)[0])
    return result

def propagateDifferencesThroughSin(state):
    """
    Returns a state containing a possible input difference after applying
    the inverse S-Box for the given output difference.
    """
    result = GOST()
    for x in range(8):
        for y in range(8):
            result.setValue(x, y, sample(getValidDiffsForInputDiff(state.getValue(x, y)), 1)[0])
    return result

def computeLinearStepList():
    """
    Compute the list of all possible values for 
    (x 0 0 0 0 0 0 0) * L = (y0 y1 y2 y3 y4 y5 y6 y7)
    """
    global LinearStepList
    gost = GOST()
    for value in range(1, 256):
        gost.setValue(0, 0, value)
        LinearStepList.append(gost.L().getRow(0))

def computeInverseLinearStepList():
    """
    Compute the list of all possible values for 
    (x 0 0 0 0 0 0 0) * Linverse = (y0 y1 y2 y3 y4 y5 y6 y7)
    """    
    global InverseLinearStepList
    gost = GOST()
    for value in range(1, 256):
        gost.setValue(0, 0, value)
        InverseLinearStepList.append(gost.Linverse().getRow(0))

def computeDDT(sbox):
    """
    Compute the differential distribution table (DDT) for a given S-Box
    """
    global DDT
    DDT = [[0 for _ in range(len(sbox))] for _ in range(len(sbox))]
    for a in range(len(sbox)):
        for b in range(len(sbox)):
            DDT[a ^ b][sbox[a] ^ sbox[b]] += 1

def computeSBOX_MATCH(sbox):
    """
    Compute the valid pairs for each input/output difference.
    """
    global SBOX_MATCH
    SBOX_MATCH = [[[] for _ in range(len(sbox))] for _ in range(len(sbox))]
    for a in range(len(sbox)):
        for b in range(len(sbox)):
            SBOX_MATCH[a ^ b][sbox[a] ^ sbox[b]].append([a, b])

def getValidBytePairsForOutputDiff(outputDiff):
    """
    Get all possible pairs (a, b) such that:
    S(a) xor S(b) = outputDiff
    """
    bytePairs = []
    for i in range(len(SBOX_MATCH)):
        if(len(SBOX_MATCH[i][outputDiff]) > 0):
            bytePairs.append(SBOX_MATCH[i][outputDiff])
    return bytePairs

def getValidBytePairsForInputDiff(inputDiff):
    """
    Get all possible pairs (a, b) such that:
    Sinverse(a) xor Sinverse(b) = inputDiff
    """
    bytePairs = []
    for i in range(len(SBOX_MATCH)):
        if(len(SBOX_MATCH[inputDiff][i]) > 0):
            bytePairs.append(SBOX_MATCH[inputDiff][i])
    return bytePairs

def getValidDiffsForInputDiff(inputDiff):
    """
    Get all possible output differences for a given input difference.
    """
    global possibleOutputDifferences
    if not possibleOutputDifferences:
        possibleOutputDifferences = [set([]) for _ in range(256)]
        # Compute Table
        for diffIn in range(256):
            for diffOut in range(256):
                if(DDT[diffIn][diffOut] > 0):
                    possibleOutputDifferences[diffIn].add(diffOut)
        
    return possibleOutputDifferences[inputDiff]


def getValidDiffsForOutputDiff(outputDiff):
    """
    Get all possible input differences for a given output difference.
    """
    global possibleInputDifferences
    if not possibleInputDifferences:
        possibleInputDifferences = [set([]) for _ in range(256)]
        # Compute Table
        for diffIn in range(256):
            for diffOut in range(256):
                if(DDT[diffIn][diffOut] > 0):
                    possibleInputDifferences[diffOut].add(diffIn)
                    
    return possibleInputDifferences[outputDiff]