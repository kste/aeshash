'''
Created on Feb 26, 2014

@author: Stefan Koelbl
'''

from gost.gost import GOST
from random import choice


import rebound
import copy
import pickle


def saveStatesToFile(filename, stateList):
    """Store a list of states in a file."""
    with open(filename, 'wb') as output:
        pickle.dump(stateList, output, pickle.HIGHEST_PROTOCOL)

def loadStatesFromFile(filename):
    """Loads a list of states from a file."""
    try:
        with open(filename, 'rb') as inputfile:
            result = pickle.load(inputfile)
    except:
        result = []
    return result

def checkForPossibleTransition(dM):
    """
    Checks if a transition from 1-8 is possible for all rows in L^1
    """
    
    # Get all byte pairs which are a possible input difference for SB^2
    validDiffsForSB3 = [[rebound.rebound.getValidDiffsForOutputDiff(dM[9].getValue(x, y)) for x in range(8)] for y in range(8)]

    # Need to find bytes for L^1 s.t 1-8 holds in every row
    matchedRows = 0
    for row in range(8):
        if(row != matchedRows):
            break
        for diff in rebound.rebound.LinearStepList:
            # Test for each byte
            if(not diff[0] in validDiffsForSB3[matchedRows][0]):
                continue
            if(not diff[1] in validDiffsForSB3[matchedRows][1]):
                continue
            if(not diff[2] in validDiffsForSB3[matchedRows][2]):
                continue
            if(not diff[3] in validDiffsForSB3[matchedRows][3]):
                continue
            if(not diff[4] in validDiffsForSB3[matchedRows][4]):
                continue
            if(not diff[5] in validDiffsForSB3[matchedRows][5]):
                continue
            if(not diff[6] in validDiffsForSB3[matchedRows][6]):
                continue
            if(not diff[7] in validDiffsForSB3[matchedRows][7]):
                continue
            # Found valid row
            dM[8].setRow(diff, matchedRows)  # K^2/L^1
            matchedRows += 1
            break

    return matchedRows



def randomStartDifferenceAndPropagteByte(dM):
    """
    Sets a random difference in AK^4_(0,0) and propagates it to SB^2
    """
    # Set random difference in AK^4
    dM[16].setRandomNonZeroValue(0, 0)

    # Propagate backwards to SB^2
    dM[15] = dM[16].Copy()  # same difference after key addition
    dM[14] = dM[15].Linverse()
    dM[13] = dM[14].Pinverse()

    # Select random differences through S-Box
    dM[12] = rebound.rebound.propagateDifferencesThroughSout(dM[13])

    # Propagate backwards to SB^2
    dM[11] = dM[12].Copy()  # same difference after key addition
    dM[10] = dM[11].Linverse()
    dM[9] = dM[10].Pinverse()

def randomStartDifferenceAndPropagte(dM):
    """
    Sets random differences for the first column in AK^3 and
    propagates it to SB^2.
    """
    # Set random difference in Column AK^3
    for row in range(8):
        dM[12].setRandomNonZeroValue(0, row)

    # Propagate backwards to SB^2
    dM[11] = dM[12].Copy()  # same difference after key addition
    dM[10] = dM[11].Linverse()
    dM[9] = dM[10].Pinverse()

    return

def solveConditionsOnSBox(m1, m2, dM, sboxToSolve):
    """
    Sets the correct byte values to fulfill the conditions imposed 
    by the differences.
    """
    # Set correct values for the state (random choice possible here)
    for x in range(8):
        for y in range(8):
            if(dM[sboxToSolve].getValue(x, y) != 0):
                # pick a random solution
                solution = choice(rebound.rebound.SBOX_MATCH[dM[sboxToSolve].getValue(x, y)][dM[sboxToSolve + 1].getValue(x, y)])

                m1[sboxToSolve].setValue(x, y, solution[0])
                m2[sboxToSolve].setValue(x, y, solution[1])
                m1[sboxToSolve + 1].setValue(x, y, GOST.sub[m1[sboxToSolve].getValue(x, y)])
                m2[sboxToSolve + 1].setValue(x, y, GOST.sub[m2[sboxToSolve].getValue(x, y)])
    return

def solveConditionsOnSBoxAndPropagateForward(m1, m2, dM, sboxToSolve):
    """
    Sets the correct byte values to fulfill the conditions imposed 
    by the differences and propagates the state forward by one round.
    """
    solveConditionsOnSBox(m1, m2, dM, sboxToSolve)

    # propagate forward
    m1[sboxToSolve + 1] = m1[sboxToSolve].S()
    m1[sboxToSolve + 2] = m1[sboxToSolve + 1].P()
    m1[sboxToSolve + 3] = m1[sboxToSolve + 2].L()

    m2[sboxToSolve + 1] = m2[sboxToSolve].S()
    m2[sboxToSolve + 2] = m2[sboxToSolve + 1].P()
    m2[sboxToSolve + 3] = m2[sboxToSolve + 2].L()

    return

def solveOneEightTransitionBySwap(m1, m2, k1, k2):
    """
    Solves the 1-8 transition by swapping the byte values (a, b) -> (b, a)
    and checking if it gives a possible 1-8 transition after propagating
    backward one round.
    """
    possibleDiffsFromSwap = [[] for _ in range(8)]
    swappedStates = [[] for _ in range(8)]

    for row in range(8):
        for byteToSwap in range(256):  # all possible swaps
            # Swap bytes in m1[11]/m2[11] and check if it gets the correct difference
            tmpM1 = copy.deepcopy(m1[7])
            tmpM2 = copy.deepcopy(m2[7])

            listToSwap = [int(x) for x in bin(byteToSwap)[2:].zfill(8)]  # creates a list of the binary representation to swap

            # swap bytes in state, this preserves differences
            for i in range(8):
                if(listToSwap[i] == 1):
                    tmp = tmpM1.getValue(i, row)
                    tmpM1.setValue(i, row, tmpM2.getValue(i, row))
                    tmpM2.setValue(i, row, tmp)

            # Store this row
            swappedStates[row].append([tmpM1.getRow(row), tmpM2.getRow(row)])

            # Propagate Backward
            tmpM1propagate = tmpM1.Linverse().P().Sinverse().AK(k1[8])
            tmpM2propagate = tmpM2.Linverse().P().Sinverse().AK(k2[8])
            
            # Store difference of first byte in list
            difference = tmpM1propagate.getValue(row, 0) ^ tmpM2propagate.getValue(row, 0)
            possibleDiffsFromSwap[row].append(difference)

    max_match = 0
    # Check for each possible 1-8 transition
    finished = False
    for diff in rebound.rebound.LinearStepList:
        # Solve for each byte in the column individually
        columnsMatched = 0

        while(columnsMatched != 4):
            # print possibleDiffsFromSwap[columnsMatched]
            if(not diff[columnsMatched] in possibleDiffsFromSwap[columnsMatched]):
                columnsMatched = 0
                break
            else:
                # also set state here
                indexOfRow = possibleDiffsFromSwap[columnsMatched].index(diff[columnsMatched])
                m1[7].setRow(swappedStates[columnsMatched][indexOfRow][0], columnsMatched)
                m2[7].setRow(swappedStates[columnsMatched][indexOfRow][1], columnsMatched)
                columnsMatched += 1
                if(columnsMatched > max_match):
                    max_match = columnsMatched
        if(columnsMatched == 4):
            print "Found solution for 7 rows, ignore last one"
            finished = True
            break
    return finished


def connectPathToOneEightTransition(dM):
    """
    Connect differences from SB^2 backward.
    """
    # Get all byte pairs which are a possible input difference for SB^2
    validDiffsForSB1 = [[rebound.rebound.getValidDiffsForOutputDiff(dM[5].getValue(x, y)) for x in range(8)] for y in range(8)]

    # Need to find bytes for L^1 s.t 1-8 holds in every row
    linearOneToEightTable = rebound.rebound.LinearStepList

    isPossible = True
    for diff in linearOneToEightTable:
        # Test for each byte
        isPossible = True
        for column in range(8):
            # Check if for all the bytes a valid pair exists
            if(not diff[column] in validDiffsForSB1[0][column]):
                isPossible = False
                break
        # All bytes are possible, use this difference
        if(isPossible):
            dM[4].setRow(diff, 0)
            print "Found Solution for AK^1 - S^1"
            return True

    return False

def solveConditionsWithLinearKeyLayerBackward(m1, m2, dM, k1, k2):
    """
    Use the key to solve conditions on S-Box between S^1 - S^2
    """
    # Determine values for the sbox
    solveConditionsOnSBox(m1, m2, dM, 4)
    solveConditionsOnSBox(m1, m2, dM, 8)

    targetM1 = m1[5].Copy()

    # Values in one row are determined now
    # We can now solve 7 rows here using freedom of K^2
    for row in range(7):
        for value in range(256):
            k1[12].setValue(0, row, value)

            #propagate backward
            propagateBackwardWithKey(m1, k1, 4, 9)

            #check if byte is correct
            if(targetM1.getValue(row, 0) == m1[5].getValue(row, 0)):
                k2[12] = k1[12].Copy()

                propagateForwardKey(k1, 12, 20)
                propagateForwardWithKey(m1, k1, 8, 16)

                propagateBackwardKey(k1, 0, 13)
                propagateBackwardWithKey(m1, k1, 0, 9)

                propagateForwardKey(k2, 12, 20)
                propagateBackwardKey(k2, 0, 13)
                propagateBackwardWithKey(m2, k2, 0, 9)
                propagateForwardWithKey(m2, k2, 8, 16)
                break
    return


def solveConditionsWithLinearKeyLayerForward(m1, m2, dM, k1, k2):
    """
    Use the key to solve conditions on S-Box between S^2 - S^3
    """
    
    # Solve for each row independently
    # Keys should be all zero at this point
    solveConditionsOnSBox(m1, m2, dM, 12)

    targetM1 = m1[12].Copy()

    for row in range(8):
        #while(True):
        #    k1[15].setRandomValue(7, row)
        for value in range(256):
            k1[15].setValue(7, row, value) #test values

            propagateForwardKey(k1, 15, 16)
            propagateForwardWithKey(m1, k1, 8, 13)

            if(targetM1.getValue(0, row) == m1[12].getValue(0, row)):
                k2[15] = k1[15].Copy()

                propagateBackwardKey(k1, 0, 16)
                propagateBackwardWithKey(m1, k1, 0, 9)
                propagateForwardWithKey(m1, k1, 8, 16)

                propagateForwardKey(k2, 15, 20)
                propagateBackwardKey(k2, 0, 16)
                #propagate backward
                propagateBackwardWithKey(m2, k2, 0, 9)
                propagateForwardWithKey(m2, k2, 8, 16)
                break
    return

def propagateForwardKey(stateList, fromRound, toRound):
    """
    Apply round functions to the key.
    """
    for i in range(fromRound, toRound):
        if(i % 4 == 0):
            stateList[i + 1] = stateList[i].AddConstant((i / 4) - 1)
        if(i % 4 == 1):
            stateList[i + 1] = stateList[i].S()
        if(i % 4 == 2):
            stateList[i + 1] = stateList[i].P()
        if(i % 4 == 3):
            stateList[i + 1] = stateList[i].L()

def propagateBackwardKey(stateList, fromRound, toRound):
    """
    Apply inverse round functions to the key.
    """    
    for i in reversed(range(fromRound, toRound)):
        if(i % 4 == 0):
            stateList[i - 1] = stateList[i].Linverse()
        if(i % 4 == 3):
            stateList[i - 1] = stateList[i].Pinverse()
        if(i % 4 == 2):
            stateList[i - 1] = stateList[i].Sinverse()
        if(i % 4 == 1):
            stateList[i - 1] = stateList[i].AddConstant((i / 4) - 1)

def propagateForward(stateList, fromRound, toRound):
    """
    Apply round functions to the state without key addition.
    """    
    for i in range(fromRound, toRound):
        if(i % 4 == 0):
            stateList[i + 1] = stateList[i].S()
        if(i % 4 == 1):
            stateList[i + 1] = stateList[i].P()
        if(i % 4 == 2):
            stateList[i + 1] = stateList[i].L()
        if(i % 4 == 3):
            stateList[i + 1] = stateList[i].Copy()

def propagateBackward(stateList, fromRound, toRound):
    """
    Apply inverse round functions to the state without key addition.
    """    
    for i in reversed(range(fromRound, toRound)):
        if(i % 4 == 0):
            stateList[i - 1] = stateList[i].Copy()
        if(i % 4 == 3):
            stateList[i - 1] = stateList[i].Linverse()
        if(i % 4 == 2):
            stateList[i - 1] = stateList[i].Pinverse()
        if(i % 4 == 1):
            stateList[i - 1] = stateList[i].Sinverse()

def propagateForwardWithKey(stateList, keyList, fromRound, toRound):
    """
    Apply round functions to the state including key addition.
    """    
    for i in range(fromRound, toRound):
        if(i % 4 == 0):
            stateList[i + 1] = stateList[i].S()
        if(i % 4 == 1):
            stateList[i + 1] = stateList[i].P()
        if(i % 4 == 2):
            stateList[i + 1] = stateList[i].L()
        if(i % 4 == 3):
            stateList[i + 1] = stateList[i].AK(keyList[i + 5])

def propagateBackwardWithKey(stateList, keyList, fromRound, toRound):
    """
    Apply inverse round function to the state including key addition.
    """
    for i in reversed(range(fromRound, toRound)):
        if(i % 4 == 0):
            stateList[i - 1] = stateList[i].AK(keyList[i + 4])
        if(i % 4 == 3):
            stateList[i - 1] = stateList[i].Linverse()
        if(i % 4 == 2):
            stateList[i - 1] = stateList[i].Pinverse()
        if(i % 4 == 1):
            stateList[i - 1] = stateList[i].Sinverse()

def startFromTheMiddle(filename):
    """
    Finds a collision for 4 rounds by combining the start from the middle and
    our message modification technique.
    """
    # Initialize States used
    rounds = 4
    k1 = [GOST() for _ in range(rounds * 5 + 1)]
    k2 = [GOST() for _ in range(rounds * 5 + 1)]
    dK = [GOST() for _ in range(rounds * 5 + 1)]
    m1 = [GOST() for _ in range(rounds * 4 + 1)]
    m2 = [GOST() for _ in range(rounds * 4 + 1)]
    dM = [GOST() for _ in range(rounds * 4 + 1)]

    iterations = 0
    counter = 0
    while(True):
        maxTransitions = 0
        while(True):
            iterations += 1
            foundSolution = True

            randomStartDifferenceAndPropagteByte(dM)

            # Test for all possible 1-8 differences for L
            transitions = checkForPossibleTransition(dM)
            if(transitions > maxTransitions):
                maxTransitions = transitions
                print "Found transition for " + str(maxTransitions)
            if(not transitions == 8):
                #no transition found for all 8 rows
                foundSolution = False

            if(foundSolution):
                print "Phase 1 Completed"
                break

        dM[7] = dM[8] #only valid if there are no diffs in key
        propagateBackward(dM, 6, 8)

        # Now we need to solve the 1-8 transition by using the freedom of the byte pairs in L^1 = m1[11]
        # This can be done for each row individually and we again need to find a match
        # with the linear layer for a 1-8 transition
        if(connectPathToOneEightTransition(dM)):
            print "Phase 2 Completed"
            propagateBackward(dM, 1, 5)
            break

    # Forward Propagation
    propagateForwardWithKey(m1, k1, 0, 16)
    propagateForwardWithKey(m2, k2, 0, 16)


    #print states
    labelsKey = ["K^0", "S^0", "P^0", "L^0", "K^1", "S^1", "P^1", "L^1", "K^2", "S^2", "P^2", "L^2", "K^3", "S^3", "P^3", "L^3", "K^4"]
    labelsM   = ["AK^0", "S^0", "P^0", "L^0", "AK^1", "S^1", "P^1", "L^1", "AK^2", "S^2", "P^2", "L^2", "AK^3", "S^3", "P^3", "L^3", "AK^4"]
    for i in range(17):
        print "------------------------" + labelsKey[i] + "---------------------------------------------" + labelsM[i] + "------------------------------------"
        dM[0].printStatesHex([k1[i+4], k2[i+4], dK[i+4], m1[i], m2[i], m1[i].AK(m2[i]), dM[i]])

    saveStatesToFile(filename, [k1, k2, dK, m1, m2, dM])

    foundValidPaths = 0
    while(foundValidPaths < 300):
        # Now solve conditions AK^3 - S^3 with key freedom
        # This is done by changing bytes in last column in L^-1(K^3)
        solveConditionsWithLinearKeyLayerBackward(m1, m2, dM, k1, k2)
        solveConditionsWithLinearKeyLayerForward(m1, m2, dM, k1, k2)

        # Now Check if condition at start is fulfilled by bruteforce
        counter += 1

        if(m1[0].AK(m2[0]).NumberOfActiveBytes() == 1):
            foundValidPaths += 1
            #check if differences cancel out
            if(m1[0].AK(m2[0]).getValue(0, 0) == m1[16].AK(m2[16]).getValue(0, 0)):
                print "Found solution " + str(foundValidPaths)
                break
            else:
                print "Found 1-8-64-8-1 but bytes do not match " + str(m1[0].AK(m2[0]).getValue(0, 0)) + " != " + str(m1[16].AK(m2[16]).getValue(0, 0)) + " Paths Tested: " + str(foundValidPaths)
        else:
            print "No solution Found " + str(counter)

    if(foundValidPaths >= 300):
        # restart
        print "Restarted Attack"
        startFromTheMiddle(filename)
        return

    # Forward Propagation
    propagateForwardWithKey(m1, k1, 0, 16)
    propagateForwardWithKey(m2, k2, 0, 16)


    #print states
    labelsKey = ["K^0", "S^0", "P^0", "L^0", "K^1", "S^1", "P^1", "L^1", "K^2", "S^2", "P^2", "L^2", "K^3", "S^3", "P^3", "L^3", "K^4"]
    labelsM   = ["AK^0", "S^0", "P^0", "L^0", "AK^1", "S^1", "P^1", "L^1", "AK^2", "S^2", "P^2", "L^2", "AK^3", "S^3", "P^3", "L^3", "AK^4"]
    for i in range(17):
        print "------------------------" + labelsKey[i] + "---------------------------------------------" + labelsM[i] + "------------------------------------"
        dM[0].printStatesHex([k1[i+4], k2[i+4], dK[i+4], m1[i], m2[i], m1[i].AK(m2[i]), dM[i]])

    return

def findMessagePairFromLoadedState(filename):
    """
    Finds the message pair for a previously stored result of Phase 2.
    """
    [k1, k2, dK, m1, m2, dM] = loadStatesFromFile(filename)

    counter = 0
    while(True):
        # Now solve conditions AK^3 - S^3 with key freedom
        # This is done by changing bytes in last column in L^-1(K^3)
        solveConditionsWithLinearKeyLayerBackward(m1, m2, dM, k1, k2)
        solveConditionsWithLinearKeyLayerForward(m1, m2, dM, k1, k2)

        # Now Check if condition at start is fulfilled by bruteforce
        counter += 1

        if(m1[0].AK(m2[0]).NumberOfActiveBytes() == 1):
            #check if differences cancel out
            if(m1[0].AK(m2[0]).getValue(0, 0) == m1[16].AK(m2[16]).getValue(0, 0)):
                print "Found solution"
                break
            else:
                print "Found 1-8-64-8-1 but bytes do not match " + str(m1[0].AK(m2[0]).getValue(0, 0)) + " != " + str(m1[16].AK(m2[16]).getValue(0, 0))
        else:
            print "No solution Found " + str(counter)


    propagateBackwardWithKey(m1, k1, 0, 8)
    propagateBackwardWithKey(m2, k2, 0, 8)
    #print states
    labelsKey = ["K^0", "S^0", "P^0", "L^0", "K^1", "S^1", "P^1", "L^1", "K^2", "S^2", "P^2", "L^2", "K^3", "S^3", "P^3", "L^3", "K^4"]
    labelsM   = ["AK^0", "S^0", "P^0", "L^0", "AK^1", "S^1", "P^1", "L^1", "AK^2", "S^2", "P^2", "L^2", "AK^3", "S^3", "P^3", "L^3", "AK^4"]
    for i in range(17):
        print "------------------------" + labelsKey[i] + "---------------------------------------------" + labelsM[i] + "------------------------------------"
        dM[0].printStatesHex([k1[i+4], k2[i+4], dK[i+4], m1[i], m2[i], m1[i].AK(m2[i]), dM[i]])

    # Forward Propagation
    propagateForwardWithKey(m1, k1, 0, 16)
    propagateForwardWithKey(m2, k2, 0, 16)


    #print states
    labelsKey = ["K^0", "S^0", "P^0", "L^0", "K^1", "S^1", "P^1", "L^1", "K^2", "S^2", "P^2", "L^2", "K^3", "S^3", "P^3", "L^3", "K^4"]
    labelsM   = ["AK^0", "S^0", "P^0", "L^0", "AK^1", "S^1", "P^1", "L^1", "AK^2", "S^2", "P^2", "L^2", "AK^3", "S^3", "P^3", "L^3", "AK^4"]
    for i in range(17):
        print "------------------------" + labelsKey[i] + "---------------------------------------------" + labelsM[i] + "------------------------------------"
        dM[0].printStatesHex([k1[i+4], k2[i+4], dK[i+4], m1[i], m2[i], m1[i].AK(m2[i]), dM[i]])

    return

def printSavedStates(filename):
    """
    Prints the states from a previously stored file.
    """

    [k1, k2, dK, m1, m2, dM] = loadStatesFromFile(filename)

    #print states
    labelsKey = ["K^0", "S^0", "P^0", "L^0", "K^1", "S^1", "P^1", "L^1", "K^2", "S^2", "P^2", "L^2", "K^3", "S^3", "P^3", "L^3", "K^4"]
    labelsM   = ["AK^0", "S^0", "P^0", "L^0", "AK^1", "S^1", "P^1", "L^1", "AK^2", "S^2", "P^2", "L^2", "AK^3", "S^3", "P^3", "L^3", "AK^4"]
    for i in range(17):
        print "------------------------" + labelsKey[i] + "---------------------------------------------" + labelsM[i] + "------------------------------------"
        dM[0].printStatesHex([k1[i+4], k2[i+4], dK[i+4], m1[i], m2[i], m1[i].AK(m2[i]), dM[i]])

    # Forward Propagation
    propagateForwardWithKey(m1, k1, 0, 16)
    propagateForwardWithKey(m2, k2, 0, 16)

    #print states
    labelsKey = ["K^0", "S^0", "P^0", "L^0", "K^1", "S^1", "P^1", "L^1", "K^2", "S^2", "P^2", "L^2", "K^3", "S^3", "P^3", "L^3", "K^4"]
    labelsM   = ["AK^0", "S^0", "P^0", "L^0", "AK^1", "S^1", "P^1", "L^1", "AK^2", "S^2", "P^2", "L^2", "AK^3", "S^3", "P^3", "L^3", "AK^4"]
    for i in range(17):
        print "------------------------" + labelsKey[i] + "---------------------------------------------" + labelsM[i] + "------------------------------------"
        dM[0].printStatesHex([k1[i+4], k2[i+4], dK[i+4], m1[i], m2[i], m1[i].AK(m2[i]), dM[i]])

    return

