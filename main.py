'''
Created on Jan 7, 2014

@author: Stefan Koelbl
'''

from gost.gost import GOST
from attacks import startinthemiddle

import rebound
import sys
import os

if __name__ == '__main__':

    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

    #Precomputations
    rebound.rebound.computeLinearStepList()
    rebound.rebound.computeInverseLinearStepList()
    rebound.rebound.computeDDT(GOST.sub)
    rebound.rebound.computeSBOX_MATCH(GOST.sub)
    
    if(len(sys.argv) > 2):
        if(sys.argv[1] == "save"):
            print "Saving states to " + sys.argv[2]
            startinthemiddle.startFromTheMiddle(sys.argv[2])
        if(sys.argv[1] == "load"):
            print "Loading states from " + sys.argv[2]
            startinthemiddle.findMessagePairFromLoadedState(sys.argv[2])
        if(sys.argv[1] == "print"):
            startinthemiddle.printSavedStates(sys.argv[2])
    else:
        startinthemiddle.startFromTheMiddle("tmp.out")
        
    pass