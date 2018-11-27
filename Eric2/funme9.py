from __future__ import print_function
import math
import mlpotential_e 

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

#### read in parameters ####
nsweeps = 9
nlayers = 5
parameterfilename = "../Eric/parameters0b"
#feidatafile       = "../Eric/h2-aimd.fei"
feidatafile       = "spring.fei"
energygradient0 = lambda x: (0.0, [0.0]*len(x))

mlpotential = mlpotential_e.mlpotential_e(nsweeps,nlayers,parameterfilename,feidatafile,energygradient0)

