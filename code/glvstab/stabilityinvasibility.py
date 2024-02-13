#!/usr/bin/env python3

__appname__ = '[stabilityinvasibility.py]'
__author__ = 'Pablo Lechon (plechon@uchicago.edu)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import pandas as pd
import numpy as np
import itertools

## CONSTANTS ##


## FUNCTIONS ##

def solveglv(A, r):
    if len(r) > 1:
        return np.dot(-np.linalg.inv(A),r)
    else:
        return -1/A*r

def jacobianglv(A, xstar):
    return xstar*A

def invasibility(xstar, Aj, rj):
    invcrit = rj + np.dot(Aj, xstar)
    return invcrit > 0

def checkallinvasions(A, r, xstar, presence):
    indabsent = np.where(presence == 0)[0]
    indpresent = np.where(presence == 1)[0]
    invasible = np.zeros(len(indabsent))
    storeind = 0
    for inv in indabsent:
        #try to invade with species inv
        Ainv = A[inv, indpresent]
        rinv = r[inv]
        invasible[storeind] = invasibility(xstar, Ainv, rinv)
    return any(invasible)



def presence_absence_matrix(n):
    combinations = itertools.product([0, 1], repeat=n)
    matrix = [list(combination) for combination in combinations if any(combination)]
    return np.array(matrix)

def pruneA(matrix, presence):
    indstodelete = np.where(presence==0)[0]
    return np.delete(np.delete(matrix, indstodelete, 0), indstodelete, 1)

def pruner(vector, presence):
    return np.array([vector[i] for i in range(len(vector)) if presence[i] == 1])

def embedequilibria(xsmall, presence):
    recipientvec = np.zeros(len(presence))
    inds = np.where(presence != 0)
    recipientvec[inds] = xsmall
    return recipientvec

def storeresults(n, sim, d, nfinal, subcom, xstar, lambdamax, invasible):
    #express vectors in long format
    nvec = n*[n]
    simvec = n*[sim]
    dvec = n*[d]
    nfinalvec = n*[nfinal]
    subcomvec = n*[subcom]
    lambdamaxvec = n*[lambdamax]
    invasiblevec = n*[invasible]
    spvec = list(range(n))
    tostore = np.stack([nvec, simvec, dvec, nfinalvec, subcomvec, \
            spvec, list(xstar), lambdamaxvec, invasiblevec])
    return np.transpose(tostore)

def main(argv):
    '''Main function'''
    #set simulation parameters
    nsims = 10
    #set model parameters
    nmax = 6
    dmax = 2
    step = 0.01
    #initialize storage
    storingmat = np.empty((0, 9))
    for sim in range(nsims):
        for n in range(2,nmax):
            print("Diversity ", n)
            d = 0
            C = -np.random.rand(n, n) #positive mean
            r = np.ones(n)
            presencemat = presence_absence_matrix(n)
            while d <= dmax:
                M = C - d * np.eye(n)
                for subcom in range(2**n - 1):
                    #get the presence absence vector
                    presencevec = presencemat[subcom, :]
                    nfinal = sum(presencevec)
                    #get parameters of subcomm
                    Msub = pruneA(M, presencevec)
                    rsub = pruner(r, presencevec)
                    #get equilibria of system
                    xsub = solveglv(Msub, rsub)
                    #if feasible, get eigenvalues
                    if (xsub >= 0).all():
                        jacmat = jacobianglv(Msub, xsub)
                        #find maximum eigenvalue
                        lambdamax = max(np.real(np.linalg.eigvals(jacmat)))
                        #if stable check for uninvasibility
                        if lambdamax < 0:
                            is_invasible = checkallinvasions(M, r, xsub,\
                                                             presencevec)
                            #embed in original population
                            xstar = embedequilibria(xsub, presencevec)
                            #create matrix to store
                            submatrix = storeresults(n, sim, d, nfinal, subcom, xstar,\
                                    lambdamax, is_invasible)
                            #add to bigger matrix
                            storingmat = np.append(storingmat, submatrix, axis = 0)
                            #otherwise, jump to next subcommunity
                    else:
                        continue
                d += step
    #save results
    df = pd.DataFrame(storingmat)
    df.to_csv('stabresults.csv', index=False, header = False)
    return 0

## CODE ##

if (__name__ == '__main__'):
   status = main(sys.argv)
   sys.exit(status)
	    

