"""
This module contains functions to set up the index arrays that define the multi-scale hierarchy.
"""

import numpy as np


def set_N_arrays (counties_pop, counties_in_state, states_in_division, rescale=1):
    """Set the number arrays N - they specify the number of sub-entities in each entity."""

    N = {}    # (m,n) entry contains an array where array[ii] gives the number of n-entities contained in
                # the ii-th m-entity

    N[(1,0)] = np.array(counties_pop/rescale, dtype=int)
    N[(2,0)] = np.array( [ np.sum(N[(1,0)][counties_in_state[ii]]) for ii in range(len(counties_in_state)) ] , dtype=int)
    N[(3,0)] = np.array( [ np.sum(N[(2,0)][states_in_division[ii]]) for ii in range(len(states_in_division)) ] , dtype=int)

    N[(4,0)] = np.array( [ np.sum(N[(3,0)]) ], dtype=int)
    N[(4,1)] = np.array( [ len(N[(1,0)]) ], dtype=int)
    N[(4,2)] = np.array( [ len(N[(2,0)]) ], dtype=int)
    N[(4,3)] = np.array( [ len(N[(3,0)]) ], dtype=int)

    Ntots = { 0: N[(4,0)][0], 1: N[(4,1)][0], 2: N[(4,2)][0], 3: N[(4,3)][0], 4: 1 }

    N[(2,1)] = np.array([ len(counties_in_state[ii]) for ii in range(Ntots[2]) ], dtype=int)
    N[(3,1)] = np.array([ np.sum( N[(2,1)][ states_in_division[ii] ] ) for ii in range(Ntots[3]) ], dtype=int)

    N[(3,2)] = np.array([ len(states_in_division[ii]) for ii in range(len(N[(3,0)])) ], dtype=int)

    N[(0,0)] = np.repeat(1, Ntots[0])
    N[(1,1)] = np.repeat(1, Ntots[1])
    N[(2,2)] = np.repeat(1, Ntots[2])
    N[(3,3)] = np.repeat(1, Ntots[3])
    N[(4,4)] = np.repeat(1, Ntots[4])

    return N, Ntots


def set_ichild_arrays (N, Ntots, counties_in_state, states_in_division, level00=False):
    """Set the child arrays - they specify for each entity the indices of its constituent sub-entities."""

    ichild = {}   # (m,n) contains a list where list[ii] contains an array with the indices of the n-entities 
                    # belonging to the ii-th m-entity

    ichild[(1,1)] = list(np.arange(Ntots[1]).reshape((Ntots[1],1)))
    ichild[(2,2)] = list(np.arange(Ntots[2]).reshape((Ntots[2],1)))
    ichild[(3,3)] = list(np.arange(Ntots[3]).reshape((Ntots[3],1)))
    ichild[(4,4)] = list(np.arange(Ntots[4]).reshape((Ntots[4],1)))

    ichild[(2,1)] = [ counties_in_state[ii] for ii in range(Ntots[2]) ]
    ichild[(3,2)] = [ states_in_division[ii] for ii in range(Ntots[3]) ]

    ichild[(3,1)] = [ np.concatenate([ ichild[(2,1)][jj] for jj in ichild[(3,2)][ii] ])  for ii in range(Ntots[3]) ]

    ichild[(4,1)] = [ np.concatenate( ichild[(2,1)] ) ]
    ichild[(4,2)] = [ np.concatenate( ichild[(3,2)] ) ]
    ichild[(4,3)] = [ np.arange(Ntots[3]) ]

    # A lot of memory - not needed unless working at m=n=0 level
    if level00:
        ichild[(0,0)] = list(np.arange(Ntots[0]).reshape((Ntots[0],1)))     

        lseps = np.array([0] + [ np.sum(N[(1,0)]*(np.arange(Ntots[1])<=ii)) for ii in range(Ntots[1]) ])
                    # List marking the indices of individuals separating different clusters
        temp = np.arange(Ntots[0])
        ichild[(1,0)] = [ temp[lseps[cc]:lseps[cc+1]] for cc in range(Ntots[1]) ]    ## A LOT OF MEMORY
        temp = None

        ichild[(2,0)] = [ np.concatenate([ ichild[(1,0)][jj] for jj in ichild[(2,1)][ii] ])  for ii in range(Ntots[2]) ]
        ichild[(3,0)] = [ np.concatenate([ ichild[(1,0)][jj] for jj in ichild[(3,1)][ii] ])  for ii in range(Ntots[3]) ]
        ichild[(4,0)] = [ np.concatenate( ichild[(1,0)] ) ]

    return ichild


def set_iparent_arrays (N, Ntots, ichild, level00=False):
    """Set the parent arrays - they specify for each sub-entity the index of the entity to which they belong."""

    iparent = {}    # Inverse of ichild. (m,n) contains an array where array[ii] contains the index of the 
                    # m-entity to which the ii-th n-entity belongs

    iparent[(1,1)] = np.arange(Ntots[1])
    iparent[(2,2)] = np.arange(Ntots[2])
    iparent[(3,3)] = np.arange(Ntots[3])
    iparent[(4,4)] = np.arange(Ntots[4])

    iparent[(2,1)] = np.array([ np.arange(Ntots[2])[[ ii in ichild[(2,1)][jj] for jj in range(Ntots[2])]][0]\
                                for ii in range(Ntots[1]) ])
    iparent[(3,2)] = np.array([ np.arange(Ntots[3])[[ ii in ichild[(3,2)][jj] for jj in range(Ntots[3])]][0]\
                                for ii in range(Ntots[2]) ])
    iparent[(4,3)] = np.zeros(Ntots[3])

    iparent[(3,1)] = np.array([ iparent[(3,2)][iparent[(2,1)][ii]] for ii in range(Ntots[1]) ])
    iparent[(4,2)] = np.zeros(Ntots[2])

    iparent[(4,1)] = np.zeros(Ntots[1])

    # A lot of memory - not needed unless working at m=n=0 level
    if level00:
        iparent[(0,0)] = np.arange(Ntots[0])
        iparent[(1,0)] = np.repeat(np.arange(Ntots[1]), N[(1,0)] )
        iparent[(2,0)] = np.repeat(iparent[(2,1)][range(Ntots[1])], N[(1,0)] )
        iparent[(3,0)] = np.repeat(iparent[(3,1)][range(Ntots[1])], N[(1,0)] )
        iparent[(4,0)] = np.zeros(Ntots[0])


    return iparent

    