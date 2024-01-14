"""
This module contains the Observables class to compute and store different observable quantities.
"""

import numpy as np

class Observables:
    """
    This class is for computing and storing different observables related to the infection vector,
    probability matrix and policies.
    """
    
    REQUIRED_PARAMS = { 'nlevels', 'ml', 'nl', 'mpol1', 'mpol2', 'mpol3' }

    def __init__(self, N, Ntots, ichild, **params):
        """Save parameters and set up arrays for storage of observables."""
        
        self._validate_params(params)

        self.N = N
        self.Ntots = Ntots
        self.ichild = ichild

        self.nlevels = params['nlevels']
        self.ml = params['ml']
        self.nl = params['nl']
        self.mpol1 = params['mpol1']
        self.mpol2 = params['mpol2']
        self.mpol3 = params['mpol3']
        
        self.times = []
        self.NItot = []     # Total infections (same as self.NI[self.nlevels-1])
        
        # Number of infections at different levels (except 0 level). Lists of levels m<ml stay empty.
        self.NI = {aa: [] for aa in range(1,self.nlevels)}  
        
        self.policy = {aa: [] for aa in range(1,4)}   # Assumes maximal 3-scale policy
        self.numRED = {aa: [] for aa in range(1,4)}   # Assumes maximal 3-scale policy
        self.peopleT = []    # Cumulative number of people under local restrictions
        
        self.R0 = []
        self.R1 = 0
        self.R2 = 0
        self.td_county = 0
        self.td_state = 0


    def _validate_params(self, params):
        """Check that all required paramaters were passed."""
        missing_params = self.REQUIRED_PARAMS - params.keys()
        if missing_params:
            raise ValueError(f"Missing required parameters: {missing_params}")
        
        
    def save_observables(self, time, Ivec, pol, Tmat):
        """Compute and save observables."""
        
        self.times.append(time)
        
        NItot = np.sum(Ivec)
        self.NItot.append(NItot)
        
        # Infections per region (nn,self.nl)-level
        for nn in range(np.max([1,self.ml]),self.nlevels):
            
            NI_nn = np.array([ np.sum( Ivec[self.ichild[(nn,self.ml)][cc]] ) for cc in range(self.Ntots[nn]) ])
            self.NI[nn].append( NI_nn )
        
        # Note: if policy has less than 3 scales then the remaining arrays will be empty
        self.numRED[1].append( np.sum( pol.current_policy1=='r' ) )
        self.numRED[2].append( np.sum( pol.current_policy2=='r' ) )
        self.numRED[3].append( np.sum( pol.current_policy3=='r' ) )
        
        self.policy[1].append( np.array( pol.current_policy1=='r', dtype=int) )
        self.policy[2].append( np.array( pol.current_policy2=='r', dtype=int) )
        self.policy[3].append( np.array( pol.current_policy3=='r', dtype=int) )
        
        previous = 0
        if len(self.peopleT)>0:
            previous = self.peopleT[-1]
            #self.peopleT.append( self.peopleT[-1] + np.sum(self.N[(self.mpol1,0)]*(pol.current_policy1=='r')) )
            #self.peopleT.append( np.sum(self.N[(self.mpol1,0)]*(pol.current_policy1=='r')) )
        self.peopleT.append( previous + np.sum(self.N[(self.mpol1,0)]*(pol.current_policy1=='r')) )
        
        # R0
        if self.ml==1 and self.nl==0:
            self.R0.append( np.sum(Tmat * self.N[(1,0)][:,None] * self.N[(1,0)][None,:] )/self.Ntots[0] )
        
        
    def save_R1_part(self, ind, Ivec, Tmat):
        """
        Compute and save contribution to R1 from county 'ind'.

        Specifically, the mean number of people is computed that the county with
        index 'ind' (must be valid index) will infect in next time step.
        This is an approximation to R1.
        """
        
        if self.ml==1 and self.nl==0:

            Tmat_ind = Tmat[:,ind].copy()
            Tmat_ind[ind] = 0
            pleaving = np.sum( Tmat_ind * self.N[(1,0)] )
            
            self.R1 += Ivec[ind] * pleaving

        else:
            raise NotImplementedError(f"compute_R1_part function doesnt work for ml!=1 or nl!=0.")
            
    def save_R2_part(self, ind, Ivec, Tmat):
        """
        Compute and save contribution to R2 from state 'ind'.

        Specifically, the mean number of people is computed that the state with
        index 'ind' (must be valid index) will infect in next time step.
        This is an approximation to R2.
        """

        if self.ml==1 and self.nl==0:
            state_counties = self.ichild[(2,1)][ind]
            
            Tmat_ind = Tmat[:,state_counties].copy()
            Tmat_ind[state_counties,:] = 0
            pleaving = np.sum( Tmat_ind * self.N[(1,0)][:,None], axis=0 )
            
            self.R2 += np.sum(Ivec[state_counties] * pleaving)

        else:
            raise NotImplementedError(f"compute_R2_part function doesnt work for ml!=1 or nl!=0.")
            
            
        
    def save_td_county(self, time):
        """Save current time as td_county."""
        
        self.td_county = time
        
    def save_td_state(self, time):
        """Save current time as td_state."""
        
        self.td_state = time
        
        
    def transform_to_numpy(self):
        """Transform observable lists into numpy arrays."""
        
        self.times = np.array(self.times)
        self.NItot = np.asarray(self.NItot)
        
        for key in self.NI:
            self.NI[key] = np.asarray(self.NI[key])
        
        for key in self.numRED:
            self.numRED[key] = np.asarray(self.numRED[key])
            
        for key in self.policy:
            self.policy[key] = np.asarray(self.policy[key])
        
        self.peopleT = np.asarray(self.peopleT)
        self.R0 = np.array(self.R0)
        
        
        