"""
This module contains the Dynamics class, which sets up the system and runs the time evolution.
"""

import numpy as np
from observables import Observables
from policy import Policy
from initial_conditions import Initial_Condition

class Dynamics:
    """
    This class is for bringing all aspects of the disease spread system together.
    It stores the T-matrix, initial condition, observables, and policy information.
    It contains modules to initialize infections and run the dynamics.
    """
    
    REQUIRED_PARAMS = { 'nlevels', 'ml', 'nl', 'seed', 'noise', 'IC', 'random_seed', 'Nt',
                        'beta', 'pR', 'distrib', 'sspread', 'model', 'rescale_usa_fac',
                        'PolicyChoice', 'mpol1', 'mpol2', 'mpol3', 'RtoGthreshold', 'REDdelay',
                        'GtoRthreshold', 'GREENdelay', 'reductionLocal', 'reductionTravel' ,
                        'adjust00size'}


    def __init__(self, N, Ntots, ichild, iparent, Tmat, **params):
        """Validate, check, and store parameters defining the system."""
        
        self._validate_params(params)
        self._check_params(params)

        self.params = params
        
        self.N = N
        self.Ntots = Ntots
        self.ichild = ichild
        self.iparent = iparent
        self.beta = params['beta']
        self.mpol1 = params['mpol1']
        self.mpol2 = params['mpol2']
        self.mpol3 = params['mpol3']
        self.ml = params['ml']
        self.nl = params['nl']
        self.pR = params['pR']
        self.model = params['model']
        self.distrib = params['distrib']
        self.sspread = params['sspread']
        self.Nt = params['Nt']
        self.noise = params['noise']
        self.adjust00size = params['adjust00size']
        
        
        # Random number generator
        self.gen = np.random.default_rng(params['random_seed'])
        
        self.Tmat = Tmat * self.beta
        self.Tmat0 = self.Tmat.copy()  # Keep copy of original Tmat
        
        # We assume that for (ml,nl) dynamics we are given T^(nl,nl). When nl=0 Tmat is saved as Ntots[1]xNtots[1] tough
        self.Tnl = self.nl
        if self.nl==0:
            self.Tnl=1
    

    def _validate_params(self, params):
        """Check that all required paramaters were passed."""
        missing_params = self.REQUIRED_PARAMS - params.keys()
        if missing_params:
            raise ValueError(f"Missing required parameters: {missing_params}")
        
    def _check_params(self, params):
        """Check that some important parameters are within bounds."""
        
        if params['ml']<params['nl']:
            raise ValueError(f"Levels requirement: ml >= nl")
        
        if params['mpol1']<params['ml']:
            raise ValueError(f"Must have: mpol1 >= ml")
        
        if (params['PolicyChoice'] in ['2scales','2scales-nointra']) and params['mpol1']>=params['mpol2']:
            raise ValueError(f"2scales policy requires: mpol1 < mpol2")
        
        if (params['PolicyChoice'] in ['3scales']) and (params['mpol1']>=params['mpol2'] or params['mpol2']>=params['mpol3']):
            raise ValueError(f"3scales policy requires: mpol2 < mpol3")
        
        
    
    def initialize(self):
        """Initialize main objects required for dynamics."""
        
        self.Tmat = self.Tmat0.copy() # Restart Tmat in case of repeated call to function
        
        self.ic = Initial_Condition(self.gen, self.N, self.Ntots, self.ichild, **self.params)
        self.pol = Policy(self.N, self.Ntots, self.ichild, self.iparent, self.Tmat0, **self.params)
        self.obs = Observables(self.N, self.Ntots, self.ichild, **self.params)
        
        
    def compute_dynamics(self):
        """Compute time evolution of number of infections Ivec."""
        
        self.initialize()
        
        # Initial condition
        self.Ivec = self.ic.generate_IC()
        self.t = 0
        
        # Update policy
        self.Tmat = self.pol.apply_policy(self.Ivec,self.Tmat)
        
        # Save t=0 observables
        self.obs.save_observables(self.t, self.Ivec, self.pol, self.Tmat)
        
        for tt in range(self.Nt):
            
            # Compute next step
            self.evolve_one_step()
            
            # Increase time
            self.t = self.t+1
            
            # Save observables
            self.obs.save_observables(self.t, self.Ivec, self.pol, self.Tmat)
            
        self.obs.transform_to_numpy()
    
    
    def compute_NIcounty_during_infection(self):
        """
        Calculate the number of infections in a county before recovery.
        
        Assumes that Tmat is the diagonal part of T^(1,0) only.
        It calculates the total number of infections that a county will accumulate from infection
        till recovery by assuming that it is dominated by local spread. That is, it ignores
        second-order infections where the county infects another county and is infected back.
        """
        
        self.initialize()
        
        # Initial condition
        self.Ivec = self.ic.generate_IC_allcounties()
        self.t = 0
        
        # Update policy
        self.Tmat = self.pol.apply_policy(self.Ivec,self.Tmat)
        
        # Save t=0 observables
        self.obs.save_observables(self.t, self.Ivec, self.pol, self.Tmat)
        
        # Dynamics until no infections in county chosen
        while (np.sum(self.Ivec)>0 and self.t<self.Nt):
            
            # Compute next step
            self.evolve_one_step()
            
            # Increase time
            self.t = self.t+1
            
            # Save observables
            self.obs.save_observables(self.t, self.Ivec, self.pol, self.Tmat)
            
        self.obs.transform_to_numpy()
        
    
    def compute_R1(self,ind):
        """
        Compute reproduction number R1 for county 'ind'.

        It approximates R1 by computing the average number of people infected outside county 'ind'
        while the county stays infected.
        
        Note: input 'ind' must be a valid county index.
        """
        
        self.initialize()
        
        # Initial condition
        self.Ivec = self.ic.generate_IC_onecounty(ind)
        self.t = 0
        
        # Update policy
        self.Tmat = self.pol.apply_policy(self.Ivec,self.Tmat)
        
        # Save t=0 observables, including R1 term
        self.obs.save_R1_part(ind, self.Ivec, self.Tmat)
        self.obs.save_observables(self.t, self.Ivec, self.pol, self.Tmat)
        
        # Dynamics until no infections in county chosen
        while (self.Ivec[ind]>0 and self.t<self.Nt):
            
            # Compute next step
            self.evolve_one_step()
            
            # Increase time
            self.t = self.t+1
            
            # Save observables
            self.obs.save_R1_part(ind, self.Ivec, self.Tmat)
            self.obs.save_observables(self.t, self.Ivec, self.pol, self.Tmat)
            
        self.obs.save_td_county(self.t)
        self.obs.transform_to_numpy()
        
        
    def compute_R2(self,ind):
        """Compute reproduction number R2 for state 'ind'.

        It approximates R2 by computing the average number of people infected outside state 'ind'
        while the state stays infected.
        
        Note: input 'ind' must be a valid state index.
        """
        
        self.initialize()
        
        state_counties = self.ichild[(2,1)][ind]
        
        # Set up initial condition
        self.Ivec = self.ic.generate_IC_onestate(ind)
        self.t = 0
        
        # Update policy
        self.Tmat = self.pol.apply_policy(self.Ivec,self.Tmat)
        
        # Save t=0 observables, including R2 term
        self.obs.save_R2_part(ind, self.Ivec, self.Tmat)
        self.obs.save_observables(self.t, self.Ivec, self.pol, self.Tmat)
        
        # Dynamics until no infections in state chosen
        while (np.sum(self.Ivec[state_counties])>0 and self.t<self.Nt):
            
            # Compute next step
            self.evolve_one_step()
            
            # Increase time
            self.t = self.t+1
            
            # Save observables
            self.obs.save_R2_part(ind, self.Ivec, self.Tmat)
            self.obs.save_observables(self.t, self.Ivec, self.pol, self.Tmat)
            
        self.obs.save_td_state(self.t)
        self.obs.transform_to_numpy()
            
    
        
    def evolve_one_step(self):
        """Evolve system one time step."""
        
        # Compute probability of infection and recovery
        probI, probR = self.compute_probability_infection()
        
        # Update state
        self.update_state(probI,probR)
        
        # Add noise
        if self.noise>0: self.add_noise()
            
        # Update policy
        self.Tmat = self.pol.apply_policy(self.Ivec,self.Tmat)
        
    
    
    def add_noise(self):
        """Add random number of infections."""

        new_infections = self.gen.poisson(self.noise)
        
        if self.ml == self.nl:
            
            infected = self.gen.choice( np.arange(len(self.Ivec))[self.Ivec==0], new_infections, replace=False )
            
            self.Ivec[infected] = 1

        else:
            
            self.Ivec += self.gen.multivariate_hypergeometric(self.N[(self.ml,self.nl)], new_infections)
            
        # Cutoff to avoid more infections than number of people
        self.Ivec = np.minimum(self.Ivec, self.N[(self.ml,self.nl)])

        
        
    def compute_probability_infection(self):
        """
        Compute probability of infection for each ml-level entity.
        
        Assumes that Tmat_ij gives the probability that an infected nl-entity in ml-entity j will infect
        an nl-entity in ml-entity i. This is different from the next-generation matrix R.
        Multiplying the code's Tmat_ij by a factor of N^(m,n)_i gives the R matrix, i.e. the expected number
        of n-entities infected in i by an infection in j.
        """
        
        if self.model in ['SIS','SISreI']:
                
            if self.ml==0 and self.nl==0:   ## Assumes Ivec has length Ntots[0], but Tmat is of size Ntots[1]**2  but represents T^{(0,0)}

                NIc = np.array([ np.sum( self.Ivec[self.ichild[(1,0)][cc]] ) for cc in range(self.Ntots[1]) ]) # number of infected in each cluster
                boolIc = NIc>0
                probNotI = np.prod( (1-self.Tmat[:,boolIc])**(NIc[boolIc].reshape((1,len(NIc[boolIc])))) , axis=1)
                probNotI = np.repeat(probNotI, self.N[(1,0)]) # broadcast to N-length vector


            if self.ml==self.nl and self.ml>0:   ## Assumes Ivec has length Ntots[ml], and Tmat is of size Ntots[ml]**2

                boolIc = self.Ivec>0
                probNotI = np.prod( 1-self.Tmat[:,boolIc]*self.Ivec[None,boolIc] , axis=1)

            if self.ml>self.nl:   ## Assumes Ivec has length Ntots[ml], and Tmat is of size Ntots[ml]**2

                boolIc = self.Ivec>0
                probNotI = np.prod( (1-self.Tmat[:,boolIc])**(self.Ivec[boolIc].reshape((1,len(self.Ivec[boolIc])))) , axis=1)
            
            probI = 1-probNotI
            
        if self.model == 'SISlinear':
                
            if self.ml==0 and self.nl==0:   ## Assumes Ivec has length Ntots[0], but Tmat is of size Ntots[1]**2

                NIc = np.array([ np.sum( self.Ivec[self.ichild[(1,0)][cc]] ) for cc in range(self.Ntots[1]) ]) # number of infected in each cluster
                probI = self.Tmat @ NIc
                probI = np.repeat(probI, self.N[(1,0)]) # broadcast to N-length vector

            if self.ml>self.nl:   ## Assumes Ivec has length Ntots[ml], and Tmat is of size Ntots[ml]**2

                probI = self.Tmat @ self.Ivec
        
        return probI, self.pR
        
        
    def update_state(self,probI,probR):
        """Update state of the system Ivec at next time step."""
        
        # Compute final probability of infection
        if self.model=='SIS':
            
            Ivec_norm = self.Ivec/self.N[(self.ml,self.nl)]
            probI = probI*(1-Ivec_norm) + (1-probR)*Ivec_norm
        
        if self.model in ['SISreI','SISlinear']:
            
            probI = probI
        
        
        # Update state
        if self.ml==self.nl:
            
            self.Ivec = self.gen.binomial( n=self.N[(self.ml,self.nl)], p=probI )
            
        if self.ml>self.nl:
            
            if self.distrib == 'binomial':
                self.Ivec = self.gen.binomial( n=self.N[(self.ml,self.nl)], p=probI )

            if self.distrib == 'poisson':
                Ivec = self.gen.poisson( lam=probI * self.N[(self.ml,self.nl)] )

            if self.distrib == 'neg_binomial':
                if np.sum(self.Ivec)>0:
                    boolExposed = probI>0
                    mu = self.N[(self.ml,self.nl)][boolExposed] * probI[boolExposed]  # mean value = NI*prob
                    n = mu/self.sspread
                    p = 1/(1+self.sspread)

                    self.Ivec[boolExposed] = self.gen.negative_binomial( n=n , p=p  )
                    self.Ivec[np.logical_not(boolExposed)] = 0

        # Cutoff to avoid more infections than number of people (can happen for poisson and neg_binomial)
        self.Ivec = np.minimum(self.Ivec, self.N[(self.ml,self.nl)])
        
        