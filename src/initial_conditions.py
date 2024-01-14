"""
Initial conditions class contains a set of functions to generate the initial state of infections.
"""

import numpy as np

class Initial_Condition:
    """
    This class is for generating vectors of infections to be used as initial conditions of the dynamics.
    """

    REQUIRED_PARAMS = { 'ml', 'nl', 'seed', 'IC' }
    
    def __init__(self, gen, N, Ntots, ichild, **params):
        """Validate and save parameters."""
        
        self._validate_params(params)
        
        self.gen = gen
        self.N = N
        self.Ntots = Ntots
        self.ichild = ichild
        
        self.IC = params['IC']
        self.seed = params['seed']
        self.ml = params['ml']
        self.nl = params['nl']
        

    def _validate_params(self, params):
        """Check that all required paramaters were passed."""
        missing_params = self.REQUIRED_PARAMS - params.keys()
        if missing_params:
            raise ValueError(f"Missing required parameters: {missing_params}")
        
        
    def generate_IC (self):
        """Call initial state generator function depending on value of IC."""
    
        if self.IC=='random':
            
            return self.generate_IC_random()
            
    
    def generate_IC_random (self):
        """Generate random vector of infections."""

        if self.ml == self.nl:
            
            initial_infected = self.gen.choice(self.Ntots[self.ml],self.seed,replace=False)

            Ivec = np.zeros(self.Ntots[self.ml])   # np.repeat(0,self.Ntots[self.ml])
            Ivec[initial_infected] = 1

        else:
            
            Ivec = self.gen.multivariate_hypergeometric(self.N[(self.ml,self.nl)], self.seed)

        return Ivec
    
    
    def generate_IC_onecounty (self,ind):
        """Generate vector with infections in the county with index 'ind'."""
        
        if (self.ml,self.nl) == (1,0):
            
            Ivec = np.zeros(self.Ntots[self.ml])
            Ivec[ind] = self.seed

            return Ivec
        
        if (self.ml,self.nl) == (1,1):
            
            Ivec = np.zeros(self.Ntots[self.ml])
            Ivec[ind] = 1
            
            return Ivec
        
        raise NotImplementedError(f"generate_IC_onecounty only implemented for (ml,nl)=(1,0) and (1,1).")
            
            
    def generate_IC_onestate (self,ind):
        """Generate vector with infections in the state with index 'ind'."""
        
        if (self.ml,self.nl) == (1,0):
            
            state_counties = self.ichild[(2,1)][ind]
            
            Ivec = np.zeros(self.Ntots[self.ml])
            Ivec[state_counties] = self.gen.multivariate_hypergeometric(self.N[(1,0)][state_counties], self.seed)

            return Ivec
        
        if (self.ml,self.nl) == (2,0) or (self.ml,self.nl) == (2,1):
            
            Ivec = np.zeros(self.Ntots[self.ml])
            Ivec[ind] = self.seed

            return Ivec
        
        if (self.ml,self.nl) == (2,2):
            
            Ivec = np.zeros(self.Ntots[self.ml])
            Ivec[ind] = 1
            
            return Ivec
        
        raise NotImplementedError(f"generate_IC_onestate only implemented for (ml,nl)=(1,0), (2,0), (2,1), (2,2).")
            
            
    def generate_IC_allcounties (self):
        """Generate vector with same number of infections in all counties."""
        
        if (self.ml,self.nl) == (1,0):
            
            Ivec = np.full(self.Ntots[self.ml], self.seed)

            return Ivec
        
        raise NotImplementedError(f"generate_IC_allcounties only implemented for (ml,nl)=(1,0).")
    