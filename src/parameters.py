"""
Default parameters
"""

nlevels = 5   # Individuals | Counties | States | Divisions | Country. USA has 5 levels
ml = 1   # Level of entities to be described, i.e. T^{(ml,nl)} or R^{(ml,nl)}
nl = 0   # Level of subentities to be counted (ml=1, nl=0 means counting people infected in each county)

seed = 100     # Initial number of infections
noise = 0      # Average number of new infections imported at each time step.
IC = 'random'    # 'random': seed is randomly spread among all entities
random_seed = 232323   # Seed for random number generator

Nt = 100      # Time steps to evolve

beta = 2   # Constant by which probability matrix T is multiplied. For USA, beta=R0 (if pR=1)
pR = 1     # Recovery probability
distrib = 'binomial'     # Infection distribution: 'binomial', 'poisson', or 'neg_binomial'
sspread = 0.01    # Degree of superspreading [ variance = mean*(1+sspread) ]

model = 'SISreI'        # 'SIS'
                        # 'SISreI': SIS with the possibility of reinfection in the same time-step as recovery
                        # 'SISlinear': Linear low-density limit of SISreI
rescale_usa_fac = 1   # Rescale number of people by this factor, in case we want to do a (0,0) simulation with fewer people

PolicyChoice = '1scale'         # '1scale': RED/GREEN 1-scale strategy
                                # '2scales': R/G 2-scale strategy
                                # '2scales-nointra': Same but without intra-mpol2-level travel restrictions
                                # '3scales': R/G 3-scale strategy
mpol1 = 1     # Level 1 at which policy applied (must be mpol1>= ml)
mpol2 = 2     # For 2scales/3scales policies: Level 2 of policy
mpol3 = 3     # For 3scales policy: Level 3 of policy (must have mpol1<mpol2<mpol3)
RtoGthreshold = 0    # Number of infections below which (<=) policy can turn from RED to GREEN
REDdelay = 0         # Number of time steps to keep RED phase ON when infections are below RtoGthreshold
GtoRthreshold = 1     # Number of infections above which (>=) we switch from GREEN to RED
GREENdelay = 0        # Number of time steps to keep GREEN phase ON when infections are above GtoRthreshold
reductionLocal = 10    # Reduction of transmission probability within region
reductionTravel = 10   # Reduction of transmission probability between regions

adjust00size = True     # Leave as True for USA example
                        # If True, assumes that Tmat is Ntots[1]xNtots[1] when ml=nl=0 
                        # (implies that homogeneity is assumed within level-1 entities when ml=nl=0)


# Parameters dictionary
params_default = { 'nlevels': nlevels, 'ml': ml, 'nl': nl, 'seed': seed, 'noise': noise, 'IC': IC,
         'random_seed': random_seed, 'Nt': Nt,
         'beta': beta, 'pR': pR, 'distrib': distrib, 'sspread': sspread, 'model': model,
         'rescale_usa_fac': rescale_usa_fac, 'PolicyChoice': PolicyChoice, 'mpol1': mpol1, 'mpol2': mpol2,
         'mpol3': mpol3, 'RtoGthreshold': RtoGthreshold, 'REDdelay': REDdelay, 'GtoRthreshold': GtoRthreshold,
         'GREENdelay': GREENdelay, 'reductionLocal': reductionLocal, 'reductionTravel': reductionTravel,
         'adjust00size': adjust00size }
