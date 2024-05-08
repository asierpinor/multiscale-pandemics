"""
This module contains the Policy class to modify transmission probabilities due to policy restrictions.
"""

import numpy as np
import time

class Policy:
    """
    This class stores the current policy status of the system and contains a number of methods
    to update the policy status based on the current vector of infections, and to modify the
    probability transmission matrix T accordingly. The policies implemented can act on 1, 2, or 3
    different scales and the restrictions are divided into local and travel.
    """

    REQUIRED_PARAMS = { 'ml', 'nl', 'PolicyChoice', 'mpol1', 'mpol2', 'mpol3', 'RtoGthreshold', 'REDdelay',
                        'GtoRthreshold', 'GREENdelay', 'reductionLocal', 'reductionTravel', 'adjust00size' }
    
    def __init__(self, N, Ntots, ichild, iparent, Tmat0, **params):
        """Validate and store parameters and initialize policy arrays."""

        self._validate_params(params)

        self.N = N
        self.Ntots = Ntots
        self.ichild = ichild
        self.iparent = iparent
        self.Tmat0 = Tmat0
        
        self.ml = params['ml']
        self.nl = params['nl']
        self.mpol1 = params['mpol1']
        self.mpol2 = params['mpol2']
        self.mpol3 = params['mpol3']
        self.PolicyChoice = params['PolicyChoice']
        self.GtoRthreshold = params['GtoRthreshold']
        self.RtoGthreshold = params['RtoGthreshold']
        self.GREENdelay = params['GREENdelay']
        self.REDdelay = params['REDdelay']
        self.reductionTravel = params['reductionTravel']
        self.reductionLocal = params['reductionLocal']
        self.adjust00size = params['adjust00size']
        
        # Take into account that Tmat is Ntots[1]xNtots[1] when ml=0 (if adjust00size==True)
        self.Tml = self.ml
        if self.ml==0 and self.adjust00size:
            self.Tml=1
        
        self.set_policy_counters()
        
    
    def _validate_params(self, params):
        """Check that all required paramaters were passed."""
        missing_params = self.REQUIRED_PARAMS - params.keys()
        if missing_params:
            raise ValueError(f"Missing required parameters: {missing_params}")
        

    def set_policy_counters (self):
        """Initialize policy status arrays to all green."""
        
        # Arrays need to exist for Observables class to work well
        self.current_policy1 = np.array(['g'])
        self.current_policy2 = np.array(['g'])
        self.current_policy3 = np.array(['g'])
        
        if self.PolicyChoice in ['1scale','2scales','2scales-nointra','3scales']:
            
            self.nclusters1 = self.Ntots[self.mpol1]
            self.cluster1_children = self.ichild[(self.mpol1,self.ml)]
            self.cluster1_Tml_children = self.ichild[(self.mpol1,self.Tml)]
            
            self.REDcounter1 = np.repeat(0,self.nclusters1)
            self.GREENcounter1 = np.repeat(0,self.nclusters1)
            self.current_policy1 = np.array(['g']*self.nclusters1)
            
            
        if self.PolicyChoice in ['2scales','2scales-nointra','3scales']:
            
            self.nclusters2 = self.Ntots[self.mpol2]
            self.cluster2_children = self.ichild[(self.mpol2,self.mpol1)]
            self.cluster2_Tml_children = self.ichild[(self.mpol2,self.Tml)]  # Tmat saved as N[1]xN[1] matrix
            
            self.REDcounter2 = np.repeat(0,self.nclusters2)
            self.GREENcounter2 = np.repeat(0,self.nclusters2)
            self.current_policy2 = np.array(['g']*self.nclusters2)
            
            if self.PolicyChoice=='2scales-nointra':
                
                self.intra_reduction = False   # Only red level1 entities reduce intra-level2-travel within a red level2 entity
                self.parents2_cluster1 = self.iparent[(self.mpol2,self.mpol1)]
                
            if self.PolicyChoice=='2scales':
                
                self.intra_reduction = True    # All level1 entities within a red level2 entities reduce intra-level2-travel
                
        
        if self.PolicyChoice in ['3scales']:
            
            self.nclusters3 = self.Ntots[self.mpol3]
            self.cluster3_children = self.ichild[(self.mpol3,self.mpol2)]
            self.cluster3_Tml_children = self.ichild[(self.mpol3,self.Tml)]  # Tmat saved as N[1]xN[1] matrix
            
            self.REDcounter3 = np.repeat(0,self.nclusters3)
            self.GREENcounter3 = np.repeat(0,self.nclusters3)
            self.current_policy3 = np.array(['g']*self.nclusters3)
            
        
        self.current_policy = self.current_policy1
            
        
        
    def update_policy (self, cases, current_policy, GREENcounter, REDcounter, GREENdelay, REDdelay,
                       GtoRthreshold, RtoGthreshold):
        """Update status of policy based on current vector of infections."""

        g_clusters_before = (current_policy == 'g')
        r_clusters_before = (current_policy == 'r')

        ### Green clusters
        # Increase counter if above threshold
        above_hGtoR = ((cases>=GtoRthreshold) & g_clusters_before)
        GREENcounter[above_hGtoR] += 1
        #GREENcounter[~above_hGtoR] = 0
        
        # Change to Red if after delay
        change_GtoR = (above_hGtoR & (GREENcounter>GREENdelay))
        current_policy[change_GtoR] = 'r'
        GREENcounter[change_GtoR] = 0
        REDcounter[change_GtoR] = 0
        
        ### Red clusters
        # Increase counter if below threshold
        below_hRtoG = ((cases<=RtoGthreshold) & r_clusters_before)
        REDcounter[below_hRtoG] += 1
        #REDcounter[~below_hRtoG] = 0
        
        # Change to Green if after delay
        change_RtoG = (below_hRtoG & (REDcounter>REDdelay))
        current_policy[change_RtoG] = 'g'
        GREENcounter[change_RtoG] = 0
        REDcounter[change_RtoG] = 0
        
        return change_GtoR, change_RtoG
    
    
    def compute_indices_lower_level (self, selector, nclusters_high, children_high_to_low):
        """Combine the set of indices of all the sub-entities belonging to entities selected by the selector."""
        
        indices_clusters = np.arange(nclusters_high)[selector]
        indices = np.concatenate([children_high_to_low[cc] for cc in indices_clusters] + [np.array([],dtype=int)])
        
        return indices
    
    
        
    def apply_policy (self,Ivec,Tmat):
        """Select and apply policy based on PolicyChoice value."""
            
        if self.PolicyChoice == '1scale': return self.apply_policy_1scale_vectorized(Ivec,Tmat)
        
        if self.PolicyChoice in ['2scales','2scales-nointra']:
            return self.apply_policy_2scales_forloop(Ivec,Tmat)
        
        if self.PolicyChoice == '3scales': return self.apply_policy_3scales_simple(Ivec,Tmat)
        

    
    def apply_policy_1scale_simple (self,Ivec,Tmat_dummy):
        """"
        Update 1-scale policy status and modify transmission probabilities accordingly.

        Simple implementation where previous policy status is not used.
        """
        
        Tmat = self.Tmat0.copy()  # Tmat_dummy is not used, only included to make all functions have the same inputs
    
        NIc = np.array([ np.sum( Ivec[self.cluster1_children[cc]] ) for cc in range(self.nclusters1) ])
        
        change_GtoR, change_RtoG = self.update_policy(NIc, self.current_policy1, self.GREENcounter1,
                                                      self.REDcounter1, self.GREENdelay, self.REDdelay,
                                                      self.GtoRthreshold, self.RtoGthreshold)
        
        # New policy status
        g_clusters_after = (self.current_policy1 == 'g')
        r_clusters_after = (self.current_policy1 == 'r')
        
        
        # Compute indices of R and G Tml-entities
        indices_R = self.compute_indices_lower_level(r_clusters_after, self.nclusters1, self.cluster1_Tml_children)
        indices_G = self.compute_indices_lower_level(g_clusters_after, self.nclusters1, self.cluster1_Tml_children)
        
        # Compute indices of R and G mpol1-scale entities
        indices_clusters_R = np.arange(self.nclusters1)[r_clusters_after]

        # Update Tmat
        # Travel changes if region 1 is G (R), and region 2 is R (G), or both R
        Tmat[indices_R[:,None],indices_G[None,:]] *= 1/self.reductionTravel
        Tmat[indices_G[:,None],indices_R[None,:]] *= 1/self.reductionTravel
        
        Tmat[indices_R[:,None],indices_R[None,:]] *= 1/self.reductionTravel
        
        # Local changes if R. Compensate for diagonal part of travel changes.
        for aa in indices_clusters_R:
            
            indices_R_aa = self.cluster1_Tml_children[aa]
            Tmat[indices_R_aa[:,None],indices_R_aa[None,:]] *= self.reductionTravel/self.reductionLocal
        
        
        return Tmat
    
    
                
    def apply_policy_1scale_forloop (self,Ivec,Tmat):
        """"
        Update 1-scale policy status and modify transmission probabilities accordingly.

        Implementation using a for loop that updates policy and modifies T sequentially for each entity.
        """
        
        NIc = np.array([ np.sum( Ivec[self.cluster1_children[cc]] ) for cc in range(self.nclusters1) ])

        for cc in range(self.nclusters1):

            if self.current_policy1[cc]=='g':

                if NIc[cc]>=self.GtoRthreshold:

                    self.GREENcounter1[cc] = self.GREENcounter1[cc] + 1

                    if self.GREENcounter1[cc]>self.GREENdelay: 

                        self.current_policy1[cc] = 'r'
                        self.GREENcounter1[cc] = 0
                        self.REDcounter1[cc] = 0
                        
                        g_clusters = (self.current_policy1=='g')
                        g_clusters[cc] = False   # Not necessary, just for clarity
                        
                        indices_G = self.compute_indices_lower_level(g_clusters, self.nclusters1, self.cluster1_Tml_children)
                        indices_R_cc = self.cluster1_Tml_children[cc]
                        
                        Tmat[indices_G[:,None], indices_R_cc[None,:]] *= 1/self.reductionTravel
                        Tmat[indices_R_cc[:,None], indices_G[None,:]] *= 1/self.reductionTravel
                        
                        Tmat[indices_R_cc[:,None], indices_R_cc[None,:]] *= 1/self.reductionLocal

                else:
                    self.GREENcounter1[cc] = 0

            if self.current_policy1[cc]=='r':

                if NIc[cc]<=self.RtoGthreshold:

                    self.REDcounter1[cc] = self.REDcounter1[cc] + 1

                    if self.REDcounter1[cc]>self.REDdelay: 

                        self.current_policy1[cc] = 'g'
                        self.GREENcounter1[cc] = 0
                        self.REDcounter1[cc] = 0
                        
                        g_clusters = (self.current_policy1=='g')
                        g_clusters[cc] = False   # Exclude self
                        
                        indices_G = self.compute_indices_lower_level(g_clusters, self.nclusters1, self.cluster1_Tml_children)
                        indices_G_cc = self.cluster1_Tml_children[cc]
                        
                        Tmat[indices_G[:,None], indices_G_cc[None,:]] *= self.reductionTravel
                        Tmat[indices_G_cc[:,None], indices_G[None,:]] *= self.reductionTravel
                        
                        Tmat[indices_G_cc[:,None], indices_G_cc[None,:]] *= self.reductionLocal

                else:
                    self.REDcounter1[cc] = 0
        
        return Tmat
    
    
    def apply_policy_1scale_vectorized (self,Ivec,Tmat):
        """"
        Update 1-scale policy status and modify transmission probabilities accordingly.

        Vectorized Implementation (like simple one) using information of previous policy status to
        only change T entries that change with the new policy status.
        """
    
        g_clusters_before = (self.current_policy1 == 'g')
        r_clusters_before = (self.current_policy1 == 'r')
        
        # Level 1 policy update
        if self.mpol1==self.ml:
            NIc = Ivec
        else:
            # NOTE: Potential performance bottleneck
            NIc = np.array([ np.sum( Ivec[self.cluster1_children[cc]] ) for cc in range(self.nclusters1) ])
        
        change_GtoR, change_RtoG = self.update_policy(NIc, self.current_policy1, self.GREENcounter1,
                                                      self.REDcounter1, self.GREENdelay, self.REDdelay,
                                                      self.GtoRthreshold, self.RtoGthreshold)

        g_clusters_after = (self.current_policy1 == 'g')
        r_clusters_after = (self.current_policy1 == 'r')
        
        # Compute indices at mpol1-scale
        indices_clusters_GtoR = np.arange(self.nclusters1)[change_GtoR]
        indices_clusters_RtoG = np.arange(self.nclusters1)[change_RtoG]

        # Compute indices at Tml-scale
        if self.mpol1==self.Tml:
            indices_GtoR = indices_clusters_GtoR
            indices_RtoG = indices_clusters_RtoG
            indices_GtoG = np.arange(self.nclusters1)[g_clusters_before & g_clusters_after]
        else:
            # NOTE: Potential performance bottleneck
            indices_GtoR = self.compute_indices_lower_level(change_GtoR, self.nclusters1, self.cluster1_Tml_children)
            indices_RtoG = self.compute_indices_lower_level(change_RtoG, self.nclusters1, self.cluster1_Tml_children)
            indices_GtoG = self.compute_indices_lower_level(g_clusters_before & g_clusters_after, self.nclusters1, self.cluster1_Tml_children)

        # Update Tmat
        # Travel changes if region 1 is GtoR (RtoG), and region to is Gbefore (Gafter)
        Tmat[indices_GtoR[:,None],indices_GtoR[None,:]] *= 1/self.reductionTravel
        Tmat[indices_GtoR[:,None],indices_GtoG[None,:]] *= 1/self.reductionTravel
        Tmat[indices_GtoG[:,None],indices_GtoR[None,:]] *= 1/self.reductionTravel
        
        Tmat[indices_RtoG[:,None],indices_RtoG[None,:]] *= self.reductionTravel
        Tmat[indices_RtoG[:,None],indices_GtoG[None,:]] *= self.reductionTravel
        Tmat[indices_GtoG[:,None],indices_RtoG[None,:]] *= self.reductionTravel

        # Local changes if GtoR. Compensate for diagonal part of travel changes.
        for aa in indices_clusters_GtoR:

            indices_GtoR_aa = self.cluster1_Tml_children[aa]
            Tmat[indices_GtoR_aa[:,None],indices_GtoR_aa[None,:]] *= self.reductionTravel/self.reductionLocal

        # Local changes if RtoG. Compensate for diagonal part of travel changes.
        for aa in indices_clusters_RtoG:

            indices_RtoG_aa = self.cluster1_Tml_children[aa]
            Tmat[indices_RtoG_aa[:,None],indices_RtoG_aa[None,:]] *= self.reductionLocal/self.reductionTravel

        return Tmat
    
    
    
    def apply_policy_2scales_simple (self,Ivec,Tmat_dummy):
        """"
        Update 2-scales policy status and modify transmission probabilities accordingly.

        Simple implementation where previous policy status is not used.
        """
        
        Tmat = self.Tmat0.copy()  # Tmat_dummy is not used, only included to make all functions have the same inputs
        
        # Level 1 policy update
        NIc = np.array([ np.sum( Ivec[self.cluster1_children[cc]] ) for cc in range(self.nclusters1) ])
        
        change_GtoR1, change_RtoG1 = self.update_policy(NIc, self.current_policy1, self.GREENcounter1,
                                                        self.REDcounter1, self.GREENdelay, self.REDdelay,
                                                        self.GtoRthreshold, self.RtoGthreshold)
        
        g1_clusters_after = (self.current_policy1 == 'g')
        r1_clusters_after = (self.current_policy1 == 'r')
        
        # Level 2 policy update
        nReds = np.array([ np.sum(self.current_policy1[self.cluster2_children[cc]]=='r') for cc in range(self.nclusters2) ])
        
        change_GtoR2, change_RtoG2 = self.update_policy(nReds, self.current_policy2, self.GREENcounter2,
                                                        self.REDcounter2, self.GREENdelay, self.REDdelay,
                                                        1, 0)
        
        g2_clusters_after = (self.current_policy2 == 'g')
        r2_clusters_after = (self.current_policy2 == 'r')
        
        
        # Level 1: Compute indices of R and G Tml-entities
        indices_R1 = self.compute_indices_lower_level(r1_clusters_after, self.nclusters1, self.cluster1_Tml_children)
        indices_G1 = self.compute_indices_lower_level(g1_clusters_after, self.nclusters1, self.cluster1_Tml_children)
        
        # Level 2: Compute indices of R and G Tml-entities
        indices_R2 = self.compute_indices_lower_level(r2_clusters_after, self.nclusters2, self.cluster2_Tml_children)
        indices_G2 = self.compute_indices_lower_level(g2_clusters_after, self.nclusters2, self.cluster2_Tml_children)
        
        # Cluster indices
        indices_clusters_R1 = np.arange(self.nclusters1)[r1_clusters_after]
        indices_clusters_G1 = np.arange(self.nclusters1)[g1_clusters_after]
        indices_clusters_R2 = np.arange(self.nclusters2)[r2_clusters_after]
        
        
        # Update Tmat with/without intra-level-2 reduction
        #
        # Level 2
        Tmat[indices_R2[:,None],indices_G2[None,:]] *= 1/self.reductionTravel
        Tmat[indices_G2[:,None],indices_R2[None,:]] *= 1/self.reductionTravel

        Tmat[indices_R2[:,None],indices_R2[None,:]] *= 1/self.reductionTravel
        
        
        ### New version
        # Level 1 - For R1 (=R1R2) add local restrictions and compensate for travel factor (assumes that R1G2 is impossible)        
        for aa in indices_clusters_R1:
                
            indices_R1_aa = self.cluster1_Tml_children[aa]
            Tmat[indices_R1_aa[:,None],indices_R1_aa[None,:]] *= self.reductionTravel/self.reductionLocal
            
        
        # Level 1 - For G1R2 remove restrictions
        for aa in indices_clusters_R2:
            
            indices_clusters_1R2_aa = self.cluster2_children[aa]  # level 1 clusters in level 2 cluster cc
            indices_clusters_G1R2_aa = indices_clusters_1R2_aa[ np.isin(indices_clusters_1R2_aa, indices_clusters_G1) ]
            
            if self.intra_reduction == True:
                
                for bb in indices_clusters_G1R2_aa:

                    indices_G1R2_aabb = self.cluster1_Tml_children[bb]
                    Tmat[indices_G1R2_aabb[:,None],indices_G1R2_aabb[None,:]] *= self.reductionTravel
                    
            else:
                
                indices_G1R2_aa = np.concatenate([self.cluster1_Tml_children[bb] for bb in indices_clusters_G1R2_aa] + [np.array([],dtype=int)])
                Tmat[indices_G1R2_aa[:,None],indices_G1R2_aa[None,:]] *= self.reductionTravel
        
        
        return Tmat
    
    
    def apply_policy_2scales_forloop (self,Ivec,Tmat):
        """"
        Update 2-scales policy status and modify transmission probabilities accordingly.

        Implementation using a for loop that updates policy and modifies T sequentially for each entity.
        """
        
        ##
        ## Level 1 restrictions (mpol1): Local and intra-level2 travel
        ##
        NIc = np.array([ np.sum( Ivec[self.cluster1_children[cc]] ) for cc in range(self.nclusters1) ])
        
        for cc in range(self.nclusters1):

            if self.current_policy1[cc]=='g':

                if NIc[cc]>=self.GtoRthreshold:

                    self.GREENcounter1[cc] = self.GREENcounter1[cc] + 1

                    if self.GREENcounter1[cc]>self.GREENdelay: 

                        self.current_policy1[cc] = 'r'
                        self.GREENcounter1[cc] = 0
                        self.REDcounter1[cc] = 0
                        
                        indices_R1_cc = self.cluster1_Tml_children[cc]
                        Tmat[indices_R1_cc[:,None], indices_R1_cc[None,:]] *= 1/self.reductionLocal
                        
                        # Intra-state travel reduction for red level1 entity
                        # Apply here only if level 2 policy does not restrict all intra-level2 travel
                        # Otherwise apply at level 2 directly
                        if self.intra_reduction == False:
                        
                            # Get indices of level1 regions that are green and within same level2
                            samelevel2 = (self.parents2_cluster1 == self.parents2_cluster1[cc])
                            g1_clusters = (self.current_policy1=='g')
                            g1_clusters[cc] = False   # Not necessary, just for clarity
                            
                            indices_G1_same2 = self.compute_indices_lower_level(g1_clusters & samelevel2, self.nclusters1, self.cluster1_Tml_children)

                            # Apply intra level 1 reductions
                            Tmat[indices_G1_same2[:,None], indices_R1_cc[None,:]] *= 1/self.reductionTravel
                            Tmat[indices_R1_cc[:,None], indices_G1_same2[None,:]] *= 1/self.reductionTravel

                else:
                    self.GREENcounter1[cc] = 0

            if self.current_policy1[cc]=='r':

                if NIc[cc]<=self.RtoGthreshold:

                    self.REDcounter1[cc] = self.REDcounter1[cc] + 1

                    if self.REDcounter1[cc]>self.REDdelay: 

                        self.current_policy1[cc] = 'g'
                        self.GREENcounter1[cc] = 0
                        self.REDcounter1[cc] = 0
                        
                        indices_G1_cc = self.cluster1_Tml_children[cc]
                        Tmat[indices_G1_cc[:,None], indices_G1_cc[None,:]] *= self.reductionLocal
                        
                        # Intra-state travel reduction for red level1 entity
                        # Apply here only if level 2 policy does not restrict all intra-level2 travel
                        # Otherwise apply at level 2 directly
                        if self.intra_reduction == False:
                        
                            # Get indices of level1 regions that are green and within same level2
                            samelevel2 = (self.parents2_cluster1 == self.parents2_cluster1[cc])
                            g1_clusters = (self.current_policy1=='g')
                            g1_clusters[cc] = False   # Here it's necessary
                            
                            indices_G1_same2 = self.compute_indices_lower_level(g1_clusters & samelevel2, self.nclusters1, self.cluster1_Tml_children)

                            # Apply intra level 1 reductions
                            Tmat[indices_G1_same2[:,None], indices_G1_cc[None,:]] *= self.reductionTravel
                            Tmat[indices_G1_cc[:,None], indices_G1_same2[None,:]] *= self.reductionTravel

                else:
                    self.REDcounter1[cc] = 0
        
        ##
        ## Level 2 restrictions (mpol2): Inter and intra-level2 travel
        ##
        
        nReds = [ np.sum(self.current_policy1[self.cluster2_children[cc]]=='r') for cc in range(self.nclusters2) ]
        
        for cc in range(self.nclusters2):

            if self.current_policy2[cc]=='g':

                if nReds[cc] >= 1:

                    self.GREENcounter2[cc] = self.GREENcounter2[cc] + 1

                    if self.GREENcounter2[cc]>self.GREENdelay: 

                        self.current_policy2[cc] = 'r'
                        self.GREENcounter2[cc] = 0
                        self.REDcounter2[cc] = 0
                        
                        g2_clusters = (self.current_policy2=='g')
                        g2_clusters[cc] = False   # Not necessary, just for clarity
                        
                        indices_R2_cc = self.cluster2_Tml_children[cc]
                        indices_G2 = self.compute_indices_lower_level(g2_clusters, self.nclusters2, self.cluster2_Tml_children)
                        
                        # Inter-state
                        Tmat[indices_G2[:,None], indices_R2_cc[None,:]] *= 1/self.reductionTravel
                        Tmat[indices_R2_cc[:,None], indices_G2[None,:]] *= 1/self.reductionTravel
                        
                        # Intra-state
                        if self.intra_reduction == True:
                            
                            Tmat[indices_R2_cc[:,None], indices_R2_cc[None,:]] *= 1/self.reductionTravel

                            # Compensate for local part
                            indices_clusters_1R2_cc = self.cluster2_children[cc]

                            for dd in indices_clusters_1R2_cc:

                                indices_1R2_ccdd = self.cluster1_Tml_children[dd]
                                Tmat[indices_1R2_ccdd[:,None], indices_1R2_ccdd[None,:]] *= self.reductionTravel
                        

                else:
                    self.GREENcounter2[cc] = 0

            if self.current_policy2[cc]=='r':

                if nReds[cc] == 0:

                    self.REDcounter2[cc] = self.REDcounter2[cc] + 1

                    if self.REDcounter2[cc]>self.REDdelay: 

                        self.current_policy2[cc] = 'g'
                        self.GREENcounter2[cc] = 0
                        self.REDcounter2[cc] = 0
                        
                        g2_clusters = (self.current_policy2=='g')
                        g2_clusters[cc] = False   # Here it's necessary
                        
                        indices_G2_cc = self.cluster2_Tml_children[cc]
                        indices_G2 = self.compute_indices_lower_level(g2_clusters, self.nclusters2, self.cluster2_Tml_children)
                        
                        # Inter-state
                        Tmat[indices_G2[:,None], indices_G2_cc[None,:]] *= self.reductionTravel
                        Tmat[indices_G2_cc[:,None], indices_G2[None,:]] *= self.reductionTravel
                        
                        # Intra-state
                        if self.intra_reduction == True:
                            
                            Tmat[indices_G2_cc[:,None], indices_G2_cc[None,:]] *= self.reductionTravel

                            # Compensate for local part
                            indices_clusters_1G2_cc = self.cluster2_children[cc]

                            for dd in indices_clusters_1G2_cc:

                                indices_1G2_ccdd = self.cluster1_Tml_children[dd]
                                Tmat[indices_1G2_ccdd[:,None], indices_1G2_ccdd[None,:]] *= 1/self.reductionTravel

                else:
                    self.REDcounter2[cc] = 0
        
        
        return Tmat
    
    
    
    def apply_policy_2scales_vectorized (self,Ivec,Tmat):
        """"
        Update 2-scales policy status and modify transmission probabilities accordingly.

        Vectorized Implementation (like simple one) using information of previous policy status to
        only change T entries that change with the new policy status.
        """
        
        g1_clusters_before = (self.current_policy1 == 'g')
        r1_clusters_before = (self.current_policy1 == 'r')
        
        g2_clusters_before = (self.current_policy2 == 'g')
        r2_clusters_before = (self.current_policy2 == 'r')
        
        # Level 1 policy update
        NIc = np.array([ np.sum( Ivec[self.cluster1_children[cc]] ) for cc in range(self.nclusters1) ])
        
        change_G1toR1, change_R1toG1 = self.update_policy(NIc, self.current_policy1, self.GREENcounter1,
                                                          self.REDcounter1, self.GREENdelay, self.REDdelay,
                                                          self.GtoRthreshold, self.RtoGthreshold)
        
        g1_clusters_after = (self.current_policy1 == 'g')
        r1_clusters_after = (self.current_policy1 == 'r')
        
        # Level 2 policy update
        nReds = np.array([ np.sum(self.current_policy1[self.cluster2_children[cc]]=='r') for cc in range(self.nclusters2) ])
        
        change_G2toR2, change_R2toG2 = self.update_policy(nReds, self.current_policy2, self.GREENcounter2,
                                                          self.REDcounter2, self.GREENdelay, self.REDdelay,
                                                          1, 0)
        
        g2_clusters_after = (self.current_policy2 == 'g')
        r2_clusters_after = (self.current_policy2 == 'r')
        
        
        # Level 1: Compute indices of R and G Tml-entities
        indices_G1_before = self.compute_indices_lower_level(g1_clusters_before, self.nclusters1, self.cluster1_Tml_children)
        indices_G1_after = self.compute_indices_lower_level(g1_clusters_after, self.nclusters1, self.cluster1_Tml_children)
        
        indices_R1_before = self.compute_indices_lower_level(r1_clusters_before, self.nclusters1, self.cluster1_Tml_children)
        indices_R1_after = self.compute_indices_lower_level(r1_clusters_after, self.nclusters1, self.cluster1_Tml_children)
        
        indices_G1toR1 = self.compute_indices_lower_level(change_G1toR1, self.nclusters1, self.cluster1_Tml_children)
        indices_R1toG1 = self.compute_indices_lower_level(change_R1toG1, self.nclusters1, self.cluster1_Tml_children)
        
        # Level 2: Compute indices of R and G Tml-entities
        indices_G2_before = self.compute_indices_lower_level(g2_clusters_before, self.nclusters2, self.cluster2_Tml_children)
        indices_G2_after = self.compute_indices_lower_level(g2_clusters_after, self.nclusters2, self.cluster2_Tml_children)
        
        indices_R2_before = self.compute_indices_lower_level(r2_clusters_before, self.nclusters2, self.cluster2_Tml_children)
        indices_R2_after = self.compute_indices_lower_level(r2_clusters_after, self.nclusters2, self.cluster2_Tml_children)
        
        # Combined G1 and G2 (could use np.union1d too, but probably slower due to sorting)
        indices_G1G2_before = indices_G2_before[ np.isin(indices_G2_before,indices_G1_before) ]
        indices_G1G2_after = indices_G2_after[ np.isin(indices_G2_after,indices_G1_after) ]
        indices_GtoG = indices_G1G2_before[ np.isin(indices_G1G2_before,indices_G1G2_after) ]
        
        # Combined R = (R1 and G2) OR (G1 and R2) OR (R1 and R2) before/after
        indices_R_after = np.append(indices_R1_after, indices_R2_after[~np.isin(indices_R2_after,indices_R1_after)])
        indices_R_before = np.append(indices_R1_before, indices_R2_before[~np.isin(indices_R2_before,indices_R1_before)])
        
        # Regions that changed from GG to R(any) or from R(any) to GG
        indices_G1G2toR = indices_G1G2_before[ np.isin(indices_G1G2_before,indices_R_after) ]
        indices_RtoG1G2 = indices_G1G2_after[ np.isin(indices_G1G2_after,indices_R_before) ]
        
        
        
        # Indices of clusters level 1 (same as above, maybe can combine)
        indices_clusters_G1_before = np.arange(self.nclusters1)[g1_clusters_before]
        indices_clusters_G1_after = np.arange(self.nclusters1)[g1_clusters_after]
        
        indices_clusters_R1_before = np.arange(self.nclusters1)[r1_clusters_before]
        indices_clusters_R1_after = np.arange(self.nclusters1)[r1_clusters_after]
        
        indices_clusters_G1toR1 = np.arange(self.nclusters1)[change_G1toR1]
        indices_clusters_R1toG1 = np.arange(self.nclusters1)[change_R1toG1]
        
        indices_clusters_1G2_before = self.compute_indices_lower_level(g2_clusters_before, self.nclusters2, self.cluster2_children)
        indices_clusters_1G2_after = self.compute_indices_lower_level(g2_clusters_after, self.nclusters2, self.cluster2_children)
        
        indices_clusters_1R2_before = self.compute_indices_lower_level(r2_clusters_before, self.nclusters2, self.cluster2_children)
        indices_clusters_1R2_after = self.compute_indices_lower_level(r2_clusters_after, self.nclusters2, self.cluster2_children)
        
        indices_clusters_G1G2_before = indices_clusters_1G2_before[ np.isin(indices_clusters_1G2_before,indices_clusters_G1_before) ]
        indices_clusters_G1G2_after = indices_clusters_1G2_after[ np.isin(indices_clusters_1G2_after,indices_clusters_G1_after) ]
        
        indices_clusters_R_before = np.append(indices_clusters_R1_before, indices_clusters_1R2_before[~np.isin(indices_clusters_1R2_before,indices_clusters_R1_before)])
        indices_clusters_R_after = np.append(indices_clusters_R1_after, indices_clusters_1R2_after[~np.isin(indices_clusters_1R2_after,indices_clusters_R1_after)])
        
        indices_clusters_G1G2toR = indices_clusters_G1G2_before[ np.isin(indices_clusters_G1G2_before,indices_clusters_R_after) ]
        indices_clusters_RtoG1G2 = indices_clusters_G1G2_after[ np.isin(indices_clusters_G1G2_after,indices_clusters_R_before) ]
        
        indices_clusters_G1toG1 = np.arange(self.nclusters1)[g1_clusters_before & g1_clusters_after]
        
        
        
        

        # Update Tmat
        # Travel changes if region 1 is R_after (R_before), and region 2 is G1G1_before (G1G2_after)
        Tmat[indices_G1G2toR[:,None],indices_G1G2toR[None,:]] *= 1/self.reductionTravel
        Tmat[indices_G1G2toR[:,None],indices_GtoG[None,:]] *= 1/self.reductionTravel
        Tmat[indices_GtoG[:,None],indices_G1G2toR[None,:]] *= 1/self.reductionTravel
        
        Tmat[indices_RtoG1G2[:,None],indices_RtoG1G2[None,:]] *= self.reductionTravel
        Tmat[indices_RtoG1G2[:,None],indices_GtoG[None,:]] *= self.reductionTravel
        Tmat[indices_GtoG[:,None],indices_RtoG1G2[None,:]] *= self.reductionTravel
        
        # Diagonal compensation due to travel changes above.
        for aa in indices_clusters_G1G2toR:
            
            indices_G1G2toR_aa = self.cluster1_Tml_children[aa]
            Tmat[indices_G1G2toR_aa[:,None],indices_G1G2toR_aa[None,:]] *= self.reductionTravel
            
        for aa in indices_clusters_RtoG1G2:
            
            indices_RtoG1G2_aa = self.cluster1_Tml_children[aa]
            Tmat[indices_RtoG1G2_aa[:,None],indices_RtoG1G2_aa[None,:]] *= 1/self.reductionTravel
        

        # Level 1 - Local part
        for aa in indices_clusters_G1toR1:
            
            indices_G1toR1_aa = self.cluster1_Tml_children[aa]
            Tmat[indices_G1toR1_aa[:,None],indices_G1toR1_aa[None,:]] *= 1/self.reductionLocal
            
        for aa in indices_clusters_R1toG1:
            
            indices_R1toG1_aa = self.cluster1_Tml_children[aa]
            Tmat[indices_R1toG1_aa[:,None],indices_R1toG1_aa[None,:]] *= self.reductionLocal
        
        
        
        # Intra-level-2
        if self.intra_reduction == False:
            
            # Indices of G1 inside an R2
            indices_clusters_G2toR2 = np.arange(self.nclusters2)[g2_clusters_before & r2_clusters_after]
            indices_clusters_R2toG2 = np.arange(self.nclusters2)[r2_clusters_before & g2_clusters_after]
            
            indices_G1toG1 = indices_G1_before[ np.isin(indices_G1_before, indices_G1_after) ]
            
            # Remove intra-level-2 restrictions for G1s inside R2s
            for cc in indices_clusters_G2toR2:
                
                indices_G2toR2_cc = self.cluster2_Tml_children[cc]
                indices_G1G2toG1R2_cc = indices_G2toR2_cc[ np.isin(indices_G2toR2_cc, indices_G1toG1) ]
            
                Tmat[indices_G1G2toG1R2_cc[:,None],indices_G1G2toG1R2_cc[None,:]] *= self.reductionTravel
                
                # Compensate block diagonal
                indices_clusters_1G2to1R2_cc = self.cluster2_children[cc]
                indices_clusters_G1G2toG1R2_cc = indices_clusters_1G2to1R2_cc[ np.isin(indices_clusters_1G2to1R2_cc, indices_clusters_G1toG1) ]                
                for dd in indices_clusters_G1G2toG1R2_cc:

                    indices_G1G2toG1R2_ccdd = self.cluster1_Tml_children[dd]
                    Tmat[indices_G1G2toG1R2_ccdd[:,None],indices_G1G2toG1R2_ccdd[None,:]] *= 1/self.reductionTravel
                
                #Tmat[indices_G1G2toG1R2_cc,indices_G1G2toG1R2_cc] *= 1/self.reductionTravel
                
            for cc in indices_clusters_R2toG2:
                
                indices_R2toG2_cc = self.cluster2_Tml_children[cc]
                
                indices_G1R2toG1G2_cc = indices_R2toG2_cc[ np.isin(indices_R2toG2_cc, indices_G1toG1) ]
            
                Tmat[indices_G1R2toG1G2_cc[:,None],indices_G1R2toG1G2_cc[None,:]] *= 1/self.reductionTravel
                
                # Compensate block diagonal
                indices_clusters_1R2to1G2_cc = self.cluster2_children[cc]
                indices_clusters_G1R2toG1G2_cc = indices_clusters_1R2to1G2_cc[ np.isin(indices_clusters_1R2to1G2_cc, indices_clusters_G1toG1) ]                
                for dd in indices_clusters_G1R2toG1G2_cc:

                    indices_G1R2toG1G2_ccdd = self.cluster1_Tml_children[dd]
                    Tmat[indices_G1R2toG1G2_ccdd[:,None],indices_G1R2toG1G2_ccdd[None,:]] *= self.reductionTravel
                
             
            # Indices of R2 before or after
            indices_clusters_R2_after = np.arange(self.nclusters2)[r2_clusters_after]
            indices_clusters_R2_before = np.arange(self.nclusters2)[r2_clusters_before]
            
            # Remove intra-level-2 restrictions for R1s to G1s that stay within R2 after / G1 to R1 with R2 before
            for cc in indices_clusters_R2_after:
                
                indices_R2_after_cc = self.cluster2_Tml_children[cc]
                
                indices_R1toG1R2_cc = indices_R2_after_cc[ np.isin(indices_R2_after_cc, indices_R1toG1) ]
                indices_G1toG1R2_cc = indices_R2_after_cc[ np.isin(indices_R2_after_cc, indices_G1toG1) ]
                
                Tmat[indices_R1toG1R2_cc[:,None],indices_G1toG1R2_cc[None,:]] *= self.reductionTravel
                Tmat[indices_G1toG1R2_cc[:,None],indices_R1toG1R2_cc[None,:]] *= self.reductionTravel
                
                Tmat[indices_R1toG1R2_cc[:,None],indices_R1toG1R2_cc[None,:]] *= self.reductionTravel
                
                # Compensate block diagonal
                indices_clusters_1R2_after_cc = self.cluster2_children[cc]
                indices_clusters_R1toG1R2_cc = indices_clusters_1R2_after_cc[ np.isin(indices_clusters_1R2_after_cc, indices_clusters_R1toG1) ]                
                for dd in indices_clusters_R1toG1R2_cc:

                    indices_R1toG1R2_cc = self.cluster1_Tml_children[dd]
                    Tmat[indices_R1toG1R2_cc[:,None],indices_R1toG1R2_cc[None,:]] *= 1/self.reductionTravel
                
            for cc in indices_clusters_R2_before:
                
                indices_R2_before_cc = self.cluster2_Tml_children[cc]
                
                indices_G1R2toR1_cc = indices_R2_before_cc[ np.isin(indices_R2_before_cc, indices_G1toR1) ]
                indices_G1R2toG1_cc = indices_R2_before_cc[ np.isin(indices_R2_before_cc, indices_G1toG1) ]
                
                Tmat[indices_G1R2toG1_cc[:,None],indices_G1R2toR1_cc[None,:]] *= 1/self.reductionTravel
                Tmat[indices_G1R2toR1_cc[:,None],indices_G1R2toG1_cc[None,:]] *= 1/self.reductionTravel
                
                Tmat[indices_G1R2toR1_cc[:,None],indices_G1R2toR1_cc[None,:]] *= 1/self.reductionTravel
                
                # Compensate block diagonal
                indices_clusters_1R2_before_cc = self.cluster2_children[cc]
                indices_clusters_G1toR1R2_cc = indices_clusters_1R2_before_cc[ np.isin(indices_clusters_1R2_before_cc, indices_clusters_G1toR1) ]                
                for dd in indices_clusters_G1toR1R2_cc:

                    indices_G1toR1R2_cc = self.cluster1_Tml_children[dd]
                    Tmat[indices_G1toR1R2_cc[:,None],indices_G1toR1R2_cc[None,:]] *= self.reductionTravel
            
        
        return Tmat
    
    
    
    
    def apply_policy_3scales_simple (self,Ivec,Tmat_dummy):
        """"
        Update 3-scales policy status and modify transmission probabilities accordingly.

        Simple implementation where previous policy status is not used.
        """
        
        Tmat = self.Tmat0.copy()  # Tmat_dummy is not used, only included to make all functions have the same inputs
        
        # Level 1 policy update
        NIc = np.array([ np.sum( Ivec[self.cluster1_children[cc]] ) for cc in range(self.nclusters1) ])
        
        change_GtoR1, change_RtoG1 = self.update_policy(NIc, self.current_policy1, self.GREENcounter1,
                                                        self.REDcounter1, self.GREENdelay, self.REDdelay,
                                                        self.GtoRthreshold, self.RtoGthreshold)
        
        g1_clusters_after = (self.current_policy1 == 'g')
        r1_clusters_after = (self.current_policy1 == 'r')
        
        # Level 2 policy update
        nReds1 = np.array([ np.sum(self.current_policy1[self.cluster2_children[cc]]=='r') for cc in range(self.nclusters2) ])
        
        change_GtoR2, change_RtoG2 = self.update_policy(nReds1, self.current_policy2, self.GREENcounter2,
                                                        self.REDcounter2, self.GREENdelay, self.REDdelay,
                                                        1, 0)
        
        g2_clusters_after = (self.current_policy2 == 'g')
        r2_clusters_after = (self.current_policy2 == 'r')
        
        
        # Level 3 policy update
        nReds2 = np.array([ np.sum(self.current_policy2[self.cluster3_children[cc]]=='r') for cc in range(self.nclusters3) ])
        
        change_GtoR3, change_RtoG3 = self.update_policy(nReds2, self.current_policy3, self.GREENcounter3,
                                                        self.REDcounter3, self.GREENdelay, self.REDdelay,
                                                        1, 0)
        
        g3_clusters_after = (self.current_policy3 == 'g')
        r3_clusters_after = (self.current_policy3 == 'r')
        
        
        # Cluster indices of R and G entities of corresponding mpolX-level (X=1,2,3)
        indices_clusters_R3 = np.arange(self.nclusters3)[r3_clusters_after]
        indices_clusters_R2 = np.arange(self.nclusters2)[r2_clusters_after]
        indices_clusters_G2 = np.arange(self.nclusters2)[g2_clusters_after]
        indices_clusters_R1 = np.arange(self.nclusters1)[r1_clusters_after]
        indices_clusters_G1 = np.arange(self.nclusters1)[g1_clusters_after]
        
        # Level 3: Compute indices of R and G Tml-entities
        indices_R3 = self.compute_indices_lower_level(r3_clusters_after, self.nclusters3, self.cluster3_Tml_children)
        indices_G3 = self.compute_indices_lower_level(g3_clusters_after, self.nclusters3, self.cluster3_Tml_children)
        
        
        # Update Tmat
        #
        # Level 3
        Tmat[indices_R3[:,None],indices_G3[None,:]] *= 1/self.reductionTravel
        Tmat[indices_G3[:,None],indices_R3[None,:]] *= 1/self.reductionTravel

        Tmat[indices_R3[:,None],indices_R3[None,:]] *= 1/self.reductionTravel
        
        # Level 2
        for aa in indices_clusters_R3:
            
            indices_clusters_2R3_aa = self.cluster3_children[aa]  # level 2 clusters in level 3 cluster aa
            
            indices_clusters_G2R3_aa = indices_clusters_2R3_aa[ np.isin(indices_clusters_2R3_aa, indices_clusters_G2) ]
            indices_clusters_R2R3_aa = indices_clusters_2R3_aa[ np.isin(indices_clusters_2R3_aa, indices_clusters_R2) ]
            
            # Remove restrictions from G2R3 entities
            for bb in indices_clusters_G2R3_aa:
                
                indices_G2R3_aabb = self.cluster2_Tml_children[bb]
                
                Tmat[indices_G2R3_aabb[:,None],indices_G2R3_aabb[None,:]] *= self.reductionTravel
                
            
            # Level 1 (leave restrictions for R2R3 entities)
            for bb in indices_clusters_R2R3_aa:
                
                indices_clusters_1R2R3_aabb = self.cluster2_children[bb]  # level 1 clusters in level 2 cluster aabb
                
                indices_clusters_G1R2R3_aabb = indices_clusters_1R2R3_aabb[ np.isin(indices_clusters_1R2R3_aabb, indices_clusters_G1) ]
                indices_clusters_R1R2R3_aabb = indices_clusters_1R2R3_aabb[ np.isin(indices_clusters_1R2R3_aabb, indices_clusters_R1) ]
                
                # Remove restrictions from G1R2R3 entities
                for cc in indices_clusters_G1R2R3_aabb:
                
                    indices_G1R2R3_aabbcc = self.cluster1_Tml_children[cc]

                    Tmat[indices_G1R2R3_aabbcc[:,None],indices_G1R2R3_aabbcc[None,:]] *= self.reductionTravel
                
                # Remove travel and add local restrictions to R1R2R3 entities
                for cc in indices_clusters_R1R2R3_aabb:
                
                    indices_R1R2R3_aabbcc = self.cluster1_Tml_children[cc]

                    Tmat[indices_R1R2R3_aabbcc[:,None],indices_R1R2R3_aabbcc[None,:]] *= self.reductionTravel/self.reductionLocal
        
        
        return Tmat
