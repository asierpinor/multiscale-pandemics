# Theory: Multiscale SIS dynamics

For more details on the  model we refer to Ref. [PAPER].

## Multi-scale structure and indices

We consider a partition of the system into a nested hierarchy of regions of different sizes. Per default we use the USA as an example with individuals ($n=0$), counties ($n=1$), states ($n=2$), divisions ($n=3$), USA ($n=4$). In the paper [PAPER], we don't use divisions and instead consider $n=3$ to be the USA.

The model requires specifying two levels $(m,n)$, which correspond to the parameters `ml` and `nl` in `parameters.py`. This means that we are using a representation describing the number of infected $n$-level regions within each $m$-level region.
For example, $(1,0)$ means that we keep track of the number of infected individuals in each county (but not on the specific individuals).

It is useful to define $N^{(n)}$ the total number of $n$-level regions in the system, and $N^{(m,n)}_i$ as the number of $n$-level regions within the $m$-level region $i$ (note $N^{(n)}\equiv N^{(4,n)}$). These quantities  correspond to the dictionaries `Ntots` and `N`, respectively, defined in the `set_N_arrays` function inside `indices.py`)


In the code we also define the dictionaries `iparent` and `ichild` using functions inside `indices.py`. The parent function $\ell: \mathbb{N} \rightarrow \mathbb{N}$ maps an index $i$ of an $n$-level region to an index $\ell(i)$ of the $(n+1)$-level region that $i$ belongs to. We will often write $\ell_i \equiv \ell(i)$.
Similarly, the child function $c$ maps an index $i$ of an $n$-level region to the set of indices $c(i)=\{c_0(i),c_1(i),\ldots\}$ of the $(n-1)$-level regions that belong to $i$.

_Note:_ The indices ($i, j, \ldots$) can refer to any level and it should be clear from the context which level is meant.




## Stochastic SIS Model

_Note:_ For more details on the underlying model expressed in terms of probability vectors see Ref.~[PAPER]. Here we mostly focus on the numerical implementation.

We simulate a stochastic discrete-time SIS model for any choice of representation level $(m,n)$. 
At each time $t$ the system is described by a vector of infections $\vec{I}^{(m,n)}_{k}$ (`Ivec` in `dynamics.py`).
When $m=n$, we have $I^{(m,m)}_{k}=1$ if the $m$-level entity $k$ is infected, and $I^{(m,m)}_{k}=0$ otherwise.
When $m>n$, then $I^{(m,n)}_{k}$ is a natural number and counts the number of infected $n$-level entities that are within the $m$-level entity $k$.
For example, $I^{(1,0)}_{k}$ gives the number of infected people in county $k$.

The spread of the disease is determined by a transmission probability matrix $T^{(m,n)}$ (`Tmat` in `dynamics.py` and the rest of the code) and the probability of recovery $p_R$ within one time step (`pR` in `parameters.py`).
In the code, when $m=n$, $T^{(m,m)}_{ij}$ denotes the probability that $m$-level entity $j$ will infect $m$-level entity $i$, if $j$ is infected.
When $m>n$, $T^{(m,n)}_{ij}$ denotes the probability that an infected $n$-level entity belonging to $m$-level entity $j$ will infect an $n$-level entity belonging to $m$-level entity $i$.
Note that in the paper [PAPER] we use instead the generation matrix $R^{(m,n)}$, which can in principle be computed from $T^{(m,n)}$.
Note also that all $n$-level entities within an $m$-level entity are assumed to have the same probabilities of infecting or being infected.
For example, $T^{(1,0)}_{ij}$ denotes the probability that an infected person in county $j$ will infect a person in county $i$.

We use $T^{(m,n)}$ and $p_R$ to compute the probability of being infected in the next time step. For this three models are implemented: `'SIS'`, `'SISreI'`, and `'SISlinear'` (see `model` parameter in `parameters.py`). The default and the one used in the paper [PAPER] is `'SISreI'`. They differ slightly in the way this probability is calculated, but in the limit of low infections they are practically indistinguishable.

Basically, `'SISreI'` assumes that $p_R=1$ and that an infected entity can recover and be re-infected in one time step for simplicity; in other words, it's as if at each time step all people start susceptible and can be infected with some probability.
Using $T^{(m,n)}$ we calculate the probability $P^{(m,n)}_i(t+1)$ that every single $n$-level entity belonging to $m$-level entity $i$ has of getting infected in the next time step (this corresponds to the variable `probI` in `dynamics.py`).
When $m=n$, this probability is given by
$$
P^{(m,m)}_i(t+1) = 1 - \prod_{j}(1-T^{(m,m)}_{ij}\,I^{(m,m)}_j(t)),
$$
where the index $j$ runs over all entities.
The product term is just the probability of not getting infected, which is a product of the probabilities of not getting infected by any infected entity.
When $m>n$, this probability can be calculated as (this expression actually holds for $m=n$ too in the stochastic model)
$$
P^{(m,n)}_i(t+1) = 1 - \prod_{j}(1-T^{(m,n)}_{ij})^{I^{(m,n)}_j(t)}.
$$
This expression is essentially the same as the previous one, but it uses the homogeneity within each $m$-entity to simplify the calculation.
Note that we have included the possibility of self-infection, which is negligible anyway.

In the `'SIS'` model, we include the possibility of having $p_R<1$ by combining it with the above probability as $P_i^{(m,n)} \rightarrow (1-I_i^{(m,n)}/N_i^{(m,n)})P_i^{(m,n)} + (I_i^{(m,n)}/N_i^{(m,n)}) (1-p_R)$, omitting the $t$ dependence for simplicity.
In the `'SISlinear'` model, we simplify the computation of the probability by making a linear approximation as $P^{(m,n)}_i = \sum_{j} T^{(m,n)}_{ij} I^{(m,n)}_j$.





## Sampling infections

We use the above probability $P^{(m,n)}_i(t+1)$ to compute the vector of infections $\vec{I}^{(m,n)}(t+1)$ at the next time step.
We do this by sampling from a probability distribution, either binomial (`'binomial'`), Poisson (`'poisson'`), or negative binomial (`'neg_binomial'`), see `distrib` parameter in `parameters.py`.
In the case of the binomial distribution, we obtain $I^{(m,n)}_i(t+1)$ by drawing numbers from $\text{Binom}(n=N^{(m,n)}_i,p=P^{(m,n)}_i(t+1))$. If the number of infections is small compared to $N^{(m,n)}_i$, this can be approximated by a Poisson distribution, $\text{Poisson}(\lambda=N^{(m,n)}_i P^{(m,n)}_i(t+1))$.

The negative binomial distribution includes the possibility of superspreading [CITE].
We parametrize the distribution, $\text{NegBin}(\mu,\kappa)$, in terms of the mean $\mu$ and a superspreading parameter $\kappa$, which corresponds to `sspread` in `parameters.py`.
This parametrization looks like
$$
    \text{NB}(x;\mu,\kappa) = \begin{pmatrix} x+r-1 \\ x \end{pmatrix} p^r (1-p)^x ,
$$
with $r = \mu/\kappa$ and $p = 1 / (1+\kappa)$.
Note that the variance of this distribution is $\sigma^2 = \mu (1+ \kappa)$.
This parametrization ensures that adding up negative binomials with the same $\kappa$ but different $\mu_i$ results in another negative binomial with mean $\sum_i \mu_i$ and same $\kappa$.






# Multiscale policies

We consider policies (see `policy.py`) applied at either one (`'1scale'`), two (`'2scales'`), or three scales (`'3scales'`), as selected by the `PolicyChoice` parameter in `parameters.py`. These policies can be applied at any level $m_1$, $m_2$, $m_3$ as specified by the parameters `mpol1`, `mpol2` and `mpol3`, see `parameters.py`.

**Note:** In the paper [PAPER] we refer to simulations with `'1scale'` as being 2-scale, and simulations with `'2scales` as being 3-scale. The reason is that in the paper we consider the individual scale as an additional scale, which in simulations is implicitly taken into account by the choice of local transmission probability $R_0$ (`beta` in `parameters.py`). We will stick to the code's notation in the context of the code's documentation.

All policies are based on categorizing regions of the corresponding level as either red ($R$) or green ($G$), and applying two types of restrictions: local restrictions within a region (lockdowns, masks, etc.), and travel restrictions between regions (travel reductions, quarantine, etc.). The general rules for single and pairs of regions are as follows:

- $G\ \longrightarrow$ no local restrictions
- $R\ \longrightarrow$ local restrictions
- $GG\ \longrightarrow$ no travel restrictions
- $GR,RG,RR\ \longrightarrow$ travel restrictions

Here, $GR$ means a region is $G$ and the other is $R$.
The policy status of each region is stored in the `current_policy` arrays inside the class `Policy`.

To determine the $R/G$ status of a region and how restrictions are implemented we use the following parameters (see `parameters.py`): $h_{G\rightarrow R}$ (`GtoRthreshold`), $h_{R\rightarrow G}$ (`RtoGthreshold`), $d_G$ (`GREENdelay`), $d_R$ (`REDdelay`), $r_L$ (`reductionLocal`) and $r_T$ (`reductionTravel`).

For simulations with **1-scale** policies we proceed as follows:

- Let region $i$ have policy $G$. If the number of infections rises above a threshold, $N^I_i \geq h_{G\rightarrow R}$, for a time duration $\Delta t > d_G$, the region changes its policy to $R$.
    
- Let region $i$ have policy $R$. If the number of infections sinks below a threshold, $N^I_i \leq h_{R\rightarrow G}$ for a time duration $\Delta t > d_R$, the region changes its policy to $G$.

- If a pair of regions changes policy from $GG$ to $RG$, $GR$, or $RR$, we rescale the transmission matrix by $T^{(0)}_{ij} \rightarrow T^{(0)}_{ij}/r_L$ when $i,j$ belong to the same region, and $T^{(0)}_{ij} \rightarrow T^{(0)}_{ij}/r_T$ when $i,j$ belong to different regions.

- If a pair of regions changes policy from $RG$, $GR$, or $RR$ to $GG$ we rescale the transmission matrix back by $T^{(0)}_{ij} \rightarrow T^{(0)}_{ij} r_L$ when $i,j$ belong to the same region, and $T^{(0)}_{ij} \rightarrow T^{(0)}_{ij} r_T$ when $i,j$ belong to different regions.

In the simulations done for the paper [PAPER] we made some simplifying assumptions to reduce the number of parameters:

1. Note that $d_G$ ($d_R$) and $h_{G\rightarrow R}$ ($h_{R\rightarrow G}$) have a similar role. If the policy is kept constant, we can expect exponential growth or decline with some approximately constant rate. This means that waiting a time $d_G$ ($d_R$) is to some extent similar to setting a higher (lower) threshold $h_{G\rightarrow R}$ ($h_{R\rightarrow G}$). Thus, we set $d_G=0$ and vary $h_{G\rightarrow R}$. 

2. Since it doesn't make sense to lift restrictions unless we know that the exponential growth isn't going to start immediately again, we assumed that a region doesn't go back to $G$ until it reaches zero cases and that it then changes immediately, i.e.~$h_{R\rightarrow G} = 0$, $d_R=0$.


For simulations with 2-scale and 3-scale policies we apply the above procedure self-similarly. The main difference is how the local restrictions are applied in the case of larger-scale regions. As an example, let's consider a 2-scale policy at the county (`mpol=1`) and state (`mpol=2`) levels. When a state turns $R$ we reduce the inter-state transmission probability by $r_T$, and the intra-state county-to-county transmission probability by $r_T$ as well, instead of $r_L$ (this is because it corresponds to travel restrictions for the constituent counties). However, the intra-county transmission probability is determined by the given county's policy. If the county is $G$ then the intra-county transmission remains unaffected, but if it's $R$ we reduce it by $r_L$.

The second important difference is that higher-level entities change policy to $R$ as soon as one of its lower-level entities changes policy to $R$. In the previous example, a state turns $R$ as soon as one of its counties turns $R$. This implies that a county's transmission probability to other counties is always determined by the policy of the states they belong to.

