# Transmission probability from USA data

The procedure to derive the transmission probability $T^{(0,0)}_{ij}$ (which is equivalent to $R^{(0,0)}_{ij}$ in the paper) is described in Ref. [PAPER].

All the necessary processes are implemented as functions in the file `usa_data_Tij.py`. The master function `compute_Tij_from_data` calls in order all the necessary functions to load and manipulate the data, and derive from it the $T^{(0,0)}$ matrix.

Some remarks:

- This file computes $T^{(0,0)}_{ij}$ assuming a `DEFAULT_BETA=1` value for $R_0$. This is why in simulations we initially multiply $T^{(0,0)}_{ij}$ by the `beta` value chosen in `parameters.py` (see constructor of `Dynamics` in `dynamics.py`).

- The values of `DEFAULT_AIR_PREFACTOR` and `DEFAULT_COMM_PREFACTOR` were used for simulations in Ref. [PAPER]. However, they can be changed to control the relative importance of airtravel and commuting with respect to local transmission.





