# Multi-Scale SIS Dynamics

## Description

This repository contains the code and simulation notebooks used in the paper [...].
It implements the discrete-time dynamics of an epidemiological SIS model subject to policy interventions.
The dynamics and the policy responses can be applied at multiple scales (e.g. people, counties, states...).

In the paper and the notebooks provided we use the USA as a working example, and divide it into 5 levels: individuals, counties, states, divisions, and the whole country.
For the simulations we use publicly available data on the US population distribution from the US Census, flight data from the FAA, and commuting data from the US Census.


## Repository structure

### Documentation

The `docs` folder contains a more detailed description of the model being simulated (`model.md`), the public sources used for the data (`data_sources.md`) and the transformations made to derive from it a matrix of probabilities (`data_manipulations.md`).


### Code

All the modules used are in the folder `src`.
Here's a brief description:

- `dynamics.py`: Contains the main `Dynamics` class where the system is initialized and the time evolution happens.

- `indices.py`: Contains functions to create the arrays of indices describing the hierarchy of partitions.

- `initial_conditions.py`: Contains the `Initial_Condition` class with modules to select the initial state of intections.

- `observables.py`: Contains the `Observables` class with modules to compute and save different quantities.

- `parameters.py`: Contains the default values of the parameters required to initialize the class objects.

- `policy.py`: Contains the `Policy` class with modules to change the policy status of different regions and implement restrictions by reducing/increasing the corresponding transmission probabilities.

- `usa_data_Tij.py`: Contains a number of functions to load and manipulate the USA data, and to use it to compute the probability matrix $T^{(m,n)}_{ij}$.



### Data

The `data` folder contains:

To-Do

(Which data should we include?? All??)



### Notebooks

The `notebooks` folder contains some usage examples, as well as the notebooks used to generate the data contained in the paper [...]. Specifically:

- To-Do





## Usage

See notebooks description above and in `notebooks` folder.




## Dependencies and Prerequisites

To-Do





## References

To-Do

(References to your paper or other related works)





