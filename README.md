# Multiscale Stochastic Simulations of Disease Spread

## Description

This repository contains the code and simulation notebooks used in the paper [PAPER].
It implements the discrete-time dynamics of an epidemiological SIS model subject to policy interventions.
The dynamics and the policy responses can be applied at multiple scales (e.g. people, counties, states...).

In the code, we use the USA as a working example and divide it into 5 levels: individuals, counties, states, divisions, and the whole country.
For the simulations we use publicly available data on the US population distribution from the US Census, flight data from the FAA, and commuting data from the US Census.


## Repository structure

### Documentation

The primary source of information is the paper [PAPER].

The `docs` folder contains some additional information on the model being simulated (`model.md`), the public sources used for the data (`data_sources.md`) and the transformations made to derive from it a matrix of probabilities (`data_to_Tmatrix.md`).


### Code

All the modules used are in the folder `src`.
Here's a brief description:

- `dynamics.py`: Contains the main `Dynamics` class where the system is initialized and the time evolution happens.

- `indices.py`: Contains functions to create the arrays of indices describing the hierarchy of partitions.

- `initial_conditions.py`: Contains the `Initial_Condition` class with modules to select the initial state of infections.

- `observables.py`: Contains the `Observables` class with modules to compute and save different quantities.

- `parameters.py`: Contains the default values of the parameters required to initialize the class objects.

- `policy.py`: Contains the `Policy` class with methods to change the policy status of different regions and implement restrictions by reducing/increasing the corresponding transmission probabilities.

- `usa_data_Tij.py`: Contains a number of functions to load and manipulate the USA data, and to use it to compute the probability matrix $T^{(0,0)}_{ij}$.



### Data

The data sources are given in `docs/data_sources.md`. The `data` folder contains:

- Commuting and population data (`UScensus_data`), airflight data (`FAA_data`), and other US data (`Simplemaps_data` and `StatesFIPS_data`) used to generate the $T_{ij}$ matrix used for simulations.

- Geodata of different regions (`geodata`) used for plots.



### Notebooks

The `notebooks` folder contains the notebooks used to generate the data and figures used in the paper [PAPER]. These notebooks exemplify how to use the source code to run simulations (specifically, `fig3.ipynb` and `fig4.ipynb`). In principle, the notebooks can be run cell by cell, unless otherwise stated. Some `.shp` geodata files might be missing from the repository, but they can be downloaded from the links specified in `docs/data_sources.md`. The notebooks include:

- `fig1.ipynb`: Plots to illustrate self-similar framework. (No simulations)
- `fig2.ipynb`: Plots to illustrate multiscale $T_{ij}$ framework. (No simulations)
- `fig3.ipynb`: Simulations and plots showing elimination dynamics for different policies.
- `fig4.ipynb`: Simulations and plots showing the cost and effectiveness of different policies for constant importation rate.


## Dependencies

The dependencies are listed in `requirements.txt`. This is the full list of packages in the enviroment that I'm using to run the code. You can create a Python virtual environment and then you should be able to install them using `pip install -r requirements.txt`.

However, all packages in `requirements.txt` are almost certainly not needed (see packages imported in the source code). The most important package required to run the code is `numpy`, `scipy` and `pandas`, and to run notebooks you need also `geopandas`, and `matplotlib`.


## References

[PAPER] _Self-similarity in pandemic spread and fractal containment policies,_ A. F. Siegenfeld, A. Piñeiro Orioli, R. Na, B. Elias, Y. Bar-Yam.





