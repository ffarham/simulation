# Attacking Social Networks

This repository contains the scripts used in the practical component of my third year project. Stochastic Block Model with parameters number of cliques, clique size and re-wiring probability is used to generate a social network that has allows different opinions of agents to persist as well as intermix at the same time. the settings used for 100 clique are: clique size of 5 and a re-wiring probability of 0.05.

To perform your own simulations, initialise `results` directory to store simulation results. The only prerequisite is to have [Python3](https://www.python.org/downloads/) installed in your computer.

- `graph.py`: Graph class to generate a Stochastic Block Model with the following parameters: number of cliques, clique size and re-wiring probability.
- `simulation.py`: script to run the main simulation on specified settings. Results stored in *results/simulation.txt*.
- `structure.py`: script to determine the structure (num. of terminals and non-terminals) of the networks generated across varying clique size and re-wiring probability for some fixed number of cliques. Results are stored seperately in *results/structure_T.txt* for terminals and in *results/structure_nonT.txt* for non-terminals.
- `plot_simualtion.py`: plots the initial active set against the resulting number of activations from the script *simulation.py*.
- `plot_structure.py`: plots the re-wiring probability + clique size against the number of terminals/non-terminals from the script *structure.py*. Need to specify which file to plot at the start. 
- `plot_model_features.py`: plots the the number of terminals and non-terminals for a particular clique size against varying re-wiring probability. Linked with the results from the script *structure.py*.