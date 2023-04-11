import sys
import logging
import json
import logging
import numpy as np
import os
import cProfile

from graph import Graph

def main():

    # NOTE: increase the limit to allow the DFS to run on large networks
    sys.setrecursionlimit(5000)
    logging.basicConfig(level=logging.INFO)

    # NOTE: initialise directory path to store results in
    dir_path = "./results/"
    assert os.path.exists(dir_path), "Directory to store results in is not defined"

    # NOTE: initialise settings
    ps = np.linspace(0,1,101)   # re-wiring probabilities
    clique_sizes = [1,2,3,4,5,6,7,8,9,10]
    num_of_cliques = 100
    p_iterations = 100  # number of iterations to re-generate the graph

    terminal_results, non_terminal_results = dict(), dict()
    for clique_size in clique_sizes:
        logging.info("Clique size: " + str(clique_size))

        terminal_struct, non_terminal_struct = dict(), dict()
        for p in ps:
            if p % 10 == 0: logging.info("Re-wiring probability: ", p)

            for _ in range(p_iterations):
                G = Graph(num_of_cliques, clique_size, p)
                G.reShapeGraph()

                terminals, non_terminals, _ = G.network_structure()
                if p not in terminal_struct: terminal_struct[p] = 0
                terminal_struct[p] += len(terminals)
                if p not in non_terminal_struct: non_terminal_struct[p] = 0
                non_terminal_struct[p] += len(non_terminals)

        for p in terminal_struct: 
            terminal_struct[p] /= p_iterations
            non_terminal_struct[p] /= p_iterations
        
        terminal_results[clique_size] = terminal_struct
        non_terminal_results[clique_size] = non_terminal_struct
    
    with open(dir_path+"structure_T.txt", 'w') as f:
        f.write(json.dumps(terminal_results))
        f.write('\n')
        
    with open(dir_path+"structure_nonT.txt", 'w') as f:
        f.write(json.dumps(non_terminal_results))
        f.write('\n')

if __name__ == '__main__':
    # cProfile.run('main()')
    main()
