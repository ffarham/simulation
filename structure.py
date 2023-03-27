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
    ps = np.linspace(0,1,101)
    clique_sizes = [1,2,3,4,5,6,7,8,9,10]
    num_of_cliques = 100
    p_iterations = 100       
    

    csResults = dict()
    for clique_size in clique_sizes:
        logging.info("Clique size: " + str(clique_size))

        net_struct_results = dict()
        for p in ps:
            # logging.info("p: " + str(p))

            p_results_greedy, p_results_aug_greedy, p_results_degree, p_results_random = dict(), dict(), dict(), dict()
            p_results_greedy[0], p_results_aug_greedy[0], p_results_degree[0] , p_results_random[0] = 0, 0, 0, 0        # NOTE: this is assuming the initial belief vector is a zero vector

            for _ in range(p_iterations):
                G = Graph(num_of_cliques, clique_size, p)
                G.reShapeGraph()

                terminals, _ = G.network_structure()
                if p not in net_struct_results: net_struct_results[p] = 0
                net_struct_results[p] += len(terminals)

        for p in net_struct_results: 
            net_struct_results[p] /= p_iterations
        
        csResults[clique_size] = net_struct_results
    
    with open(dir_path+"structure.txt", 'w') as f:
        f.write(json.dumps(csResults))

if __name__ == '__main__':
    cProfile.run('main()')
