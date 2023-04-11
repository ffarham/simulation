import os
import sys
import logging
import copy
import random
import json
import cProfile

from graph import Graph

def main():
    # NOTE: increase the limit to allow the DFS to run on large networks
    sys.setrecursionlimit(5000)
    logging.basicConfig(level=logging.INFO)

    dir_path = "./results/"
    assert os.path.exists(dir_path), "Directory to store results in is not defined"

    # NOTE: initialise settings
    p = 0.05
    num_of_cliques = 100
    clique_size = 5
    target_size = 30    # max size of initial active set

    p_iterations = 10    # average over different re-wiring probabilities
    k_iterations = 100    # average over different threshold instances in one particular rewiring instance

    p_results_greedy, p_results_aug_greedy, p_results_degree, p_results_random = dict(), dict(), dict(), dict()
    p_results_greedy[0], p_results_aug_greedy[0], p_results_degree[0] , p_results_random[0] = 0, 0, 0, 0,         # NOTE: this is assuming the initial belief vector is a zero vector
    for p_iter in range(p_iterations):
        logging.info(p_iter)
        G = Graph(num_of_cliques, clique_size, p)  
        G.reShapeGraph()

        k_results_greedy, k_results_aug_greedy, k_results_degree, k_results_random = dict(), dict(), dict(), dict()
        for _ in range(k_iterations):

            # reset candidates and vertices at the start of each k iteration as the values get poped
            candidates_A = copy.deepcopy(G.getCandidates())
            candidates_AandCB = copy.deepcopy(candidates_A)
            vertices = copy.deepcopy(G.getVertices())
            degree_ordered_nodes = copy.deepcopy(G.get_degree_ordered_nodes())

            GA_0, GAA_0, DA_0, RA_0 = [], [], [], []
            for k in range(1, target_size+1): 
         
                # hill-climbing
                if len(candidates_A) == 0: 
                    if k not in k_results_greedy: k_results_greedy[k] = 0
                    k_results_greedy[k] += G.sigma(GA_0)
                else:
                    (best_candidate_A, total_activations_A) = G.hill_climbing(GA_0, candidates_A)
                    candidates_A.remove(best_candidate_A)
                    GA_0.append(best_candidate_A)
                    if k not in k_results_greedy: k_results_greedy[k] = 0
                    k_results_greedy[k] += total_activations_A

                # augmented hill-climbing 
                if len(candidates_AandCB) == 0: 
                    if k not in k_results_aug_greedy: k_results_aug_greedy[k] = 0
                    k_results_aug_greedy[k] += G.sigma(GAA_0)
                else:
                    (best_candidate_AA, total_activations_AA) = G.hill_climbing_augmented(GAA_0, candidates_AandCB)
                    candidates_AandCB.remove(best_candidate_AA)
                    GAA_0.append(best_candidate_AA)
                    if k not in k_results_aug_greedy: k_results_aug_greedy[k] = 0
                    k_results_aug_greedy[k] += total_activations_AA

                # degree heuristic to determine next best node
                if len(degree_ordered_nodes) > 0:
                    # NOTE: adds randomness to selection
                    random.shuffle(degree_ordered_nodes[0])
                    degree_candidate = degree_ordered_nodes[0].pop(0)
                    if degree_ordered_nodes[0] == []: degree_ordered_nodes.pop(0)
                    DA_0.append(degree_candidate)      
                if k not in k_results_degree: k_results_degree[k] = 0
                k_results_degree[k] += G.sigma(DA_0)

                # random selection of initial active set
                randNode = random.choice(vertices)
                vertices.remove(randNode)
                RA_0.append(randNode)
                if k not in k_results_random: k_results_random[k] = 0
                k_results_random[k] += G.sigma(list(RA_0))
            
            # randomise the thresholds
            G.reset()

        # average the results over k iterations
        for k in k_results_greedy: 
            if k not in p_results_greedy: p_results_greedy[k] = 0
            p_results_greedy[k] += k_results_greedy[k] // k_iterations
            if k not in p_results_aug_greedy: p_results_aug_greedy[k] = 0
            p_results_aug_greedy[k] += k_results_aug_greedy[k] // k_iterations
            if k not in p_results_degree: p_results_degree[k] = 0
            p_results_degree[k] += k_results_degree[k] // k_iterations
            if k not in p_results_random: p_results_random[k] = 0
            p_results_random[k] += k_results_random[k] // k_iterations

    # average the results over p iterations
    for k in p_results_greedy: 
        p_results_greedy[k] = p_results_greedy[k] // p_iterations
        p_results_aug_greedy[k] = p_results_aug_greedy[k] // p_iterations
        p_results_degree[k] = p_results_degree[k] // p_iterations
        p_results_random[k] = p_results_random[k] // p_iterations
    
    # save results to file
    with open(dir_path + "simulation.txt", 'w') as f:
        f.write(json.dumps(p_results_greedy))
        f.write('\n')
        f.write(json.dumps(p_results_aug_greedy))
        f.write('\n')
        f.write(json.dumps(p_results_degree))
        f.write('\n')
        f.write(json.dumps(p_results_random))
        f.write('\n')

if __name__ == '__main__':
    # cProfile.run("main()")
    main()