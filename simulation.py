import random
import numpy as np
import logging
import sys
import time
import json

"""
    HRCGraph - class for Homophilic Relaxed Caveman Graph
"""
class HRCGraph:

    """ constructor - generate a graph
        :param l - number of cliques
        :param s - size of each clique
        :param p - rewiring probability
    """
    def __init__(self, l, s, p):

        self.n = l * s
        self.m = s * (s-1) * s * l
        self.initial_vertices = [i for i in range(1, (l*s) + 1)]
        self.vertices = [i for i in range(1, (l*s) + 1)]
        self.edges = dict()
        self.in_degree = dict()
        self.degree_heuristics = []
        self.thresholds = dict()
        self.initBeliefs = dict()
        self.convergentBeliefs = dict()
        self.split = False
        self.reShape = False
        self.convergentT = None
        self.mapSVtoVs = dict()
        self.mapSVtoSs = dict()
        self.terminating_groups = []
        self.terminals = []
        self.non_terminals = []

        # create edge set for clique with self loops
        for i in range(1, len(self.vertices)+1):
            a = ((i-1) // s) * s if i % s == 0 else (i // s) * s
            b = a + s
            N = dict()
            for j in range(a+1, b+1):
                N[j] = random.uniform(0.01, 1.0)

            # normalise the weights
            total = sum([y for _,y in N.items()])                
            for m in N:
                N[m] = round(N[m]/total, 4)
            
            # ensure the weights sum to 1
            total = sum([y for _,y in N.items()])
            if total != 1:
                e = 1 - total
                N[i] += e
            
            self.edges[i] = N

        for node in self.vertices:
            self.thresholds[node] = random.uniform(0,1)
            self.initBeliefs[node] = 0 # random.uniform(0,1)
        
        # rewire the edges
        for u in self.edges:
            queue = list(self.edges[u].keys())
            while(len(queue) > 0):
                v = queue.pop(0)
                if v not in self.edges[u]: continue

                # do not rewire the self-loops
                if u == v: continue
                
                x = random.randint(1, len(self.vertices))
                while(x == v):
                    x = random.randint(1, len(self.vertices))

                # rewire with probability p
                if random.uniform(0,1) <= p:
                    if x in self.edges[u]: self.edges[u][x] += self.edges[u][v]
                    else: self.edges[u][x] = self.edges[u][v]

                    del self.edges[u][v]
        
        # determine the in-degree of each node
        for u in self.edges:
            for v in self.edges[u]:
                if v not in self.in_degree: self.in_degree[v] = 0
                self.in_degree[v] += 1
        sorted_degrees = sorted(self.in_degree.items(), key=lambda item: item[1], reverse=True)
        temp = []
        prev = None
        for item in sorted_degrees:
            if prev is None: temp.append(item[0])
            else:
                if item[1] != prev: 
                    self.degree_heuristics.append(temp)
                    temp = []
                temp.append(item[0])
            prev = item[1]
        self.degree_heuristics.append(temp)


    """ transpose - method to reverse the edges in the graph
        - creates a new dictionary of edges
    """
    def __transpose(self):
        if len(self.edges) == 0:
            logging.warning("found empty edge set when atttempting to transpose")
            return 

        newE = {}
        for i in self.edges:
            for j in self.edges[i]:  
                if j not in newE: newE[j] = {}

                newE[j][i] = self.edges[i][j]

        self.edges = newE

    
    """ dfs - perform DFS to determine minimal closed groups
    """
    def __dfs(self, u, visited, currentMCG):
        visited.add(u)
        currentMCG.append(u)
        for v in self.edges[u]:
            if v not in visited: self.__dfs(v, visited, currentMCG)
    

    """ dfsNumbering - number the nodes in the network by adding them to the stack
        - used to get the order of nodes to perform the DFS that determines minimal closed groups
    """
    def __dfsNumbering(self, u, visited, stack):
        visited.add(u)
        for v in self.edges[u]:
            if v not in visited: self.__dfsNumbering(v, visited, stack)
        stack = stack.append(u)


    """ getMCGs - determines the minimal closed groups in the graph with no outgoing edges
        - Kosaraju's algorithm: first perform DFS numbering on nodes and use that numbering as an order of nodes to traverse using DFS, each traversal traverses a minimal closed group
    """
    def __getMCGs(self):
        stack = []
        visited = set()

        for u in range(1, len(self.vertices) + 1):
            if u not in visited: self.__dfsNumbering(u, visited, stack)
        
        self.__transpose()

        MCGs = []
        visited = set()
        while len(stack) > 0:
            currentMCG = []
            u = stack.pop()
            if u not in visited:
                self.__dfs(u, visited, currentMCG)
                MCGs.append(currentMCG)

        self.__transpose()
        return MCGs


    """ splitGraph - determines the terminating groups and non-terminals
    """
    def __splitGraph(self):
        if self.split: return

        # determine the terminating groups: minimal closed groups with no outgoing edges
        MCGs = self.__getMCGs()
        for group in MCGs:
            flag = True
            for u in group:
                for v in self.edges[u]:
                    if v not in group:
                        flag = False
                        break
                if not flag: break

            if flag: self.terminating_groups.append(group)
            else: self.non_terminals += group
        
        self.split = True

        
    """ contraction - contract the graph on terminating group (minimal closed groups with no outgoing edges)
        :param tg - termianting group to contract
    """
    def __contract(self, tg):

        # remove the edges within the minimal closed group
        for u in tg:
            for v in tg:
                if v in self.edges[u]: del self.edges[u][v]
        
        # create supervertex
        sv = self.n + 1
        self.n += 1

        # add self-loop to the supervertex
        self.edges[sv] = dict()
        self.edges[sv][sv] = 1

        # rewire any edge connected to minimal closed group to the supervertex
        for u in self.non_terminals:
            for v in tg:
                # NOTE: by property of TG, there are no outgoing edges in the minimal closed group
                # incoming edges to minimal closed group: (u,v)
                if v in self.edges[u]:
                    # check if an edge to supervertex already exists
                    if sv not in self.edges[u]: self.edges[u][sv] = 0
                    
                    self.edges[u][sv] += self.edges[u][v]

                    del self.edges[u][v]          

        # clean up the edge set
        cuNodes = list(self.edges.keys())
        for node in cuNodes:
            if not self.edges[node]: del self.edges[node]

        # remove each node in the minimal closed group from the vertex set
        for node in tg:
            self.vertices.remove(node)

        # add supervertex to the vertices
        self.vertices.append(sv)

        return sv

    
    """ get_influence_vector - function to calculate the influence vector of the given terminating group
        :param tg - terminating group to calculate the infleunce vector of
        :return the influence vector
    """
    def __getInfluenceVector(self, tg):
        # calculate the consesnus reached by each terminating group
        tg.sort()
        beliefVector = np.array([self.initBeliefs[i] for i in tg])
        T = []
        for i in range(len(tg)):
            temp = [0] * len(tg)
            T.append(temp)

        for i in range(len(tg)):
            for j in range(len(tg)):
                if tg[j] in self.edges[tg[i]]: T[i][j] = self.edges[tg[i]][tg[j]]
        
        # ensure matrix T is row stochastic
        for i, row in enumerate(T):
            total = sum(row)
            for j, element in enumerate(row):
                T[i][j] = round(element / total, 4)
            newTotal = sum(T[i])
            if newTotal != 1:
                e = 1 - newTotal
                T[i][i] += e
        
        # get eigen vector and corresponding eigen value
        T = np.array(T)
        eigenValues, eigenVectors = np.linalg.eig(T.T)

        # find eigen vector corresponding to eigen value 1
        index = np.where(np.isclose(eigenValues, 1))[0][0]
        leftEigenVector = eigenVectors[:, index].real

        # Normalize the left eigenvector
        leftEigenVector = leftEigenVector / np.sum(leftEigenVector)
        
        return leftEigenVector


    """ reShape - re-shape the graph to fascilitate computing convergence
        :return the non-terminals 
    """
    def reShapeGraph(self):
        # if graph is already reshaped, return the determined non-terminals
        if self.reShape: return 

        # split the graph into terminating groups and non terminals if required
        if not self.split:
            self.__splitGraph()

        # contract each terminating group 
        for group in self.terminating_groups:
            s = self.__getInfluenceVector(group)
            sv = self.__contract(group)
            self.terminals.append(sv)
            self.mapSVtoVs[sv] = group
            self.mapSVtoSs[sv] = s

        self.reShape = True

        # converge the interaction matrix
        # define a matrix of co-efficients
        A = []
        for i in range(len(self.vertices)):
            temp = [0] * len(self.vertices)
            A.append(temp)
        
        for i, u in enumerate(self.vertices):
            for j, v in enumerate(self.vertices):
                if v in self.edges[u]: A[i][j] = self.edges[u][v]
        A = np.array(A)

        # converge A to some epsilon accuracy
        # TODO: test the cimulation on some fixed num. of iterations for all instances of A
        counter, epsilon = 0, 0.0001
        prev = A
        A = np.linalg.matrix_power(A, 2)
        while (np.sum(np.abs(A - prev)) > epsilon):
            if counter > 100: 
                logging.warning("counter limit hit during computing convergence")
                break
            prev = A
            A = np.linalg.matrix_power(A, 2)
            counter += 1
        
        self.convergentT = A

    """ hill_climbing - determine node that increases the convergent belief the most for each node in the network
        TODO: pick highest convergnce after highest activation
    """
    def hill_climbing(self, A_0, candidates):
        if not self.reShape: self.reShapeGraph()

        # TODO: randomise of convergent belief is also the same
        bestCandidate = None
        for c in candidates:
            total_activations = 0
            temp_activation = A_0 + [c]

            # determine convergent belief of all supervertices
            for sv in self.terminals:
                s = self.mapSVtoSs[sv]
                nodes = self.mapSVtoVs[sv]
                # create belief vector of desired nodes + any activations
                p = [self.initBeliefs[x] if x not in temp_activation else 1 for x in nodes]

                consensus = np.dot(s, p)
                self.convergentBeliefs[sv] = consensus
                # determine the nodes activated in the respective terminating group
                for node in nodes:
                    if consensus >= self.thresholds[node]: total_activations += 1
            
            # if all nodes belong to some minimal closed group
            if len(self.non_terminals) == 0: 
                value = sum(list(self.convergentBeliefs.values()))
                if bestCandidate is None:
                    bestCandidate = (c, value, total_activations)
                else:
                    if total_activations > bestCandidate[2]:
                        bestCandidate = (c, value, total_activations)
                    elif total_activations == bestCandidate[2]:
                        bestCandidate = bestCandidate if bestCandidate[1] > value else (c, value, total_activations)

            else:
                # still need to calculate convergent belief of each non terminals
                A = self.convergentT        
                b = np.array([self.initBeliefs[x] if x in self.initBeliefs else self.convergentBeliefs[x] for x in self.vertices])
                
                # calculate the vector of convergent beliefs in G'
                res = A.dot(b.T)
                value = sum(res)

                 # count the number of vertices in R that get activated
                for i, r in enumerate(self.non_terminals):
                    if res[i] >= self.thresholds[r]: total_activations += 1

                if bestCandidate is None:
                    bestCandidate = (c, value, total_activations)
                else:
                    if total_activations > bestCandidate[2]:
                        bestCandidate = (c, value, total_activations)
                    elif total_activations == bestCandidate[2]:
                        bestCandidate = bestCandidate if bestCandidate[1] > value else (c, value, total_activations)
        
        return bestCandidate

    """ sigma - influence function
        :param A_0 - initial set of active nodes
        :return the total number of active nodes at convergence
    """
    def sigma(self, A_0):
        
        if not self.reShape: self.reShapeGraph()

        total_activations = 0

        # determine convergent belief of all supervertices
        for sv in self.terminals:
            s = self.mapSVtoSs[sv]
            nodes = self.mapSVtoVs[sv]
            # create belief vector of desired nodes + any activations
            p = [self.initBeliefs[x] if x not in A_0 else 1 for x in nodes]

            consensus = np.dot(s, p)
            self.convergentBeliefs[sv] = consensus
            # determine the nodes activated in the respective terminating group
            for node in nodes:
                if consensus >= self.thresholds[node]: total_activations += 1
        
        # if all nodes belong to some minimal closed group
        if len(self.non_terminals) == 0: return total_activations

        # still need to calculate convergent belief of each non terminals
        A = self.convergentT        
        b = np.array([self.initBeliefs[x] if x in self.initBeliefs else self.convergentBeliefs[x] for x in self.vertices])
        
        # calculate the vector of convergent beliefs in G'
        c = A.dot(b.T)

        # count the number of vertices in R that get activated
        for i, r in enumerate(self.non_terminals):
            if c[i] >= self.thresholds[r]: total_activations += 1
        
        return total_activations
    
    """ updateThresholds - method to handle threhold randomisation 
        - also calculates the new subset of nodes that get activated with new threshold values
        :param A_0 - the new initial active set of nodes
    """
    def updateThresholds(self):
        
        for node in self.thresholds:
            self.thresholds[node] = random.uniform(0,1)


    """ getCandidates - get the candidates to activate
        :return the nodes in the terminating groups
    """
    def getCandidates(self):
        if not self.split:
            self.splitGraph()
        
        candidates = []
        for group in self.terminating_groups:
            candidates += group
        
        return candidates


    """ getTopDegreeHeuristic - get the top nodes with maximum degree
        :param i - top i nodes
    """
    def getTopDegreeHeuristic(self):
        
        if self.degree_heuristics == []: return []

        random.shuffle(self.degree_heuristics[0])

        output = self.degree_heuristics[0].pop(0)
        if self.degree_heuristics[0] == []: self.degree_heuristics.pop(0)

        return output

    """ getVertices - return all vertices in the initial network
    """
    def getVertices(self):
        return self.initial_vertices

    def __str__(self):
        return "Vertices: " + str(self.vertices) + "\nEdges: " + str(self.edges)


def main():
    sys.setrecursionlimit(5000)
    logging.basicConfig(level=logging.INFO)

    ps = [0, 0.2, 0.4, 0.6, 0.8, 1]
    clique_size = 5
    num_of_cliques = 20
    targetSize = 30
    iterations = 10
    
    for p in ps:
        # store results of greedy, degree and random selections
        resultsG = dict()
        resultsD = dict()
        resultsR = dict()
        # NOTE: this is assuming the initial belief vector is a zero vector
        resultsG[0] , resultsD[0] , resultsR[0] = 0, 0, 0

        # sizes of initial active sets to go through
        for _ in range(iterations):
            G = HRCGraph(num_of_cliques, clique_size, p)
            G.reShapeGraph()
            candidates = G.getCandidates()
            vertices = G.getVertices()

            GA_0, DA_0, RA_0 = [], [], []
            prev_activations_greedy = None
            for k in range(1, targetSize+1):  

                # hill climbing to determine next best node 
                if len(candidates) > 0:
                    best_candidateG, _, total_activations = G.hill_climbing(GA_0, candidates)
                    candidates.remove(best_candidateG)
                    GA_0.append(best_candidateG)
                    if k not in resultsG: resultsG[k] = 0
                    resultsG[k] += total_activations
                    prev_activations_greedy = total_activations
                else:
                    resultsG[k] += prev_activations_greedy                    

                # degree heuristic to determine next best node
                best_candidateD = G.getTopDegreeHeuristic()
                DA_0.append(best_candidateD)
                if k not in resultsD: resultsD[k] = 0
                resultsD[k] += G.sigma(DA_0)

                # random selection of initial active set
                randNode = random.choice(vertices)
                vertices.remove(randNode)
                RA_0.append(randNode)
                if k not in resultsR: resultsR[k] = 0
                resultsR[k] += G.sigma(list(RA_0))


        # average the results
        for k in resultsD: 
            resultsD[k] = resultsD[k] // iterations
            resultsG[k] = resultsG[k] // iterations
            resultsR[k] = resultsR[k] // iterations

        # save results to file
        with open("./results2/" + str(p) + ".txt", 'w') as f:
            f.write(json.dumps(resultsG))
            f.write('\n')
            f.write(json.dumps(resultsD))
            f.write('\n')
            f.write(json.dumps(resultsR))
            f.write('\n')


if __name__ == '__main__':
    main() 
    