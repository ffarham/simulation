import random
import numpy as np
import logging

"""
    HRCGraph - class for Homophilic Relaxed Caveman Graph
"""
class Graph:

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
        self.convergent_beliefs_HC = dict()
        self.convergent_beliefs_AHC = dict()
        self.convergent_beliefs_S = dict()
        self.split = False
        self.reShape = False
        self.convergentT = None
        self.mapSVtoVs = dict()
        self.mapSVtoS= dict()
        self.mapVtoSV = dict()
        self.terminating_groups = []
        self.terminals = set()
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
            # NOTE: initialising the initial vector of beliefs to be a zero vector
            self.initBeliefs[node] = 0 
        
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
        
        # determine the in-degree of each node -> used for degree heuristic selection
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

        T = np.array([[self.edges[tg[i]].get(tg[j], 0) for j in range(len(tg))] for i in range(len(tg))])
        T /= np.sum(T, axis=1, keepdims=True)
        np.fill_diagonal(T, np.diag(T) + (1 - np.sum(T, axis=1)))

        # # get eigen vector and corresponding eigen value
        eigenValues, eigenVectors = np.linalg.eig(T.T)

        # # find eigen vector corresponding to eigen value 1
        index = np.where(np.isclose(eigenValues, 1))[0][0]
        leftEigenVector = eigenVectors[:, index].real

        # # Normalize the left eigenvector
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

            self.terminals.add(sv)
            self.mapSVtoVs[sv] = group
            self.mapSVtoS[sv] = s
            self.convergent_beliefs_HC[sv] = 0
            self.convergent_beliefs_AHC[sv] = 0
            self.convergent_beliefs_S[sv] = 0
            for node in group:
                self.mapVtoSV[node] = sv

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
        # TODO: test the cumulation on some fixed num. of iterations for all instances of A
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
    """
    def hill_climbing(self, A_0, candidates):

        # keep track of best candidates that purely increase the resulting size of active set the most
        best_candidates = []
        best_value = -1
        mapCtoCB = dict()

        # iterate through all candidates to determine the best ones
        for c in candidates:
            total_activations = 0
            temp_activation = A_0 + [c]

            # calculate the consensus reached by the minimal closed group that c is part of
            sv = self.mapVtoSV[c]
            s = self.mapSVtoS[sv]
            nodes = self.mapSVtoVs[sv]
            p_sv = [ 1 if x in temp_activation else self.initBeliefs[x] for x in nodes ]
            current_consensus = self.convergent_beliefs_HC[sv]
            new_consensus = np.dot(s, p_sv)
            self.convergent_beliefs_HC[sv] = new_consensus
            mapCtoCB[c] = (sv, new_consensus)

            # calculate the convergent belief of each node
            A = self.convergentT
            p = np.array([ self.convergent_beliefs_HC[x] if x in self.terminals else self.initBeliefs[x] for x in self.vertices ])
            convergent_p = A.dot(p.T)
            self.convergent_beliefs_HC[sv] = current_consensus

            # calcualte the number of nodes activated
            for i, v in enumerate(self.vertices):
                if v in self.terminals:
                    for node in self.mapSVtoVs[v]:
                        if convergent_p[i] >= self.thresholds[node]: total_activations += 1
                else:
                    if convergent_p[i] >= self.thresholds[v]: total_activations += 1


            # keep track of the best candidate
            if total_activations > best_value: best_candidates, best_value = [(c, total_activations)], total_activations
            elif total_activations == best_value: best_candidates.append((c, total_activations))

        # update the map of convergent beliefs 
        selected = random.choice(best_candidates)
        selected_SV, consensus = mapCtoCB[selected[0]]
        self.convergent_beliefs_HC[selected_SV] = consensus
        return selected


    def hill_climbing_augmented(self, A_0, candidates):

        # keep track of best candidates that increase the resulting size of active set and also increase the total belief in the network the most
        best_candidates = []
        best_value, best_convergent_belief_value = -1, -1
        mapCtoCB = dict()

        # iterate through all candidates to determine the best one
        for c in candidates:
            total_activations = 0
            temp_activation = A_0 + [c]

            # calculate the consensus reached by the minimal closed group that c is part of
            sv = self.mapVtoSV[c]
            s = self.mapSVtoS[sv]
            nodes = self.mapSVtoVs[sv]
            p_sv = [ 1 if x in temp_activation else self.initBeliefs[x] for x in nodes ]
            current_consensus = self.convergent_beliefs_AHC[sv]
            new_consensus = np.dot(s, p_sv)
            self.convergent_beliefs_AHC[sv] = new_consensus
            mapCtoCB[c] = (sv, new_consensus)

            # calculate the convergent belief of each node
            A = self.convergentT
            p = np.array([ self.convergent_beliefs_AHC[x] if x in self.terminals else self.initBeliefs[x] for x in self.vertices ])
            convergent_p = A.dot(p.T)
            self.convergent_beliefs_AHC[sv] = current_consensus
            value = sum(convergent_p)

            # calcualte the number of nodes activated
            for i, v in enumerate(self.vertices):
                if v in self.terminals:
                    for node in self.mapSVtoVs[v]:
                        if convergent_p[i] >= self.thresholds[node]: total_activations += 1
                else:
                    if convergent_p[i] >= self.thresholds[v]: total_activations += 1
            
            # keep track of the best candidate
            if total_activations > best_value: best_candidates, best_value, best_convergent_belief_value = [(c, total_activations)], total_activations, value
            elif total_activations == best_value:
                if value > best_convergent_belief_value: best_candidates, best_value, best_convergent_belief_value = [(c, total_activations)], total_activations, value
                elif value == best_convergent_belief_value: best_candidates.append((c, total_activations))

        # update the map of convergent beliefs 
        selected = random.choice(best_candidates)
        selected_SV, consensus = mapCtoCB[selected[0]]
        self.convergent_beliefs_AHC[selected_SV] = consensus
        return selected


    """ sigma - influence function
        :param A_0 - initial set of active nodes
        :return the total number of active nodes at convergence
    """
    def sigma(self, A_0):
        
        if not self.reShape: self.reShapeGraph()

        total_activations = 0

        # determine convergent belief of all supervertices
        for sv in self.terminals:
            s = self.mapSVtoS[sv]
            nodes = self.mapSVtoVs[sv]
            # create belief vector of desired nodes + any activations
            p = [self.initBeliefs[x] if x not in A_0 else 1 for x in nodes]

            consensus = np.dot(s, p)
            self.convergent_beliefs_S[sv] = consensus
            # determine the nodes activated in the respective terminating group
            for node in nodes:
                if consensus >= self.thresholds[node]: total_activations += 1
        
        # if all nodes belong to some minimal closed group
        if len(self.non_terminals) == 0: return total_activations

        # still need to calculate convergent belief of each non terminals
        A = self.convergentT        
        b = np.array([self.initBeliefs[x] if x in self.initBeliefs else self.convergent_beliefs_S[x] for x in self.vertices])
        
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
    def reset(self):
        for node in self.thresholds:
            self.thresholds[node] = random.uniform(0,1)
        for node in self.convergent_beliefs_HC:
            self.convergent_beliefs_HC[node] = 0
            self.convergent_beliefs_AHC[node] = 0

    """ getCandidates - get the candidates to activate
        :return the nodes in the terminating groups
    """
    def getCandidates(self):
        if not self.split: self.__splitGraph()
        
        candidates = []
        for group in self.terminating_groups:
            candidates += group
        
        return candidates


    """ get_degree_ordered_nodes - return the nodes ordered in largest incoming degree
    """
    def get_degree_ordered_nodes(self):
        return self.degree_heuristics


    """ getVertices - return all vertices in the initial network
    """
    def getVertices(self):
        return self.initial_vertices
    
    def network_structure(self):
        return [self.terminals , self.non_terminals, [ (u,v) for u in self.edges for v in self.edges[u] if self.edges[u][v] > 0] ]


    def __str__(self):
        return "Vertices: " + str(self.vertices) + "\nEdges: " + str(self.edges)

