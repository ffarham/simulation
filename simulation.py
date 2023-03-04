import random
import numpy as np
import logging

"""
    HRCGraph - class for Homophilic Relaxed Caveman Graph
"""
class HRCGraph:

    """ constructor - generate a graph
        :param l - number of cliques
        :param s - size of each clique
        :param p - rewire probability
        :param h - homophily factor
    """
    def __init__(self, l, s, p, h):

        self.n = l * s
        self.m = s * (s-1) * s * l
        self.vertices = [i for i in range(1, (l*s) + 1)]
        self.edges = dict()
        self.thresholds = dict()
        self.beliefs = dict()
        self.mcgs = []

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

        
        # initialise the thresholds + beliefs
        for i in range(1, len(self.vertices) + 1):
            self.thresholds[i] = random.uniform(0,1)
            self.beliefs[i] = random.uniform(0,1)
        
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
                
                r = None
                # if both nodes u and v have the same belief
                if (self.sameColour(u, v)): r = p * h
                else: r = p * (1 - h)

                # rewire with probability r
                if random.uniform(0,1) <= r:
                    if x in self.edges[u]: self.edges[u][x] += self.edges[u][v]
                    else: self.edges[u][x] = self.edges[u][v]

                    del self.edges[u][v]
    

    """ sameColour - determine if 2 nodes have the same opinion
        - that is their belief either both exceeds their respective thresholds or they dont
    """
    def sameColour(self, u, v):
        return (self.beliefs[u] >= self.thresholds[u] and self.beliefs[v] >= self.thresholds[v]) or (self.beliefs[u] < self.thresholds[u] and self.beliefs[v] < self.thresholds[v])


    """ transpose - method to reverse the edges in the graph
        - creates a new dictionary of edges
    """
    def transpose(self):
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
    def dfs(self, u, visited, currentMCG):
        visited[u] = True
        currentMCG.append(u)
        for v in self.edges[u]:
            if not visited[v]: self.dfs(v, visited, currentMCG)
    

    """ dfsNumbering - number the nodes in the network by adding them to the stack
        - used to get the order of nodes to perform the DFS that determines minimal closed groups
    """
    def dfsNumbering(self, u, visited, stack):
        visited[u] = True
        for v in self.edges[u]:
            if not visited[v]: self.dfsNumbering(v, visited, stack)
        stack = stack.append(u)
        return


    """ mcg - determines the minimal closed groups in the graph with no outgoing edges
        - Kosaraju's algorithm: first perform DFS numbering on nodes and use that numbering as an order of nodes to traverse using DFS, each traversal traverses a minimal closed group
    """
    def mcg(self):
        logging.info("determining minimal closed groups in the network")
        stack = []
        visited = [False] * (len(self.vertices) + 1)
        visited[0] = None

        for u in range(1, len(self.vertices) + 1):
            if not visited[u]: self.dfsNumbering(u, visited, stack)
        
        self.transpose()

        visited = [False] * (len(self.vertices) + 1)
        while len(stack) > 0:
            currentMCG = []
            u = stack.pop()
            if not visited[u]:
                self.dfs(u, visited, currentMCG)
                self.mcgs.append(currentMCG)

        self.transpose()

        return 

        
    """ contraction - contract the graph on terminating group (minimal closed groups with no outgoing edges)
        :param tg - termianting group to contract
        :param consensus - the consensus belief reached by the terminating group
        :param R - nodes outside of all terminating groups 
    """
    def contract(self, tg, consensus, R):
        logging.info("contracting terminating groups")
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
        for u in R:
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

        # update the vector of beliefs 
        self.beliefs[sv] = consensus
        for node in tg:
            del self.beliefs[node]

        # remove each node in the minimal closed group
        while len(tg) > 0:
            node = tg.pop(0)
            self.vertices.remove(node)

        # add supervertex to the vertices
        self.vertices.append(sv)

    
    """ calculateConsensus - function to calculate the consensus belief reached by a minimal closed group
        :param mcg - minimal closed group to calculate the consensus belief of
        :param A_0 - the initial active set of nodes
        :return the consensus belief 
    """
    def calculateConsensus(self, mcg, A_0):
        # activate the nodes in the terminating groups
        A_0 = set(A_0)
        for node in mcg:
            if node in A_0: self.beliefs[node] = 1
        
        # calculate the consesnus reached by each terminating group
        mcg.sort()
        beliefVector = np.array([self.beliefs[i] for i in mcg])
        T = []
        for i in range(len(mcg)):
            temp = [0] * len(mcg)
            T.append(temp)

        for i in range(len(mcg)):
            for j in range(len(mcg)):
                if mcg[j] in self.edges[mcg[i]]: T[i][j] = self.edges[mcg[i]][mcg[j]]
        
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
        
        # calculate the consensus reached by the minimal closed group
        consensus = np.dot(beliefVector, leftEigenVector)

        return consensus


    """ sigma - influence function
        :param A_0 - initial set of active nodes
        :return the total number of active nodes at convergence
    """
    def sigma(self, A_0):
        logging.info("activation function called")
        output = 0

        # calculate minimal closed groups if required
        if len(self.mcgs) == 0 and len(self.vertices) > 0: self.mcg()    

        # determine the terminating groups with no outgoing edges
        TG = []
        for group in self.mcgs:
            flag = True
            for u in group:
                # determine if there is an edge going outside the minimal closed group
                for v in self.edges[u]:
                    if v not in group:
                        flag = False
                        break
                if not flag: break

            if flag: TG.append(group)

        # determine the nodes not in a terminating group
        temp = set([item for group in TG for item in group])
        R = [i for i in range(1, len(self.vertices) + 1) if i not in temp]

        
        for group in TG:
            # calculate consensus reached by each terminating group + keep track of activated nodes withing
            consensus = self.calculateConsensus(group, A_0)

            # keep track of any activated nodes
            for node in group:
                if consensus >= self.thresholds[node]: output += 1

            # contract the graph on the minimal closed group
            self.contract(group, consensus, R)

        # if all nodes belong to some minimal closed group
        if len(R) == 0: return output

        # still need to calculate convergent belief of each node in R -> solve system of linear equations
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
        
        b = np.array(list(self.beliefs.values()))
        
        # calculate the vector of convergent beliefs in G'
        c = A.dot(b.T)

        # count the number of vertices in R that get activated
        for i, r in enumerate(R):
            if c[i] >= self.thresholds[r]: output += 1
        
        return output

    
    def getVertices(self):
        return self.vertices

    def getEdges(self):
        return self.edges
    
    def getThresholds(self):
        return self.thresholds.values()
    
    def getBeliefs(self):
        return self.beliefs.values()
    
    def getMCGs(self):
        if len(self.mcgs) == 0 and len(self.vertices) > 0:
            self.mcg()
        return self.mcgs

    def __str__(self):
        return "Vertices: " + str(self.vertices) + "\nEdges: " + str(self.edges)


def main():
    logging.basicConfig(level=logging.INFO)
    
    # TODO: maybe test the results for different rewiring probability and homophily factor
    G = HRCGraph(10, 10, 1, 0)
    sig = G.sigma([])
    print("\nTotal activations: ", sig)
    return

if __name__ == '__main__':
    main() 