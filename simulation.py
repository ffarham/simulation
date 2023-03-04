import random
import numpy as np

"""
    HRCGraph - class for Homophilic Relaxed Caveman Graph
"""
class HRCGraph:

    """ constructor
        :param l - number of cliques
        :param k - size of each clique
        :param p - rewire probability
        :param h - homophily factor
    """
    def __init__(self, l, k, p, h):

        self.n = l * k
        self.m = k * (k-1) * k * l
        self.vertices = [i for i in range(1, (l*k) + 1)]
        self.edges = dict()
        self.thresholds = dict()
        self.beliefs = dict()
        self.mcgs = []

        # create edge set for clique with self loops
        for i in range(1, len(self.vertices)+1):
            a = ((i-1) // k) * k if i % k == 0 else (i // k) * k
            b = a + k
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
                if ((self.beliefs[u] >= self.thresholds[u] and self.beliefs[v] >= self.thresholds[v]) or (self.beliefs[u] < self.thresholds[u] and self.beliefs[v] < self.thresholds[v])):
                    r = p * h
                else:
                    r = p * (1 - h)

                # rewire with probability r
                if random.uniform(0,1) <= r:
                    if x in self.edges[u]:
                        self.edges[u][x] += self.edges[u][v]
                    else:
                        self.edges[u][x] = self.edges[u][v]

                    del self.edges[u][v]
    

    """ transpose - method to transpose the edges in the graph
        - creates a new dictionary of edges
    """
    def transpose(self):
        newE = {}
        for i in self.edges:
            for j in self.edges[i]:
                
                if j not in newE:
                    newE[j] = {}

                newE[j][i] = self.edges[i][j]

        self.edges = newE

    
    """ dfs
    """
    def dfs(self, u, visited, currentMCG):
        visited[u] = True
        currentMCG.append(u)
        for v in self.edges[u]:
            if not visited[v]:
                self.dfs(v, visited, currentMCG)
    

    """ numbering
    """
    def numbering(self, u, visited, stack):
        visited[u] = True
        for v in self.edges[u]:
            if not visited[v]:
                self.numbering(v, visited, stack)
        stack = stack.append(u)
        return


    """ mcg - determine the minimal closed groups in the graph with no outgoing edges
    """
    def mcg(self):
        stack = []
        visited = [False] * (len(self.vertices) + 1)
        visited[0] = None

        for u in range(1, len(self.vertices) + 1):
            if not visited[u]:
                self.numbering(u, visited, stack)
        
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

        
    """ contraction - contract the graph on minimal closed groups with no outgoing edges
    """
    def contraction(self):
        return


    """ sigma - influence function
        :param A_0 - initial set of active nodes
        :return the total number of active nodes at convergence
    """
    def sigma(self, A_0):
        output = 0

        # calculate minimal closed groups if required
        if len(self.mcgs) == 0 and len(self.vertices) > 0:
            self.mcg()    

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
                if not flag:
                    break
            if flag:
                TG.append(group)

        # determine the nodes not in a terminating group
        temp = set([item for group in TG for item in group])
        R = [i for i in range(1, len(self.vertices) + 1) if i not in temp]

        # activate the nodes in the terminating groups
        A_0 = set(A_0)
        for group in TG:
            for node in group:
                if node in A_0:
                    self.beliefs[node] = 1
        
        # calculate the consesnus reached by each terminating group
        for group in TG:
            group.sort()
            beliefVector = np.array([self.beliefs[i] for i in group])
            T = []
            for i in range(len(group)):
                temp = [0] * len(group)
                T.append(temp)

            for i in range(len(group)):
                for j in range(len(group)):
                    if group[j] in self.edges[group[i]]:
                        T[i][j] = self.edges[group[i]][group[j]]
            
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
            
            # keep track of any activated nodes
            for node in group:
                if consensus >= self.thresholds[node]:
                    output += 1

            # contract the minimal closed group
            # remove the edges within the minimal closed group
            for u in group:
                for v in group:
                    if v in self.edges[u]:
                        del self.edges[u][v]
            
            # create supervertex
            sv = self.n + 1
            self.n += 1

            # add self-loop to the supervertex
            self.edges[sv] = dict()
            self.edges[sv][sv] = 1

            # rewire any edge connected to minimal closed group to the supervertex
            for u in R:
                for v in group:
                    # NOTE: by property of TG, there are no outgoing edges in the group
                    # incoming edges to minimal closed group: (u,v)
                    if v in self.edges[u]:
                        # check if an edge to supervertex already exists
                        if sv not in self.edges[u]:
                            self.edges[u][sv] = 0
                        
                        self.edges[u][sv] += self.edges[u][v]

                        del self.edges[u][v]          

            # clean up the edge set
            cuNodes = list(self.edges.keys())
            for node in cuNodes:
                if not self.edges[node]:
                    del self.edges[node]

            # update the vector of beliefs 
            self.beliefs[sv] = consensus
            for node in group:
                del self.beliefs[node]

            # remove each node in the group
            while len(group) > 0:
                node = group.pop(0)
                self.vertices.remove(node)

            # add supervertex to the vertices
            self.vertices.append(sv)

        if len(R) == 0: 
            return output

        # still need to calculate convergent belief of each node in R -> solve system of linear equations
        # define a matrix of co-efficients
        A = []
        for i in range(len(self.vertices)):
            temp = [0] * len(self.vertices)
            A.append(temp)
        
        for i, u in enumerate(self.vertices):
            for j, v in enumerate(self.vertices):
                if v in self.edges[u]:
                    A[i][j] = self.edges[u][v]
        A = np.array(A)

        # converge A to some epsilon accuracy
        counter, epsilon = 0, 0.0001
        prev = A
        A = np.linalg.matrix_power(A, 2)
        while (np.sum(np.abs(A - prev)) > epsilon):
            if counter > 100: 
                print("WARNING: hit counter limit")
                break
            prev = A
            A = np.linalg.matrix_power(A, 2)
            counter += 1
        
        b = np.array(list(self.beliefs.values()))
        
        # calculate the vector of convergent beliefs in G'
        c = A.dot(b.T)

        # count the number of vertices in R that get activated
        for i, r in enumerate(R):
            if c[i] >= self.thresholds[r]:
                output += 1
        
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
    # TODO: maybe test the results for different rewiring probability and homophily factor
    G = HRCGraph(3, 2, 0.5, 0.5)
    sig = G.sigma([])

    print("Total activations: ", sig)
    return

if __name__ == '__main__':
    main()