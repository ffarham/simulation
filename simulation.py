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
        for i in range(1, self.n+1):
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

            # ensure they sum to 1

        
        # initialise the thresholds + beliefs
        for i in range(1, self.n + 1):
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
                
                x = random.randint(1, self.n)
                while(x == v):
                    x = random.randint(1, self.n)
                
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
    """
    def transpose(self):
        visited = set()
        for i in self.edges:
            queue =  list(self.edges[i].keys())
            while(len(queue) > 0):
                j = queue.pop(0)

                if j in visited: continue

                w_ij = self.edges[i][j]
                w_ji = self.edges[j][i] if i in self.edges[j] else None 
                
                self.edges[j][i] = w_ij

                if w_ji is not None:
                    self.edges[i][j] = w_ji
                else:
                    del self.edges[i][j]

            visited.add(i)
        return
    
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
        visited = [False] * (self.n + 1)
        visited[0] = None

        for u in range(1, self.n +1):
            if not visited[u]:
                self.numbering(u, visited, stack)
        
        self.transpose()

        visited = [False] * (self.n + 1)
        while len(stack) > 0:
            currentMCG = []
            u = stack.pop(0)
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
        :param A - initial set of active nodes
        :return the total number of active nodes at convergence
    """
    def sigma(self, A):
        # calculate minimal closed groups if required
        if len(self.mcgs) == 0 and self.n > 0:
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
        R = [i for i in range(1, self.n + 1) if i not in temp]

        # activate the nodes in the terminating groups
        A = set(A)
        for group in TG:
            for node in group:
                if node in A:
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
            
            T = np.array(T)
            eigvals, eigvecs = np.linalg.eig(T.T)

            idx = np.where(np.isclose(eigvals, 1))[0][0]
            left_eigvec = eigvecs[:, idx].real

            # Normalize the left eigenvector
            left_eigvec = left_eigvec / np.sum(left_eigvec)

            print("Left eigenvector:", left_eigvec)
            
            consensus = np.dot(beliefVector, left_eigvec)
            print("beliefs: ", beliefVector)
            print("consensus: ", consensus)

        # solve system of linear equations

        return 0
    
    def getVertices(self):
        return self.vertices

    def getEdges(self):
        return self.edges
    
    def getThresholds(self):
        return self.thresholds.values()
    
    def getBeliefs(self):
        return self.beliefs.values()
    
    def getMCGs(self):
        if len(self.mcgs) == 0 and self.n > 0:
            self.mcg()
        return self.mcgs

    def __str__(self):
        return "Vertices: " + str(self.vertices) + "\nEdges: " + str(self.edges)


def main():
    # TODO: maybe test the results for different rewiring probability and homophily factor
    G = HRCGraph(3, 2, 0.5, 0.5)
    print(G)
    print(G.getMCGs())
    G.sigma([])
    return

if __name__ == '__main__':
    main()