import random

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
                N[m] = round(N[m]/total,3)
            
            self.edges[i] = N
        
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

    """ MCG - determine the minimal closed groups in the graph with no outgoing edges
    """
    def MCG(self):
        return 
        
    """ sigma - influence function
        :param A - initial set of active nodes
        :return the total number of active nodes at convergence
    """
    def sigma(self, A):
        


        return 0
    
    def getVertices(self):
        return self.vertices

    def getEdges(self):
        return self.edges
    
    def getThresholds(self):
        return self.thresholds.values()
    
    def getBeliefs(self):
        return self.beliefs.values()

    def __str__(self):
        return "Vertices: " + str(self.vertices) + "\nEdges: " + str(self.edges)


def main():
    # TODO: maybe test the results for different rewiring probability and homophily factor
    G = HRCGraph(3, 2, 0.5, 0.5)
    print(G)
    G.transpose()
    print(G)

    return

if __name__ == '__main__':
    main()