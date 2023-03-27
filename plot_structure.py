import json
import matplotlib.pyplot as plt
import numpy as np

def main():
    data = {}
    with open("./results/structure.txt", 'r') as f:
        data = json.loads(f.readline())
    
    fig = plt.figure(figsize = (12,10))
    ax = plt.axes(projection='3d')

    Z = np.array([[ data[clique_size][p] for p in data[clique_size] ] for clique_size in data ])
    x = [ int(xd) for xd in data.keys() ]
    y = [ float(yd) for yd in data[str(x[0])].keys() ]
    X, Y = np.meshgrid(y, x)

    ax.plot_surface(X,Y,Z, cmap=plt.cm.cividis)
    ax.set_xlabel("Re-wiring Probability p", labelpad=20)
    ax.set_ylabel("Clique Size", labelpad=20)
    ax.set_zlabel("No. of Terminating Groups", labelpad=20)
    plt.yticks([ y[0] for y in Y])

    plt.title("")

    plt.show()


if __name__ == "__main__":
    main()
