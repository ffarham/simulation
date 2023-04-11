import json
import matplotlib.pyplot as plt
import logging

def main():
    logging.basicConfig(level=logging.INFO)
    
    # NOTE: assign the clique to plot features of
    clique_size = "5"
    
    with open('./results/structure_T.txt', 'r') as f:
        data = json.loads(f.readline())[clique_size] 
        x = [ float(xd) for xd in data.keys() ]
        y = [ float(yd) for yd in data.values() ]
        labelValue = "Terminals"
        plt.plot(x,y, label=labelValue)
    
    with open('./results/structure_nonT.txt', 'r') as f:
        data = json.loads(f.readline())[clique_size] 
        x = [ float(xd) for xd in data.keys() ]
        y = [ float(yd) for yd in data.values() ]
        labelValue = "Non-Terminals"
        plt.plot(x,y, label=labelValue)
        
     

    # plt.title("Re-wiring probability p = "+str(p))
    plt.legend(loc="upper right")
    plt.xlabel("Re-wiring Probability")
    plt.ylabel("No. of Terminals or Non-Terminals")

    plt.show() 

if __name__ == "__main__":
    main()