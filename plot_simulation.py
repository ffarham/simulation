import json
import matplotlib.pyplot as plt
import logging

def main():
    logging.basicConfig(level=logging.INFO)

    p = 0.05
    
    with open('./results/simulation.txt', 'r') as f:
        
        # plot the results
        data = json.loads(f.readline()) 
        x = [ int(xd) for xd in data.keys() ]
        y = [ int(yd) for yd in data.values() ]
        labelValue = "Greedy"
        plt.plot(x,y, label=labelValue)
        
        data = json.loads(f.readline()) 
        x = [ int(xd) for xd in data.keys() ]
        y = [ int(yd) for yd in data.values() ]
        labelValue = "Augmented Greedy"
        plt.plot(x,y, label=labelValue)

        data = json.loads(f.readline())   
        x = [ int(xd) for xd in data.keys() ]
        y = [ int(yd) for yd in data.values() ]
        labelValue = "Degree"
        plt.plot(x,y,label=labelValue)

        data = json.loads(f.readline()) 
        x = [ int(xd) for xd in data.keys() ]
        y = [ int(yd) for yd in data.values() ]
        labelValue = "Random"
        plt.plot(x,y,label=labelValue)

        # plt.title("Re-wiring probability p = "+str(p))
        plt.legend(loc="upper left")
        plt.xlabel("Size of Initial Active Set")
        plt.ylabel("Size of Resulting Active Set")

        plt.show() 

if __name__ == "__main__":
    main()