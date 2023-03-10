import os
import json
import matplotlib.pyplot as plt

def main():
    ps = [0, 0.2, 0.4, 0.6, 0.8, 1]
    for p in ps:
        file_path = "./results2/" + str(p)+".txt"
        if os.path.exists(file_path):
            with open("results2/" + str(p)+".txt", 'r') as f:

                # plot the results
                data = json.loads(f.readline()) 
                x = [ int(xd) for xd in data.keys() ]
                y = [ int(yd) for yd in data.values() ]
                labelValue = "Greedy"
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

                plt.title("Re-wiring probability p = "+str(p))
                plt.legend(loc="upper left")
                plt.xlabel("Size of Initial Active Set")
                plt.ylabel("Size of Resulting Active Set")

                plt.show() 

if __name__ == "__main__":
    main()