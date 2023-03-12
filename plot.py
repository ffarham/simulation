import os
import json
import matplotlib.pyplot as plt
import logging

def main():
    logging.basicConfig(level=logging.INFO)

    # NOTE: set the re-wiring probabilities to load
    ps = [0, 0.2, 0.4, 0.6, 0.8, 1]

    for p in ps:
        # NOTE: define directory to load results from
        dir_path = "./results2/"
        assert os.path.exists(dir_path), "Directory to load results from is not defined"

        file_path = dir_path + str(p)+".txt"
        if os.path.exists(file_path):
            # fetch simulated data from the file
            with open(file_path, 'r') as f:

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
        else:
            logging.info("file " + file_path + " not found")

if __name__ == "__main__":
    main()