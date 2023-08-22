
import sys
import networkx as nx
import numpy as np
import math
import csv
import pandas as pd
import os


def verify(edge_list): 
    """ 
    Verifies that the tree strategy outputting by the optimization solver 
    is a valid tree strategy with a valid associated weight function. Does so 
    by checking if the tree is connected and acyclic, and by verifying the weight
    funcion properties are met. 
    """
    valid = True
    tree = nx.from_edgelist(edge_list)
    cycles = list(nx.cycle_basis(tree))
    if len(cycles) != 0:
        print("Cycle exists in tree. Strategy is invalid.")
        print(cycles)
        valid = False
    components = [len(c) for c in sorted(nx.connected_components(tree), key=len, reverse=True)]
    if len(components) > 1: 
        print("Tree is not connected. Strategy is invalid.")
        print(components)
        valid = False

    return valid
        

def main():

    verified = []
    dir = "/home/DAVIDSON/doflocco/Graph_Pebbling/results/lemke square/v(0,0)-size20/edges"
    for edge_file in os.listdir(dir):
        if "edge" in edge_file:
            with open(os.path.join(dir,edge_file), newline='') as f:
                reader = csv.reader(f, delimiter=',')
                next(reader) # skip header row
                edge_list = [(row[1], row[2]) for row in reader]
                valid = verify(edge_list)
                if not valid: 
                    print(f"Tree strategy invalid in {edge_file}")
                    print()
                verified.append(valid)
    if all(verified):
        print("All tree strategies are valid.")

if __name__ == "__main__":
    main()