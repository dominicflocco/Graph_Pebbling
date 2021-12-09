import networkx as nx
import numpy as np
import math
import pandas as pd
import scipy.special

from Optimizer import Optimizer 
from TreeStrategy import TreeStrategy 
from TreeStrategy import NonTreeStrategy
from PebblingGraph import PebblingGraph 

def countElementaryCirtuits(n):
   
    count = 0 
    for i in range(1, n):
        coef = scipy.special.binom(n, n-i+1)
        count += coef*math.factorial(n-i)
    return count

def boundLemkeSquare(lemkeSquare): 

    nodes = lemkeSquare.nodes
    bounds = {}
    visited = set()
    for r in nodes:
        if r not in visited: 
            lemkeSquare.root = r
            
            optimizer = Optimizer(solver = "Gurobi", 
                                threads = 32, 
                                logFile = "lemke_square_v-" + str(r),
                                lpFile = "lemke_square_v-" + str(r),
                                objGap=False,
                                timeLimit=False, 
                                cert=True)
            
            res = optimizer.generateTreeStrategiesSym(lemkeSquare, 20, 10)
            bounds[r] = res[1]
        visited.add(r)
        r_rev = "("+ r[4]+ ", " + r[1] +")"
        visited.add(r_rev)
    
    pd.DataFrame.from_dict(bounds, orient="index").to_csv("lemke_square_bounds.csv")

def runNontreeOpt(graph): 

    optimizer = Optimizer(solver = "Gurobi", 
                                threads = 32, 
                                logFile = "lemke_maxUnsolvable.log",
                                lpFile = "lemke_maxUnsolvable.lp",
                                objGap=False,
                                timeLimit = 3600, 
                                cert=True)
    
    lolipops = graph.generateLolipops()  
    strategies = {}
    i = 0
    for pop in lolipops: 
        strategies[i] = NonTreeStrategy(graph, list(pop.edges))
        print(strategies[i].edges)
        print(strategies[i].weights)
        i += 1
    res = optimizer.maxUnsolvable(strategies, graph, graph.root)

    #print(res)
def runHybridOpt(graph):

    optimizer = Optimizer(solver = "Gurobi", 
                                threads = 32, 
                                logFile = "lemke_hybrid.log",
                                lpFile = "lemke_hybrid.lp",
                                objGap=False,
                                timeLimit = False, 
                                cert=True)
    
    lollipops = graph.generateLolipops(cycleLen=4)

    strategies = []
    for pop in lollipops: 
        strategies.append(NonTreeStrategy(graph, list(pop.edges)))
        
    optimizer.hybridOpt(graph, strategies, 5, 5)
def runMaxUnsolvable(graph): 
    optimizer = Optimizer(solver = "Gurobi", 
                                threads = 32, 
                                logFile = "3cube-maxunsolvable.log",
                                lpFile = "3cube-maxunsolvable.lp",
                                objGap=False,
                                timeLimit = False, 
                                cert=True)
    
    strategies = graph.generateLollipops(cycleLen=4)
    bound = optimizer.maxUnsolvable(strategies, graph)
    print("root = " + str(graph.root) + " -> " + str(bound))

def main(): 
    # Lemke 
    lemke_edges = [(0, 1), (0, 2), (1, 3), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6), (3, 7), (4, 7), (5, 7), (6, 7)]
    r = 2
    lemke = PebblingGraph(lemke_edges, r)
    #runMaxUnsolvable(lemke)
    lemke_nodes = lemke.nodes
    for r in lemke_nodes:
        lemke = PebblingGraph(lemke_edges, r)
        #runMaxUnsolvable(lemke)
    #strategies = lemke.generateLollipops()
    
    # # Lemke Square 
    # lemkeSquare = nx.cartesian_product(lemke.graph, lemke.graph)
    # lemkeSquare = list(lemkeSquare.edges())
    # r = (0,0)
    # #print(countElementaryCirtuits(64))
    # lemkeSquare = PebblingGraph(lemkeSquare, r)
    
    # # Bruhat 
    # # bruhat = [(0,1), (0,2), (0,4), (1,2), (1,7), (2,3), (2,20), (3,23), (4, 5), (4, 8), (5,6), (5,9), (6,7), (6,10), (7,11),
    # #           (8,9), (8,16), (9,12), (10,11), (10,13), (11, 19), (12, 13), (12,14), (13,15), (14, 15), (14,17), (15, 18),
    # #           (16, 17), (16,20), (17, 21), (18,19), (18,22), (19,23), (20, 21), (21,22), (22,23)]
    # # r = 1
    # # bruhat = PebblingGraph(bruhat, r)

    # # n-cube
    # n = 3
    # hypercube = [(0,1), (0,2), (0,4), (1,3), (1,5), (2,3), (2,6), (3,7), (4,5), (4,6), (5,7), (6,7)]
    hypercube = [(0,1), (0,2), (0,4), (1,3), (1,5), (2,3), (2,6), (3,7), (4,5), (4,6), (5,7), (6,7)]
    
    # # hypercube = nx.hypercube_graph(n)
    # # hypercube = hypercube.edges
    r = 1
    hypercube = PebblingGraph(hypercube, r)
    
    strategies = hypercube.generateLollipops(cycleLen=4)
    for s in strategies:
        print(s.weights)
    runMaxUnsolvable(hypercube)
    # # print(len(hypercube.edges))
    # # print(len(hypercube.nodes))
    # runNontreeOpt(hypercube)

    # Petersen 
    # petersen = [(0,1), (1,2), (2,3), (3,4), (4,0), (0,6), (1,7), (2,8), (3,9), (4,5), (5,7), (5,8), (6,9), (6,8), (7,9)]
    # r = 5
    # petersen = PebblingGraph(petersen, r)

    # Lolipop 
    # lolipop = [(0,1), (1,2), (2,3), (3,4), (4,1)]
    # lolipop = PebblingGraph(lolipop, r)

    #lolipops = lemke.generateLolipops(cycleLen=4) 
    # strategies = {}
    # i = 0

    # for pop in lolipops:
    #     #print(pop)
    #     strategies[i] = NonTreeStrategy(lemke, pop.edges)
    #     print(strategies[i].weights)
    #     i+=1
    #boundLemkeSquare(lemkeSquare)

    # optimizer = Optimizer(solver = "Gurobi", 
    #                             threads = 32, 
    #                             logFile = "lemke_square_v-" + str(r),
    #                             lpFile = "lemke_square_v-" + str(r),
    #                             objGap=False,
    #                             timeLimit=False, 
    #                             cert=True)
    
    # res = optimizer.generateTreeStrategiesSym(lemkeSquare, 20, 10)

if __name__ == "__main__":
    main()